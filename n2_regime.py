#!/usr/bin/env python3
"""
n2_regime.py  v0.1.0   N2 failure-regime diagnostic (specification N2_REGIME_DIAGNOSTIC_SPEC_v0.3.md).

Exploratory, descriptive, one run, on frozen outputs only. Reads:
  * n2_blockfolds.csv            (N2 replay; 4,760 block-folds)                    sha256 b7934523...
  * gate3a_blocks.jsonl          (Gate 3A per-block output; G1 fold agreement)     sha256 6decb9fb...
  * n0_records.csv               (N0 regularization parity; O-only kernel edge)    sha256 6d3b0341...
  * paper4_design_models/*/design_model_manifest.json  (n_design, block feature names, design statistics)
It reads no canonical matrix, no design row, no evaluation row, fits no selector and no predictor.

Usage
  python3 n2_regime.py selftest
  python3 n2_regime.py run --blockfolds PATH --gate3a-blocks PATH --n0-records PATH --manifests DIR \
                           --spec N2_REGIME_DIAGNOSTIC_SPEC_v0.3.md --out n2_regime_out
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, sys, time
from pathlib import Path
import numpy as np

VERSION = "0.1.1"   # v0.1.0 crashed in A2 (LinAlgError, SVD did not converge) on the real run before any output; see §8 of the spec
SEED = 20261103
B_BOOT = 10_000
ALPHA_TWO_SIDED = 0.01          # Bonferroni over five covariates: 99 % two-sided intervals
MIN_SUBGROUP_RECORDS = 30
FLOOR_REL = 0.0029              # Gate 2 MDE at 512 clusters
T_DEPART = 2.0                  # NOVA-S frozen threshold
PRIMARY_EXPECTED = 0.0011178    # N2 primary, equal-record mean (reproduction gate, tolerance 1e-6)
EXPECTED_HASHES = {
    "blockfolds": "b7934523cfc124510fdb72c0469089eecd8b8f4f3676b68552eb28618957d3bf",
    "gate3a_blocks": "6decb9fb0b9f568571a5b24f12c9dad88ef9f31be9653681e0877b2bf36e3747",
    "n0_records": "6d3b0341b5fed7f81e6b0935ced5b0a0c7c8d0398b2e450a52258c11968b75ff",
}
COVARIATES = ["log_n", "panel_class", "skew_O", "edge_O", "g1_stability"]
# prespecified signs for A5 subgroup side: +1 = high side, -1 = low side, 0 = two-sided (larger-D side)
PRESPECIFIED_SIDE = {"log_n": +1, "panel_class": 0, "skew_O": 0, "edge_O": +1, "g1_stability": -1}

# C2 lexical lipid rule (spec v0.3 §5): frozen before the run; no claim of RefMet equivalence.
LIPID_SHORTHAND = re.compile(
    r"^(PC|PE|PI|PS|PG|PA|LPC|LPE|LPI|LPS|TG|TAG|DG|DAG|MG|MAG|CE|Cer|HexCer|Hex2Cer|LacCer|GlcCer|SM|FA|CAR|CL|PIP|PIP2|PIP3)"
    r"(\s|\(|\d|-)", re.I)
LIPID_WORDS = re.compile(
    r"(phosphatidyl|ceramide|sphingomyelin|sphingosine|sphinganine|triglycer|triacylglycer|diglycer|diacylglycer|"
    r"monoacylglycer|cholesteryl|cholesterol ester|lyso|acylcarnitine|carnitine\s*\(|cardiolipin|ganglioside|"
    r"\(d1[68]:[01]/|\(t18:|/\d{1,2}:\d\))", re.I)


def sha_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def is_lipid(name: str) -> bool:
    s = name.strip()
    return bool(LIPID_SHORTHAND.match(s) or LIPID_WORDS.search(s))


def study_of(record: str) -> str:
    return record.split("_AN")[0]


# ------------------------------------------------------------------ loading
def load_blockfolds(p: Path):
    rows = []
    with p.open() as fh:
        for r in csv.DictReader(fh):
            rows.append({
                "record": r["record"], "block": int(r["block"]), "fold": int(r["fold"]),
                "n_cells": float(r["n_cells"]), "jG1": int(r["jG1"]), "jV": int(r["jV"]), "jS": int(r["jS"]),
                "t": float(r["t_inner"]), "departed": r["departed"].strip().lower() == "true",
                "sq_G1": float(r["sq_G1"]), "sq_jV": float(r["sq_jV"]),
            })
    return rows


def load_gate3a_c5(p: Path):
    """record -> (agreement count, fold count) summed over blocks."""
    out = {}
    with p.open() as fh:
        for line in fh:
            if not line.strip():
                continue
            rec = json.loads(line)
            a = sum(int(b.get("_g1_full_design_agreement", 0)) for b in rec["blocks"])
            f = sum(int(b.get("_g1_folds", 0)) for b in rec["blocks"])
            out[rec["key"]] = (a, f)
    return out


def load_n0_edge(p: Path):
    """record -> O-only relative bounded-loss improvement of kernel ridge (wide grid) over log-linear ridge (wide grid)."""
    out = {}
    with p.open() as fh:
        for r in csv.DictReader(fh):
            p1, p3 = float(r["P1_wide_bd"]), float(r["P3_wide_bd"])
            out[r["record"]] = 1.0 - p3 / p1 if p1 > 0 else float("nan")
    return out


def load_manifest_covariates(mdir: Path, records):
    """record -> (log_n, panel_class, skew_O, manifest sha256)."""
    out = {}
    for rec in records:
        mp = mdir / rec / "design_model_manifest.json"
        if not mp.exists():
            raise SystemExit(f"manifest missing for {rec}: {mp}")
        m = json.loads(mp.read_text())
        n_design = int(m["n_design"])
        lipid_blocks, skews = 0, []
        blocks = m["blocks"]
        for b in blocks:
            feats = list(b["O"]) + list(b["T"]) + list(b["P"])
            frac = sum(is_lipid(f) for f in feats) / max(1, len(feats))
            lipid_blocks += frac >= 0.5
            pre = b["preprocessing"]
            med, mean, sd = pre["median"], pre["mean_after_design_median_imputation"], pre["sd_ddof1_after_design_median_imputation"]
            vals = []
            for f in b["O"]:
                s = float(sd.get(f, float("nan")))
                if s and s > 0 and f in med and f in mean:
                    vals.append(3.0 * (float(mean[f]) - float(med[f])) / s)
            if vals:
                skews.append(float(np.mean(vals)))
        panel = 1 if blocks and lipid_blocks / len(blocks) >= 0.5 else 0
        out[rec] = (math.log(n_design), panel, float(np.mean(skews)) if skews else float("nan"), sha_file(mp), n_design)
    return out


# ------------------------------------------------------------------ outcomes
def record_table(bf, c5, edge, mani):
    """Per-record outcomes and covariates. Also reproduces the N2 primary from the block-fold file."""
    by_rec = {}
    for r in bf:
        by_rec.setdefault(r["record"], []).append(r)
    # G1-plausible set per block: picks of G1 across folds
    g1set = {}
    for r in bf:
        g1set.setdefault((r["record"], r["block"]), set()).add(r["jG1"])
    recs = []
    primary_terms = []
    for rec in sorted(by_rec):
        rows = by_rec[rec]
        # NOVA-S arm reconstruction: jS == jV when departed else jG1 -> loss sq_jV if departed else sq_G1
        sS = sum(r["sq_jV"] if r["departed"] else r["sq_G1"] for r in rows)
        sG = sum(r["sq_G1"] for r in rows)
        primary_terms.append(1.0 - sS / sG)
        diff = [r for r in rows if r["jV"] != r["jG1"]]
        if not diff:
            continue
        cells = sum(r["n_cells"] for r in diff)
        D = sum(r["sq_G1"] - r["sq_jV"] for r in diff) / cells
        ratio = 1.0 - sum(r["sq_jV"] for r in diff) / sum(r["sq_G1"] for r in diff)
        tbar = float(np.mean([r["t"] for r in diff]))
        a, f = c5.get(rec, (0, 0))
        ln, panel, skew, msha, n_design = mani[rec]
        recs.append({
            "record": rec, "study": study_of(rec), "n_design": n_design, "n_blockfolds": len(rows), "n_differing": len(diff),
            "cells_differing": cells, "D": D, "ratio_diff": ratio, "tbar": tbar,
            "log_n": ln, "panel_class": panel, "skew_O": skew, "edge_O": edge.get(rec, float("nan")),
            "g1_stability": (a / f) if f else float("nan"), "g1_agree": a, "g1_folds": f, "manifest_sha256": msha,
        })
    primary = float(np.mean(primary_terms))
    return recs, primary, g1set


def wls(X, y, w):
    """Weighted least squares by normal equations (solve), falling back to lstsq. Rows with any non-finite value are dropped."""
    ok = np.isfinite(X).all(1) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    X, y, w = X[ok], y[ok], w[ok]
    XtW = X.T * w
    try:
        return np.linalg.solve(XtW @ X, XtW @ y)
    except np.linalg.LinAlgError:
        W = np.sqrt(w)[:, None]
        beta, *_ = np.linalg.lstsq(X * W, y * W[:, 0], rcond=None)
        return beta


def cluster_boot(recs, func, rng, B=B_BOOT):
    """func(list_of_recs) -> np.array of statistics. Resample study clusters with replacement."""
    clusters = {}
    for r in recs:
        clusters.setdefault(r["study"], []).append(r)
    keys = sorted(clusters)
    base = func(recs)
    draws = np.empty((B, base.size))
    for b in range(B):
        idx = rng.integers(0, len(keys), len(keys))
        sample = [r for i in idx for r in clusters[keys[i]]]
        draws[b] = func(sample)
    lo = np.nanpercentile(draws, 100 * ALPHA_TWO_SIDED / 2, axis=0)
    hi = np.nanpercentile(draws, 100 * (1 - ALPHA_TWO_SIDED / 2), axis=0)
    return base, lo, hi, len(keys)


def design_matrix(recs, outcome="D", with_t=True):
    cov = np.array([[r[c] for c in COVARIATES] for r in recs], float)
    mu, sd = np.nanmean(cov, 0), np.nanstd(cov, 0, ddof=1)
    sd[sd == 0] = 1.0
    Z = (cov - mu) / sd
    cols = [np.ones(len(recs))]
    if with_t:
        cols.append(np.array([r["tbar"] for r in recs]))
    X = np.column_stack(cols + [Z])
    y = np.array([r[outcome] for r in recs])
    w = np.array([r["cells_differing"] for r in recs], float)
    return X, y, w, mu, sd


def regression_stats(recs, outcome, with_t):
    def f(sample):
        X, y, w, _, _ = design_matrix(sample, outcome, with_t)
        return wls(X, y, w)
    return f


def primary_from_rows(rows, depart_mask):
    """Equal-record mean of per-record ratio of sums with the NOVA-S arm departing only where depart_mask is True."""
    by = {}
    for r, d in zip(rows, depart_mask):
        s = by.setdefault(r["record"], [0.0, 0.0])
        s[0] += r["sq_jV"] if d else r["sq_G1"]
        s[1] += r["sq_G1"]
    return float(np.mean([1 - a / b for a, b in by.values()]))


def analyse(bf, recs, primary, g1set, rng):
    out = {"n_records": len(recs), "n_clusters": len({r["study"] for r in recs}),
           "n2_primary_reproduced": primary}
    # ---- A1 descriptive by tertile
    a1 = {}
    for c in COVARIATES:
        vals = np.array([r[c] for r in recs], float)
        if c == "panel_class":
            groups = {"broad": [r for r in recs if r[c] == 0], "lipid": [r for r in recs if r[c] == 1]}
        else:
            q1, q2 = np.nanpercentile(vals, [100 / 3, 200 / 3])
            groups = {"T1_low": [r for r in recs if r[c] <= q1], "T2": [r for r in recs if q1 < r[c] <= q2],
                      "T3_high": [r for r in recs if r[c] > q2]}
        a1[c] = {}
        for g, sub in groups.items():
            if not sub:
                continue
            base, lo, hi, k = cluster_boot(sub, lambda s: np.array([np.average([r["D"] for r in s], weights=[r["cells_differing"] for r in s]),
                                                                    np.mean([r["ratio_diff"] for r in s])]), rng, B=2000)
            a1[c][g] = {"records": len(sub), "clusters": k, "D_cellweighted": base[0], "D_99": [lo[0], hi[0]],
                        "ratio_equal_record": base[1], "ratio_99": [lo[1], hi[1]]}
    out["A1"] = a1
    # ---- A2 partial association
    X, y, w, mu, sd = design_matrix(recs, "D", True)
    ok = np.isfinite(X).all(1) & np.isfinite(y) & np.isfinite(w)
    nonfinite = {c: int((~np.isfinite(X[:, 2 + i])).sum()) for i, c in enumerate(COVARIATES)}
    nonfinite["D"] = int((~np.isfinite(y)).sum()); nonfinite["tbar"] = int((~np.isfinite(X[:, 1])).sum())
    print("non-finite records per column (excluded from A2/A3):", nonfinite)
    Zc = X[ok][:, 2:]
    corr = np.corrcoef(Zc.T)
    cond = float(np.linalg.cond(X[ok]))
    base, lo, hi, k = cluster_boot([r for r, o in zip(recs, ok) if o], regression_stats(recs, "D", True), rng)
    names = ["intercept", "tbar"] + COVARIATES
    a2 = {n: {"coef": float(base[i]), "ci99": [float(lo[i]), float(hi[i])],
              "excludes_zero": bool(lo[i] > 0 or hi[i] < 0)} for i, n in enumerate(names)}
    out["A2"] = {"coefficients": a2, "condition_number": cond, "clusters": k, "records_used": int(ok.sum()),
                 "nonfinite_excluded": nonfinite,
                 "covariate_correlations": {COVARIATES[i]: {COVARIATES[j]: float(corr[i, j]) for j in range(5)} for i in range(5)},
                 "standardisation": {c: {"mean": float(mu[i]), "sd": float(sd[i])} for i, c in enumerate(COVARIATES)}}
    # ---- A3 covariate -> t
    base3, lo3, hi3, _ = cluster_boot([r for r, o in zip(recs, ok) if o], regression_stats(recs, "tbar", False), rng)
    names3 = ["intercept"] + COVARIATES
    out["A3"] = {n: {"coef": float(base3[i]), "ci99": [float(lo3[i]), float(hi3[i])],
                     "excludes_zero": bool(lo3[i] > 0 or hi3[i] < 0)} for i, n in enumerate(names3)}
    # ---- A4' inside / outside the G1-plausible set
    diff_rows = [r for r in bf if r["jV"] != r["jG1"]]
    inside = [r["jV"] in g1set[(r["record"], r["block"])] for r in diff_rows]
    rec_of = {r["record"]: study_of(r["record"]) for r in diff_rows}

    def side_stats(rows_sub):
        rec_rows = {}
        for r in rows_sub:
            rec_rows.setdefault(r["record"], []).append(r)
        pseudo = [{"study": rec_of[k], "rows": v} for k, v in rec_rows.items()]

        def f(sample):
            rr = [r for p in sample for r in p["rows"]]
            cells = sum(r["n_cells"] for r in rr)
            D = sum(r["sq_G1"] - r["sq_jV"] for r in rr) / cells
            dep = [r for r in rr if r["departed"]]
            Dd = (sum(r["sq_G1"] - r["sq_jV"] for r in dep) / sum(r["n_cells"] for r in dep)) if dep else float("nan")
            return np.array([D, Dd])
        base, lo, hi, k = cluster_boot(pseudo, f, rng)
        return {"block_folds": len(rows_sub), "departures": sum(r["departed"] for r in rows_sub),
                "records": len(rec_rows), "clusters": k,
                "D_all_differing": float(base[0]), "D_all_99": [float(lo[0]), float(hi[0])],
                "D_departures": float(base[1]), "D_departures_99": [float(lo[1]), float(hi[1])]}
    ins = [r for r, i in zip(diff_rows, inside) if i]
    outs = [r for r, i in zip(diff_rows, inside) if not i]
    dep_all = [r["departed"] for r in bf]
    key_inside = {(r["record"], r["block"], r["fold"]) for r in ins}
    contrib_in = primary_from_rows(bf, [d and (r["record"], r["block"], r["fold"]) in key_inside for r, d in zip(bf, dep_all)])
    contrib_out = primary_from_rows(bf, [d and (r["record"], r["block"], r["fold"]) not in key_inside for r, d in zip(bf, dep_all)])
    out["A4"] = {"inside_G1set": side_stats(ins), "outside_G1set": side_stats(outs),
                 "primary_contribution_inside_only": contrib_in, "primary_contribution_outside_only": contrib_out,
                 "primary_full": primary_from_rows(bf, dep_all)}
    # ---- A5 prespecified median split for covariates passing A2
    a5 = {}
    for c in COVARIATES:
        if not out["A2"]["coefficients"][c]["excludes_zero"]:
            continue
        vals = np.array([r[c] for r in recs], float)
        med = float(np.nanmedian(vals))
        side = PRESPECIFIED_SIDE[c]
        if c == "panel_class":
            hi_g = [r for r in recs if r[c] == 1]; lo_g = [r for r in recs if r[c] == 0]
        else:
            hi_g = [r for r in recs if r[c] > med]; lo_g = [r for r in recs if r[c] <= med]
        if side == 0:
            side = +1 if np.mean([r["D"] for r in hi_g]) >= np.mean([r["D"] for r in lo_g]) else -1
        sub = hi_g if side > 0 else lo_g
        sub_recs = {r["record"] for r in sub}
        rows_sub = [r for r in bf if r["record"] in sub_recs]
        studies = {r["record"]: study_of(r["record"]) for r in rows_sub}
        by = {}
        for r in rows_sub:
            by.setdefault(r["record"], []).append(r)
        pseudo = [{"study": studies[k], "rows": v} for k, v in by.items()]

        def f(sample):
            terms = []
            for p in sample:
                sS = sum(r["sq_jV"] if r["departed"] else r["sq_G1"] for r in p["rows"])
                sG = sum(r["sq_G1"] for r in p["rows"])
                terms.append(1 - sS / sG)
            return np.array([np.mean(terms)])
        base, lo, hi, k = cluster_boot(pseudo, f, rng)
        a5[c] = {"side": "high" if side > 0 else "low", "cut": med if c != "panel_class" else "class",
                 "records": len(sub), "clusters": k, "delta_acq_departures": float(base[0]),
                 "lower99": float(lo[0]), "upper99": float(hi[0]),
                 "meets_floor": bool(len(sub) >= MIN_SUBGROUP_RECORDS and base[0] >= FLOOR_REL and lo[0] > 0)}
    out["A5"] = a5
    # ---- A6 t >= 4 listing
    hi_t = [r for r in diff_rows if r["t"] >= 4]
    rec_cov = {r["record"]: r for r in recs}
    out["A6"] = [{"record": r["record"], "study": study_of(r["record"]), "block": r["block"], "fold": r["fold"], "t": r["t"],
                  "gain_rel": 1 - r["sq_jV"] / r["sq_G1"], "inside_G1set": r["jV"] in g1set[(r["record"], r["block"])],
                  **{c: rec_cov[r["record"]][c] for c in COVARIATES if r["record"] in rec_cov}} for r in hi_t]
    # ---- reading (mechanical, per spec §7)
    reading = []
    for c in COVARIATES:
        a2c = out["A2"]["coefficients"][c]
        a3c = out["A3"][c]
        if a2c["excludes_zero"]:
            sub = a5.get(c, {})
            sign_ok = {"log_n": a2c["coef"] > 0, "edge_O": a2c["coef"] > 0, "g1_stability": a2c["coef"] < 0}.get(c, True)
            if sub.get("meets_floor") and sign_ok:
                reading.append(f"{c}: REGIME IDENTIFIED (exploratory), subgroup {sub['side']} side, {sub['records']} records")
            elif c == "log_n" and a2c["coef"] < 0:
                reading.append("log_n: negative association; read as small-n estimation noise (GDOA-0), not a regime")
            else:
                reading.append(f"{c}: passes A2 but subgroup criterion not met (records {sub.get('records')}, "
                               f"delta {sub.get('delta_acq_departures')}, lower99 {sub.get('lower99')})")
        elif a3c["excludes_zero"]:
            reading.append(f"{c}: predicts departure only (A3), not advantage given t")
        else:
            reading.append(f"{c}: no association")
    if "g1_stability" in a5 and a5["g1_stability"].get("meets_floor"):
        a4i = out["A4"]["inside_G1set"]
        reading.append("stability reading qualifier: inside-G1set departures D_99 = %s -> %s" % (
            a4i["D_departures_99"], "stability (resolution helps)" if a4i["D_departures_99"][0] > 0 else "stability without cost (tie-breaking)"))
    out["reading"] = reading
    return out


# ------------------------------------------------------------------ run
def cmd_run(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    paths = {"blockfolds": Path(a.blockfolds), "gate3a_blocks": Path(a.gate3a_blocks), "n0_records": Path(a.n0_records)}
    hashes = {k: sha_file(p) for k, p in paths.items()}
    for k, h in hashes.items():
        if h != EXPECTED_HASHES[k]:
            raise SystemExit(f"REFUSED: {k} sha256 {h} != expected {EXPECTED_HASHES[k]}")
    spec = Path(a.spec)
    me = Path(__file__).resolve()
    bf = load_blockfolds(paths["blockfolds"])
    if len(bf) != 4760:
        raise SystemExit(f"REFUSED: expected 4,760 block-folds, found {len(bf)}")
    c5 = load_gate3a_c5(paths["gate3a_blocks"])
    edge = load_n0_edge(paths["n0_records"])
    records = sorted({r["record"] for r in bf})
    mani = load_manifest_covariates(Path(a.manifests), records)
    recs, primary, g1set = record_table(bf, c5, edge, mani)
    if abs(primary - PRIMARY_EXPECTED) > 1e-6:
        raise SystemExit(f"REFUSED: N2 primary not reproduced from block-fold file: {primary:.7f} vs {PRIMARY_EXPECTED}")
    fields = ["record", "study", "n_design", "n_blockfolds", "n_differing", "cells_differing", "D", "ratio_diff", "tbar",
              "log_n", "panel_class", "skew_O", "edge_O", "g1_stability", "g1_agree", "g1_folds", "manifest_sha256"]
    with (out / "n2_regime_records.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader()
        for r in recs:
            w.writerow({k: r[k] for k in fields})
    print(f"records table written ({len(recs)} records); N2 primary reproduced: {primary:.7f}")
    rng = np.random.default_rng(SEED)
    res = analyse(bf, recs, primary, g1set, rng)
    res.update({"script": me.name, "version": VERSION, "script_sha256": sha_file(me), "spec": spec.name,
                "spec_sha256": sha_file(spec) if spec.exists() else None, "inputs": hashes,
                "manifests_sha256_count": len({r['manifest_sha256'] for r in recs}),
                "seed": SEED, "bootstrap": B_BOOT, "alpha_two_sided": ALPHA_TWO_SIDED, "archive_use": 11,
                "classification": "exploratory, descriptive; no selector trained; no evaluation row read",
                "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())})
    (out / "n2_regime_summary.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    print(json.dumps({k: res[k] for k in ("n_records", "n_clusters", "n2_primary_reproduced", "reading")}, indent=2))
    print("A2:", json.dumps(res["A2"]["coefficients"], indent=1, default=float))
    print("A4:", json.dumps(res["A4"], indent=1, default=float))
    print("summary sha256", sha_file(out / "n2_regime_summary.json"))


# ------------------------------------------------------------------ selftest (synthetic; no real file read)
def _synth(rng, n_rec=60, planted=None, n_blocks=4, nfold=5):
    """Synthetic block-folds and covariates. planted: dict covariate -> effect on D (per unit SD)."""
    bf, c5, edge, mani = [], {}, {}, {}
    for i in range(n_rec):
        rec = f"ST{100000 + i // 2:06d}_AN{i:06d}"     # two records per study -> clusters
        cov = {"log_n": rng.normal(4, 0.5), "panel_class": int(rng.random() < 0.3), "skew_O": rng.normal(1, 0.5),
               "edge_O": rng.normal(0.05, 0.03), "g1_stability": float(np.clip(rng.normal(0.7, 0.15), 0, 1))}
        eff = sum((planted or {}).get(c, 0.0) * (cov[c] - {"log_n": 4, "panel_class": 0.3, "skew_O": 1, "edge_O": 0.05, "g1_stability": 0.7}[c])
                  / {"log_n": 0.5, "panel_class": 0.46, "skew_O": 0.5, "edge_O": 0.03, "g1_stability": 0.15}[c] for c in cov)
        for b in range(n_blocks):
            g1picks = rng.integers(0, 20, 2)
            for f in range(nfold):
                jG1 = int(rng.choice(g1picks)); jV = int(rng.integers(0, 20)); t = float(abs(rng.normal(1.2, 1)))
                n_cells = 200.0
                sqG = n_cells * (0.7 + rng.normal(0, 0.05))
                sqV = sqG - n_cells * (eff + rng.normal(0, 0.03)) if jV != jG1 else sqG
                dep = jV != jG1 and t > T_DEPART
                bf.append({"record": rec, "block": b, "fold": f, "n_cells": n_cells, "jG1": jG1, "jV": jV,
                           "jS": jV if dep else jG1, "t": t, "departed": dep, "sq_G1": sqG, "sq_jV": sqV})
        folds = n_blocks * nfold
        c5[rec] = (int(round(cov["g1_stability"] * folds)), folds)
        edge[rec] = cov["edge_O"]
        mani[rec] = (cov["log_n"], cov["panel_class"], cov["skew_O"], "synthetic", int(round(math.exp(cov["log_n"]))))
    return bf, c5, edge, mani


def cmd_selftest(a):
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("PASS " if cond else "FAIL ") + msg); ok &= bool(cond)
    rng = np.random.default_rng(1)
    # (i) planted association recovered
    bf, c5, edge, mani = _synth(rng, planted={"g1_stability": -0.02, "edge_O": 0.02})
    recs, primary, g1set = record_table(bf, c5, edge, mani)
    X, y, w, _, _ = design_matrix(recs, "D", True)
    base, lo, hi, k = cluster_boot(recs, regression_stats(recs, "D", True), np.random.default_rng(2), B=1000)
    i5, i4 = 2 + COVARIATES.index("g1_stability"), 2 + COVARIATES.index("edge_O")
    check(hi[i5] < 0 and lo[i4] > 0, f"planted C5 (-) and C4 (+) recovered: C5 {base[i5]:.4f} [{lo[i5]:.4f},{hi[i5]:.4f}], C4 {base[i4]:.4f} [{lo[i4]:.4f},{hi[i4]:.4f}]")
    check(k == 30, f"two records per study -> 30 clusters ({k})")
    # (ii) null false-positive rate across 20 replicates
    fp = 0
    for rep in range(20):
        r2 = np.random.default_rng(100 + rep)
        bf0, c50, e0, m0 = _synth(r2)
        recs0, _, _ = record_table(bf0, c50, e0, m0)
        b0, l0, h0, _ = cluster_boot(recs0, regression_stats(recs0, "D", True), np.random.default_rng(200 + rep), B=500)
        fp += int(sum((l0[2:] > 0) | (h0[2:] < 0)))
    check(fp <= 3, f"null: {fp} of 100 covariate tests exclude 0 at 99 % (allowed <= 3)")
    # (iii) cluster resampling keeps a study's records together
    seen_split = False
    clusters = {}
    for r in recs:
        clusters.setdefault(r["study"], []).append(r["record"])
    keys = sorted(clusters)
    rr = np.random.default_rng(3)
    for _ in range(200):
        idx = rr.integers(0, len(keys), len(keys))
        sample = [rec for i in idx for rec in clusters[keys[i]]]
        for s, members in clusters.items():
            cnt = [sample.count(m) for m in members]
            if len(set(cnt)) > 1:
                seen_split = True
    check(not seen_split, "cluster bootstrap never splits a study's records")
    # (iv) primary reconstruction is self-consistent on synthetic data
    full = primary_from_rows(bf, [r["departed"] for r in bf])
    check(abs(full - primary) < 1e-12, "NOVA-S primary reconstruction consistent between the two code paths")
    # (v) lexical lipid rule on known names
    pos = ["PC 16:0/18:1", "Cer(d18:0/16:0)", "Ganglioside GM3 (d18:1/23:0)", "LPC 16:0", "TG 52:2", "Sphingomyelin (d18:1/16:0)",
           "CL(20:1(11Z)/18:2(9Z,12Z)/18:1(11Z)/18:1(11Z))", "Lysophosphatidylcholine"]
    neg = ["L-Lysine", "Cinnamaldehyde", "Adenosine monophosphate", "Pivampicillin", "Cefalexin", "Fatty acid binding protein",
           "Carnitine", "Pepstatin", "Salicin", "Caffeine"]
    check(all(is_lipid(x) for x in pos) and not any(is_lipid(x) for x in neg), "lexical lipid rule on 8 positive and 10 negative names")
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    r = sub.add_parser("run")
    r.add_argument("--blockfolds", required=True); r.add_argument("--gate3a-blocks", required=True)
    r.add_argument("--n0-records", required=True); r.add_argument("--manifests", required=True)
    r.add_argument("--spec", default="N2_REGIME_DIAGNOSTIC_SPEC_v0.3.md"); r.add_argument("--out", default="n2_regime_out")
    a = ap.parse_args()
    {"selftest": cmd_selftest, "run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
