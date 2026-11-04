#!/usr/bin/env python3
"""NOVA stage N2: development evaluation of NOVA-V/U/S against G1 on the 186-record archive (design rows only).

Implements N2_PROTOCOL_v1.0.md. Per block and outer fold k (Gate 3A folds, seed 20261019):
  on D_-k only: inner 5-fold cross-fit (seed 20261019 + k) gives, for the nuisance family, per-row scores
    psi_ij = l_i(O) - l_i(O u {j})  (squared loss, floored standardization; charter Lemma 1 = the orthogonal score)
    at hyperparameters chosen by the inner-CV total of each model; V^_j = mean_i psi_ij, se_j its standard error.
    Picks: NOVA-V = argmax V^_j; NOVA-U = argmax (V^_j - lambda_u se_j); NOVA-S = NOVA-V's pick if the paired mean
    of psi_i,jV - psi_i,jG1 exceeds 2 of its standard errors, else the G1 pick; G1 re-estimated on D_-k (frozen
    Paper IV rule). Ties: lowest index.
  on D_k: the family refitted on D_-k at the chosen hyperparameters for each picked candidate; losses pooled.
Endpoint per family f and variant: Delta_acq = 1 - R_f(O u {j_variant}) / R_f(O u {j_G1}), predictor fixed in both
arms; equal-record mean; record bootstrap 10,000. Primary: NOVA-S with the RBF nuisance under the floored squared
loss. Sensitivity: the Gate 3A cap min(e^2, 25) on the same predictions. Ablations: NOVA-V, NOVA-U (rbf); the
additive nuisance for all three variants.
Loss standardization: each target divided by max(sd_t, FLOOR * median_t sd_t), sd from the archive's design-row
constants (identical in both arms; not selection dependent).

    python3 -I p6_n2.py selftest
    python3 -I p6_n2.py run --bench <...> --protocol N2_PROTOCOL_v1.0.md --protocol-sha256 <hash> [--out n2_out] [--workers 4]
Requires p5_gate3a.py v0.1.0 (SHA-256 checked) in the same directory for the archive binding, log rule, features,
folds and G1 rule.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import sys
import time

import numpy as np

VERSION = "0.1.0"
GATE3A_SHA = "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"
SEED = 20261019
K_OUT, K_IN = 5, 5
LAMBDAS = (1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4)
BWS = (0.5, 1.0, 2.0)
LAMBDA_U = 1.0
CLAIM_SE = 2.0
FLOOR = 0.01
CAP = 25.0
MDE = 0.0029
ALPHA = 0.025
B_REPS = 10_000
FAMILIES = ("rbf", "add")
VARIANTS = ("NOVA-V", "NOVA-U", "NOVA-S")
PRIMARY = ("rbf", "NOVA-S")
_G = None


def sha_file(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def gate3a():
    global _G
    if _G is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "p5_gate3a.py")
        if sha_file(path) != GATE3A_SHA:
            raise SystemExit(f"p5_gate3a.py has SHA-256 {sha_file(path)}, expected {GATE3A_SHA}")
        spec = importlib.util.spec_from_file_location("p5_gate3a_frozen", path)
        g = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(g)
        _G = g
    return _G


# ---------------------------------------------------------------- kernels
def rbf_k(A, B, bw):
    d2 = ((A[:, None, :] - B[None, :, :]) ** 2).sum(-1)
    return np.exp(-d2 / (2.0 * bw ** 2))


def med_dist(A):
    if len(A) < 2:
        return 1.0
    d2 = ((A[:, None, :] - A[None, :, :]) ** 2).sum(-1)
    m = np.sqrt(np.median(d2[np.triu_indices(len(A), 1)]))
    return m if m > 1e-12 else 1.0


def krr_all(Ktr, Ytr, Kte):
    my = Ytr.mean(axis=0)
    w, Q = np.linalg.eigh(Ktr)
    QtY = Q.T @ (Ytr - my)
    return [Kte @ (Q @ (QtY / (np.maximum(w, 0.0) + lam)[:, None])) + my for lam in LAMBDAS]


def n_hyper():
    return len(BWS) * len(LAMBDAS)


def predict_grid(family, Otr, Ote, Ytr, xtr=None, xte=None, base_tr=None, base_te=None):
    """(n_hyper, n_te, q) predictions. 'rbf': KRR on log-rule standardized (O[, x]). 'add' with a candidate:
    base + KRR on residuals with kernel K_x + K_O o K_x; without a candidate: same as 'rbf' on O."""
    g = gate3a()
    out = []
    if family == "rbf" or xtr is None:
        Ftr, Fte = g.features("P3", Otr, Ote, xtr, xte)
        md = med_dist(Ftr)
        for bw in BWS:
            out += krr_all(rbf_k(Ftr, Ftr, bw * md), Ytr, rbf_k(Fte, Ftr, bw * md))
        return np.stack(out)
    Otr_s, Ote_s = g.features("P1", Otr, Ote)
    xtr_s, xte_s = g.features("P1", xtr[:, None], xte[:, None])
    if xtr_s.shape[1] == 0:                                   # candidate dropped by the log rule: no correction
        return np.repeat(base_te[None], n_hyper(), axis=0)
    mdO, mdx = med_dist(Otr_s), med_dist(xtr_s)
    R = Ytr - base_tr
    for bw in BWS:
        KO_tr, KO_te = rbf_k(Otr_s, Otr_s, bw * mdO), rbf_k(Ote_s, Otr_s, bw * mdO)
        Kx_tr, Kx_te = rbf_k(xtr_s, xtr_s, bw * mdx), rbf_k(xte_s, xtr_s, bw * mdx)
        out += [base_te + p for p in krr_all(Kx_tr + KO_tr * Kx_tr, R, Kx_te + KO_te * Kx_te)]
    return np.stack(out)


def row_loss(P, Y):
    """(n_hyper, n) per-row squared loss and capped loss."""
    E2 = (P - Y[None]) ** 2
    return E2.sum(axis=2), np.minimum(E2, CAP).sum(axis=2)


# ---------------------------------------------------------------- one block
def inner_scores(family, O, X, Y, seed):
    """On the given rows: inner cross-fitted per-row scores psi (n, m) at the inner-CV-chosen hyper of each model,
    plus the chosen hyper indices (m + 1)."""
    n, m = len(Y), X.shape[1]
    fold = np.random.default_rng(seed).permutation(n) % K_IN
    L = np.full((m + 1, n_hyper(), n), np.nan)
    for f in range(K_IN):
        tr, te = fold != f, fold == f
        if te.sum() == 0 or tr.sum() < 3:
            continue
        P0 = predict_grid(family, O[tr], O[te], Y[tr])
        L[0][:, te] = row_loss(P0, Y[te])[0]
        base_tr = predict_grid(family, O[tr], O[tr], Y[tr]) if family == "add" else None
        for j in range(m):
            if family == "add":
                b = int(np.nansum(L[0], axis=1).argmin()) if f == K_IN - 1 else 0
                # base at a provisional hyper during inner folds: use the mid-grid hyper (bw 1, lambda 1) for the
                # base inside inner folds; the final base hyper is chosen afterwards (protocol section 3)
                bi = len(LAMBDAS) + 2
                Pj = predict_grid(family, O[tr], O[te], Y[tr], X[tr, j], X[te, j], base_tr[bi], P0[bi])
            else:
                Pj = predict_grid(family, O[tr], O[te], Y[tr], X[tr, j], X[te, j])
            L[j + 1][:, te] = row_loss(Pj, Y[te])[0]
    tot = np.nansum(L, axis=2)
    best = tot.argmin(axis=1)
    l0 = L[0][best[0]]
    psi = np.stack([l0 - L[j + 1][best[j + 1]] for j in range(m)], axis=1)
    return psi, best


def run_block(O, X, Y, subjects, seed, family):
    g = gate3a()
    fold = g.folds_for(subjects, seed)
    m = X.shape[1]
    res = {"sq": {}, "bd": {}, "n": 0.0, "picks": [], "claims": 0, "folds": 0, "agree": {v: 0 for v in VARIANTS}}
    arms = ("G1",) + VARIANTS
    for a in arms:
        res["sq"][a] = 0.0
        res["bd"][a] = 0.0
    for k in range(K_OUT):
        tr, te = fold != k, fold == k
        if te.sum() == 0:
            continue
        jg = g.g1_pick(O[tr], Y[tr], X[tr])
        psi, best = inner_scores(family, O[tr], X[tr], Y[tr], SEED + k)
        ok = ~np.isnan(psi).any(axis=1)
        psi = psi[ok]
        Vh = psi.mean(axis=0)
        se = psi.std(axis=0, ddof=1) / np.sqrt(len(psi))
        jv = int(np.flatnonzero(Vh == Vh.max())[0])
        u = Vh - LAMBDA_U * se
        ju = int(np.flatnonzero(u == u.max())[0])
        claim = False
        if jv != jg:
            d = psi[:, jv] - psi[:, jg]
            sd = d.std(ddof=1) / np.sqrt(len(d))
            claim = bool(d.mean() > CLAIM_SE * sd)
        js = jv if claim else jg
        picks = {"G1": jg, "NOVA-V": jv, "NOVA-U": ju, "NOVA-S": js}
        res["picks"].append(picks)
        res["claims"] += int(claim)
        res["folds"] += 1
        for v in VARIANTS:
            res["agree"][v] += int(picks[v] == jg)
        base_tr = predict_grid(family, O[tr], O[tr], Y[tr])[best[0]] if family == "add" else None
        P0te = predict_grid(family, O[tr], O[te], Y[tr])[best[0]] if family == "add" else None
        cache = {}
        for a, j in picks.items():
            if j not in cache:
                if family == "add":
                    P = predict_grid(family, O[tr], O[te], Y[tr], X[tr, j], X[te, j], base_tr, P0te)[best[j + 1]]
                else:
                    P = predict_grid(family, O[tr], O[te], Y[tr], X[tr, j], X[te, j])[best[j + 1]]
                E2 = (P - Y[te]) ** 2
                cache[j] = (float(E2.sum()), float(np.minimum(E2, CAP).sum()))
            res["sq"][a] += cache[j][0]
            res["bd"][a] += cache[j][1]
        res["n"] += float(te.sum() * Y.shape[1])
    return res


def block_arrays_floored(b, frame):
    O = frame[b["O"]].to_numpy(float)
    X = frame[b["P"]].to_numpy(float)
    T = frame[b["T"]].to_numpy(float)
    sd = np.array([b["target_reference"]["sd_ddof1_observed_design"][t] for t in b["T"]], float)
    sd = np.where(np.isfinite(sd) & (sd > 0), sd, 0.0)
    floor = FLOOR * np.median(sd[sd > 0]) if (sd > 0).any() else 1.0
    return O, X, T / np.maximum(sd, floor)


def run_record(task):
    bench, study = task
    g = gate3a()
    audit = g._bind(bench)
    key, m, _df, dev, _test = audit.load_record(study)
    del _df, _test
    ids = list(dev.index.to_numpy().astype(str))
    out = {"key": key, "blocks": 0, "unusable": [], "fam": {}}
    for fam in FAMILIES:
        out["fam"][fam] = {"sq": {a: 0.0 for a in ("G1",) + VARIANTS}, "bd": {a: 0.0 for a in ("G1",) + VARIANTS},
                           "n": 0.0, "claims": 0, "folds": 0, "agree": {v: 0 for v in VARIANTS}}
    for b in m["blocks"]:
        O, X, Y = block_arrays_floored(b, dev)
        if not (np.isfinite(O).all() and np.isfinite(X).all() and np.isfinite(Y).all()):
            out["unusable"].append({"block": b["block_id"], "reason": "non-finite design value"})
            continue
        try:
            rs = {fam: run_block(O, X, Y, ids, SEED, fam) for fam in FAMILIES}
        except (ValueError, np.linalg.LinAlgError) as e:
            out["unusable"].append({"block": b["block_id"], "reason": str(e)[:120]})
            continue
        out["blocks"] += 1
        for fam, r in rs.items():
            acc = out["fam"][fam]
            for loss in ("sq", "bd"):
                for a in acc[loss]:
                    acc[loss][a] += r[loss][a]
            acc["n"] += r["n"]
            acc["claims"] += r["claims"]
            acc["folds"] += r["folds"]
            for v in VARIANTS:
                acc["agree"][v] += r["agree"][v]
    return out


# ---------------------------------------------------------------- aggregation
def boot(vals, q, seed=SEED, reps=B_REPS):
    v = np.asarray(vals, float)
    rng = np.random.default_rng(seed)
    means = v[rng.integers(0, len(v), size=(reps, len(v)))].mean(axis=1)
    return float(np.quantile(means, q))


def aggregate(records):
    out = {}
    for fam in FAMILIES:
        for loss in ("sq", "bd"):
            for v in VARIANTS:
                d = [1 - r["fam"][fam][loss][v] / r["fam"][fam][loss]["G1"] for r in records
                     if r["blocks"] and r["fam"][fam][loss]["G1"] > 0]
                row = {"delta_acq": float(np.mean(d)), "lower_97.5": boot(d, ALPHA), "upper_97.5": boot(d, 1 - ALPHA),
                       "records": len(d), "records_positive": int(np.sum(np.array(d) > 0))}
                row["qualifies"] = bool((fam, v) == PRIMARY and loss == "sq" and row["delta_acq"] >= MDE
                                        and row["lower_97.5"] > 0)
                out[f"{fam}|{v}|{loss}"] = row
        folds = sum(r["fam"][fam]["folds"] for r in records)
        out[f"{fam}|claim_rate"] = sum(r["fam"][fam]["claims"] for r in records) / max(folds, 1)
        for v in VARIANTS:
            out[f"{fam}|{v}|agreement_with_G1"] = sum(r["fam"][fam]["agree"][v] for r in records) / max(folds, 1)
    return out


# ---------------------------------------------------------------- self-tests
def synth(rng, n, kind):
    g = gate3a()
    O, X, Y = g.synth_block(rng, n, "null")
    if kind == "quadratic":
        z = rng.normal(size=n)
        X[:, 5] = z                                           # negative values: identity branch; Cov(T, z) = 0
        Y += 1.5 * (z ** 2 - 1.0)[:, None]                   # on every target: a dominant nonlinear candidate
        Y = (Y - Y.mean(0)) / Y.std(0, ddof=1)
    return O, X, Y


def selftest(a):
    rng = np.random.default_rng(SEED)
    rep = {}
    # planted: NOVA-S (rbf) should pick candidate 5 in >= 80% of block-folds; G1 should not
    picks = {"NOVA-S": 0, "NOVA-V": 0, "G1": 0, "folds": 0, "claims": 0, "claims_correct": 0}
    for b in range(a.planted_blocks):
        O, X, Y = synth(rng, 150, "quadratic")
        r = run_block(O, X, Y, [f"s{i}" for i in range(len(Y))], SEED, "rbf")
        for p in r["picks"]:
            picks["folds"] += 1
            picks["NOVA-S"] += int(p["NOVA-S"] == 5)
            picks["NOVA-V"] += int(p["NOVA-V"] == 5)
            if p["NOVA-S"] != p["G1"]:
                picks["claims"] += 1
                picks["claims_correct"] += int(p["NOVA-S"] == 5)
            picks["G1"] += int(p["G1"] == 5)
    rep["planted"] = {k: v / picks["folds"] for k, v in picks.items() if k not in ("folds", "claims", "claims_correct")}
    rep["planted"]["claim_rate"] = picks["claims"] / picks["folds"]
    rep["planted"]["claims_correct"] = picks["claims_correct"] / max(picks["claims"], 1)
    # null: claim rate <= 5%, Delta_acq(NOVA-S) within 2 SE of 0
    records = []
    for i in range(a.null_records):
        acc = {"key": f"N{i}", "blocks": 2, "unusable": [], "fam": {}}
        for fam in ("rbf",):
            acc["fam"][fam] = {"sq": {x: 0.0 for x in ("G1",) + VARIANTS}, "bd": {x: 0.0 for x in ("G1",) + VARIANTS},
                               "n": 0.0, "claims": 0, "folds": 0, "agree": {v: 0 for v in VARIANTS}}
            for b in range(2):
                O, X, Y = synth(rng, 80, "null")
                r = run_block(O, X, Y, [f"s{j}" for j in range(len(Y))], SEED, fam)
                for loss in ("sq", "bd"):
                    for x in acc["fam"][fam][loss]:
                        acc["fam"][fam][loss][x] += r[loss][x]
                acc["fam"][fam]["n"] += r["n"]
                acc["fam"][fam]["claims"] += r["claims"]
                acc["fam"][fam]["folds"] += r["folds"]
                for v in VARIANTS:
                    acc["fam"][fam]["agree"][v] += r["agree"][v]
        records.append(acc)
    global FAMILIES
    FAMILIES = ("rbf",)
    agg = aggregate(records)
    FAMILIES = ("rbf", "add")
    d = agg["rbf|NOVA-S|sq"]
    se = (d["upper_97.5"] - d["lower_97.5"]) / (2 * 1.96)
    rep["null"] = {"claim_rate": agg["rbf|claim_rate"], "delta_acq_NOVA-S": d["delta_acq"], "se": se,
                   "delta_acq_NOVA-V": agg["rbf|NOVA-V|sq"]["delta_acq"]}
    checks = {"planted: NOVA-V picks the square candidate in >= 80% of block-folds": rep["planted"]["NOVA-V"] >= 0.8,
              "planted: when NOVA-S claims, the claim is the square candidate in >= 90% of claims (power is reported, not tested)":
              rep["planted"]["claims_correct"] >= 0.9 and picks["claims"] > 0,
              "planted: G1 picks it in <= 30%": rep["planted"]["G1"] <= 0.3,
              "null: claim rate <= 5%": rep["null"]["claim_rate"] <= 0.05,
              "null: Delta_acq(NOVA-S) within 2 SE of 0": abs(d["delta_acq"]) <= 2 * se + 1e-12,
              "floor: a zero-SD target is floored, not divided by zero": _floor_check()}
    for k, v in checks.items():
        print("OK  " if v else "FAIL", k)
    print(json.dumps(rep, indent=1))
    return 0 if all(checks.values()) else 1


def _floor_check():
    import pandas as pd
    b = {"O": ["o1"], "P": ["p1"], "T": ["t1", "t2"], "target_reference": {"sd_ddof1_observed_design": {"t1": 0.0, "t2": 2.0}}}
    df = pd.DataFrame({"o1": [1.0, 2.0], "p1": [1.0, 3.0], "t1": [5.0, 5.0], "t2": [1.0, 3.0]})
    _, _, Y = block_arrays_floored(b, df)
    return bool(np.isfinite(Y).all() and np.allclose(Y[:, 0], 5.0 / (FLOOR * 2.0)))


# ---------------------------------------------------------------- run
def protocol_ready(path, want):
    got = sha_file(path)
    if got != want:
        raise SystemExit(f"protocol SHA-256 is {got}; the command named {want}. Refusing.")
    if "______" in open(path, encoding="utf-8").read():
        raise SystemExit("protocol has blank fields. Refusing.")
    return got


def cmd_run(a):
    from concurrent.futures import ProcessPoolExecutor
    g = gate3a()
    out = os.path.abspath(a.out)
    if os.path.exists(os.path.join(out, "n2_summary.json")):
        raise SystemExit("n2_summary.json exists: the run has already been made.")
    proto = protocol_ready(a.protocol, a.protocol_sha256)
    bench = os.path.abspath(a.bench)
    audit = g._bind(bench)
    studies = json.loads((audit.ROOT / "paper4_design_models/p4d_design_model_freeze_manifest.json").read_text())["studies"]
    os.makedirs(out, exist_ok=True)
    manifest = {"script": f"p6_n2.py v{VERSION}", "script_sha256": sha_file(os.path.abspath(__file__)),
                "protocol_sha256": proto, "gate3a_sha256": GATE3A_SHA, "audit_sha256": g.AUDIT_SHA,
                "predictor_sha256": g.PREDICTOR_SHA, "archive_root": str(audit.ROOT), "records": len(studies),
                "seed": SEED, "lambdas": LAMBDAS, "bandwidths": BWS, "lambda_u": LAMBDA_U, "claim_se": CLAIM_SE,
                "floor": FLOOR, "cap": CAP, "mde": MDE, "primary": PRIMARY, "evaluation_rows_used": False,
                "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(manifest, open(os.path.join(out, "n2_run_manifest.json"), "w"), indent=1)
    records, unusable = [], []
    with ProcessPoolExecutor(max_workers=a.workers) as pool, open(os.path.join(out, "n2_records.jsonl"), "w") as fh:
        for i, r in enumerate(pool.map(run_record, [(bench, s) for s in studies])):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            records.append(r)
            unusable += [dict(record=r["key"], **u) for u in r["unusable"]]
            print(f"[{i + 1}/{len(studies)}] {r['key']} blocks {r['blocks']} unusable {len(r['unusable'])}", flush=True)
    res = aggregate(records)
    summary = {**{k: v for k, v in manifest.items() if k != "started_utc"}, "records_with_usable_blocks":
               sum(1 for r in records if r["blocks"]), "unusable_blocks": len(unusable), "results": res,
               "primary_qualifies": res[f"{PRIMARY[0]}|{PRIMARY[1]}|sq"]["qualifies"],
               "finished_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    json.dump(summary, open(os.path.join(out, "n2_summary.json"), "w"), indent=1)
    json.dump(unusable, open(os.path.join(out, "n2_unusable.json"), "w"), indent=1)
    with open(os.path.join(out, "n2_records.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        cols = [f"{fam}_{loss}_{a}" for fam in FAMILIES for loss in ("sq", "bd") for a in ("G1",) + VARIANTS]
        w.writerow(["record", "blocks"] + cols)
        for r in records:
            w.writerow([r["key"], r["blocks"]] + [r["fam"][fam][loss][a] for fam in FAMILIES for loss in ("sq", "bd")
                                                   for a in ("G1",) + VARIANTS])
    for k, v in res.items():
        print(k, json.dumps(v if not isinstance(v, dict) else {kk: (round(vv, 5) if isinstance(vv, float) else vv)
                                                                for kk, vv in v.items()}))
    print("primary_qualifies", summary["primary_qualifies"])
    print("n2_summary.json", sha_file(os.path.join(out, "n2_summary.json")))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("selftest")
    st.add_argument("--planted-blocks", type=int, default=3)
    st.add_argument("--null-records", type=int, default=6)
    rn = sub.add_parser("run")
    rn.add_argument("--bench", required=True)
    rn.add_argument("--protocol", required=True)
    rn.add_argument("--protocol-sha256", required=True)
    rn.add_argument("--out", default="n2_out")
    rn.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    sys.exit(selftest(a) if a.cmd == "selftest" else cmd_run(a))


if __name__ == "__main__":
    main()
