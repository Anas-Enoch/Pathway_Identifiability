#!/usr/bin/env python3
"""
nhanes_fa_run.py  v0.1.0   NHANES serum fatty-acid acquisition test (NHANES_FA_PROTOCOL_v0.3.md).

Calls the FROZEN N2 implementation (nhanes_fa/frozen/p6_n2.py 681ce852..., p5_gate3a.py 4600e2f3...) for every
selection decision and every prediction. This file adds: the block allocator (protocol s3), an XPT loader restricted
to the protocol's variables (s2), a per-participant mirror of the frozen fold loop (so survey linearisation and the
base-panel arm are available; its pooled sums are tested equal to the frozen run_block), the design-based variance
(s5), the feasibility gate (s6) and the one-run guards (s9). It contains no selection logic of its own.

  python3 -I nhanes_fa_run.py selftest
  python3 -I nhanes_fa_run.py alloc  --analytes nhanes_fa_analytes.json --out nhanes_fa_block.json
  python3 -I nhanes_fa_run.py dev    --fas FAS_G.xpt --demo DEMO_G.xpt --block nhanes_fa_block.json --out nhanes_fa_out
  python3 -I nhanes_fa_run.py ext    --fas FAS_H.xpt --demo DEMO_H.xpt --block nhanes_fa_block.json --out nhanes_fa_out
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, os, sys, time
from pathlib import Path
import numpy as np

VERSION = "0.2.0"   # v0.2.0 (protocol v0.4): K prespecified blocks pooled at participant level; nondegenerate gate; selection code unchanged
PROTOCOL = "NHANES_FA_PROTOCOL_v0.4.md"
SALT_V04 = "NHANES-FA-v0.4|"
K_BLOCKS = 5
SALT = "NHANES-FA-v0.2|"
FROZEN_DIR = Path(__file__).resolve().parent / "nhanes_fa" / "frozen"
P6_SHA = "681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e"
P5_SHA = "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"
SEED_DEV, SEED_EXT = 20261103, 20261104
MIN_AGE, MIN_N = 12, 1500
MDE_CEILING = 0.0030          # 0.30 percentage points
SUCCESS_FACTOR = 0.10
ALPHA_ONE_SIDED = 0.025
POWER = 0.80
CLASS_ORDER = ("SFA", "MUFA", "PUFA")
ROLE_ORDER = ("O", "T", "P")
DEMO_VARS = ("SEQN", "RIDAGEYR", "SDMVSTRA", "SDMVPSU")
B_BOOT = 10_000


def sha_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------- frozen module
_N2 = None


def frozen():
    """Load p6_n2.py from nhanes_fa/frozen after verifying both hashes. p6_n2 itself re-verifies p5_gate3a."""
    global _N2
    if _N2 is None:
        p6, p5 = FROZEN_DIR / "p6_n2.py", FROZEN_DIR / "p5_gate3a.py"
        for p, want in ((p6, P6_SHA), (p5, P5_SHA)):
            if not p.exists() or sha_file(p) != want:
                raise SystemExit(f"REFUSED: {p} missing or sha256 != {want}")
        spec = importlib.util.spec_from_file_location("p6_n2_frozen", p6)
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        _N2 = m
    return _N2


# ---------------------------------------------------------------- s3 allocator
def allocate(analytes, salt=SALT):
    """analytes: list of {"var": LBX..., "class": SFA|MUFA|PUFA, "label": ...}. Returns block dict (protocol s3)."""
    by = {c: [] for c in CLASS_ORDER}
    for a in analytes:
        if a["class"] not in by:
            raise SystemExit(f"unknown class {a['class']} for {a['var']}")
        by[a["class"]].append(a["var"])
    seq = []
    for c in CLASS_ORDER:
        seq += [(c, v) for v in sorted(by[c], key=lambda v: (hashlib.sha256((salt + v).encode()).hexdigest(), v))]
    roles = {r: [] for r in ROLE_ORDER}
    for k, (c, v) in enumerate(seq):
        roles[ROLE_ORDER[k % 3]].append(v)
    m = len(seq)
    sizes = {"O": math.ceil(m / 3)}
    sizes["T"] = math.ceil((m - sizes["O"]) / 2); sizes["P"] = m - sizes["O"] - sizes["T"]
    assert {r: len(v) for r, v in roles.items()} == sizes, "allocator size invariant violated"
    cls = {v: c for c, v in seq}
    return {"salt": salt, "m": m, "sizes": sizes, "O": roles["O"], "P": roles["P"], "T": roles["T"],
            "class_of": cls, "sequence": [v for _, v in seq]}


# ---------------------------------------------------------------- s2 loader
def load_cycle(fas_path: Path, demo_path: Path, block, _frames=None):
    """_frames=(fas_df, demo_df) is used only by the self-test, which cannot write XPT files."""
    import pandas as pd
    need = ["SEQN", "WTFAS2YR"] + block["O"] + block["P"] + block["T"]
    if _frames is None:
        fas = pd.read_sas(fas_path, format="xport")
        demo = pd.read_sas(demo_path, format="xport")
    else:
        fas, demo = _frames
    missing = [v for v in need if v not in fas.columns] + [v for v in DEMO_VARS if v not in demo.columns]
    if missing:
        raise SystemExit(f"REFUSED: variables missing from files: {missing}")
    fas = fas[need]; demo = demo[list(DEMO_VARS)]
    for df, name in ((fas, "FAS"), (demo, "DEMO")):
        if df["SEQN"].duplicated().any():
            raise SystemExit(f"REFUSED: duplicate SEQN in {name}")
    d = fas.merge(demo, on="SEQN", how="left")
    if d[["SDMVSTRA", "SDMVPSU"]].isna().any().any():
        raise SystemExit("REFUSED: FAS records without DEMO design variables")
    an = block["O"] + block["P"] + block["T"]
    elig = (d["RIDAGEYR"] >= MIN_AGE) & (d["WTFAS2YR"] > 0) & d[an].notna().all(axis=1)
    d = d.assign(eligible=elig.to_numpy())
    return d


def arrays(d, block):
    e = d[d["eligible"]]
    O = e[block["O"]].to_numpy(float); X = e[block["P"]].to_numpy(float); T = e[block["T"]].to_numpy(float)
    n2 = frozen()
    sd = T.std(axis=0, ddof=1)
    sd = np.where(np.isfinite(sd) & (sd > 0), sd, 0.0)
    floor = n2.FLOOR * np.median(sd[sd > 0]) if (sd > 0).any() else 1.0
    Y = T / np.maximum(sd, floor)
    ids = e["SEQN"].to_numpy().astype(np.int64).astype(str).tolist()
    return O, X, Y, ids, {"target_sd": sd.tolist(), "floor": float(floor)}


# ---------------------------------------------------------------- mirror of the frozen fold loop (per-row losses)
def mirror_block(O, X, Y, subjects, seed):
    """Same calls and same decisions as p6_n2.run_block(family='rbf'); returns per-row losses for arms
    G1, NOVA-V, NOVA-U, NOVA-S and BASE (O alone), plus picks, claims and fold index per row."""
    n2 = frozen(); g = n2.gate3a()
    fold = g.folds_for(subjects, seed)
    n, m = len(Y), X.shape[1]
    arms = ("G1",) + n2.VARIANTS + ("BASE",)
    sq = {a: np.full(n, np.nan) for a in arms}; bd = {a: np.full(n, np.nan) for a in arms}
    picks, claims, folds = [], 0, 0
    for k in range(n2.K_OUT):
        tr, te = fold != k, fold == k
        if te.sum() == 0:
            continue
        jg = g.g1_pick(O[tr], Y[tr], X[tr])
        psi, best = n2.inner_scores("rbf", O[tr], X[tr], Y[tr], n2.SEED + k)
        ok = ~np.isnan(psi).any(axis=1)
        psi = psi[ok]
        Vh = psi.mean(axis=0)
        se = psi.std(axis=0, ddof=1) / np.sqrt(len(psi))
        jv = int(np.flatnonzero(Vh == Vh.max())[0])
        u = Vh - n2.LAMBDA_U * se
        ju = int(np.flatnonzero(u == u.max())[0])
        claim = False
        if jv != jg:
            dd = psi[:, jv] - psi[:, jg]
            sdd = dd.std(ddof=1) / np.sqrt(len(dd))
            claim = bool(dd.mean() > n2.CLAIM_SE * sdd)
        js = jv if claim else jg
        pk = {"G1": jg, "NOVA-V": jv, "NOVA-U": ju, "NOVA-S": js}
        picks.append(pk); claims += int(claim); folds += 1
        cache = {}
        for a, j in pk.items():
            if j not in cache:
                P = n2.predict_grid("rbf", O[tr], O[te], Y[tr], X[tr, j], X[te, j])[best[j + 1]]
                E2 = (P - Y[te]) ** 2
                cache[j] = (E2.sum(axis=1), np.minimum(E2, n2.CAP).sum(axis=1))
            sq[a][te], bd[a][te] = cache[j]
        P0 = n2.predict_grid("rbf", O[tr], O[te], Y[tr])[best[0]]
        E0 = (P0 - Y[te]) ** 2
        sq["BASE"][te], bd["BASE"][te] = E0.sum(axis=1), np.minimum(E0, n2.CAP).sum(axis=1)
    return {"sq": sq, "bd": bd, "fold": fold, "picks": picks, "claims": claims, "folds": folds,
            "cells": np.where(np.isnan(sq["G1"]), 0.0, float(Y.shape[1]))}


# ---------------------------------------------------------------- s5 design-based inference
def t_quantile(p, df):
    """Quantile of Student t with df degrees of freedom (bisection on the CDF via the regularized incomplete beta)."""
    def betacf(a, b, x):
        MAXIT, EPS, FPMIN = 300, 3e-14, 1e-300
        qab, qap, qam = a + b, a + 1, a - 1
        c, d = 1.0, 1.0 - qab * x / qap
        d = 1.0 / (d if abs(d) > FPMIN else FPMIN); h = d
        for mm in range(1, MAXIT + 1):
            m2 = 2 * mm
            aa = mm * (b - mm) * x / ((qam + m2) * (a + m2))
            d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > FPMIN else FPMIN)
            c = 1.0 + aa / c; c = c if abs(c) > FPMIN else FPMIN
            h *= d * c
            aa = -(a + mm) * (qab + mm) * x / ((a + m2) * (qap + m2))
            d = 1.0 + aa * d; d = 1.0 / (d if abs(d) > FPMIN else FPMIN)
            c = 1.0 + aa / c; c = c if abs(c) > FPMIN else FPMIN
            dl = d * c; h *= dl
            if abs(dl - 1.0) < EPS:
                break
        return h

    def betainc(a, b, x):
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log(1 - x))
        return bt * betacf(a, b, x) / a if x < (a + 1) / (a + b + 2) else 1.0 - bt * betacf(b, a, 1 - x) / b

    def cdf(t):
        x = df / (df + t * t)
        return 1.0 - 0.5 * betainc(df / 2.0, 0.5, x) if t >= 0 else 0.5 * betainc(df / 2.0, 0.5, x)
    lo, hi = -50.0, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) < p: lo = mid
        else: hi = mid
    return 0.5 * (lo + hi)


def ratio_linearised(num_i, den_i, w, strata, psu):
    """Delta = sum(w*num)/sum(w*den); Taylor-linearised stratified with-replacement PSU variance over the full file.
    Records not in the domain contribute num=den=0. Returns (Delta, SE, df, n_strata, n_psu)."""
    W_num, W_den = float((w * num_i).sum()), float((w * den_i).sum())
    delta = W_num / W_den
    u = w * (num_i - delta * den_i) / W_den
    var, n_psu, n_str = 0.0, 0, 0
    for h in np.unique(strata):
        mh = strata == h
        ps = np.unique(psu[mh]); nh = len(ps)
        n_str += 1; n_psu += nh
        if nh < 2:
            continue
        z = np.array([u[mh & (psu == p)].sum() for p in ps])
        var += nh / (nh - 1) * ((z - z.mean()) ** 2).sum()
    return delta, math.sqrt(var), n_psu - n_str, n_str, n_psu


def inference(d, res, w_col="WTFAS2YR"):
    """Per-participant a_i, c_i over the full file (zeros outside the domain); returns the s1/s5 quantities."""
    n_all = len(d)
    elig = d["eligible"].to_numpy()
    idx = np.flatnonzero(elig)
    out = {}
    w = d[w_col].to_numpy(float)
    strata = d["SDMVSTRA"].to_numpy(); psu = d["SDMVPSU"].to_numpy()
    for loss in ("sq", "bd"):
        L = res[loss]
        full = {a: np.zeros(n_all) for a in L}
        for a in L:
            v = np.where(np.isnan(L[a]), 0.0, L[a]); full[a][idx] = v
        # Delta_ext: NOVA-S vs G1 ;  Delta_G1->0: G1 vs BASE
        dS, seS, df, nstr, npsu = ratio_linearised(full["G1"] - full["NOVA-S"], full["G1"], w, strata, psu)
        dB, seB, _, _, _ = ratio_linearised(full["BASE"] - full["G1"], full["BASE"], w, strata, psu)
        tq = t_quantile(1 - ALPHA_ONE_SIDED, df)
        # unweighted participant bootstrap (sensitivity)
        rng = np.random.default_rng(SEED_DEV)
        a_i, c_i = (full["G1"] - full["NOVA-S"])[idx], full["G1"][idx]
        draws = np.empty(B_BOOT)
        for b in range(B_BOOT):
            ii = rng.integers(0, len(idx), len(idx)); draws[b] = a_i[ii].sum() / c_i[ii].sum()
        out[loss] = {"delta_ext": dS, "se": seS, "df": df, "lower_97.5_onesided": dS - tq * seS,
                     "delta_g1_base": dB, "se_g1_base": seB, "lower_97.5_g1_base": dB - tq * seB,
                     "unweighted_delta_ext": float(a_i.sum() / c_i.sum()),
                     "unweighted_boot_lower_97.5": float(np.quantile(draws, ALPHA_ONE_SIDED)),
                     "strata": nstr, "psu": npsu}
    out["n_file"] = int(n_all); out["n_eligible"] = int(elig.sum())
    out["claims"] = res["claims"]; out["folds"] = res["folds"]; out["picks"] = res["picks"]
    out["agreement_NOVA_S_G1"] = sum(p["NOVA-S"] == p["G1"] for p in res["picks"]) / max(1, res["folds"])
    return out


# ---------------------------------------------------------------- commands
def cmd_alloc(a):
    analytes = json.loads(Path(a.analytes).read_text())
    blk = allocate(analytes)
    Path(a.out).write_text(json.dumps(blk, indent=2) + "\n")
    print(json.dumps({k: blk[k] for k in ("m", "sizes", "O", "T", "P")}, indent=1)); print("sha256", sha_file(a.out))


def cmd_alloc5(a):
    """Protocol v0.4 s3: block 1 = the frozen v0.3 block (read from file, hash recorded); blocks 2..K from salts
    'NHANES-FA-v0.4|<k>|' with the same allocator. All K are written together and frozen together."""
    analytes = json.loads(Path(a.analytes).read_text())
    b1 = json.loads(Path(a.block1).read_text())
    if b1 != allocate(analytes, SALT):
        raise SystemExit("REFUSED: block 1 file does not equal the v0.3 allocation of this analyte list")
    blocks = [dict(b1, k=1, source="v0.3 block, sha256 " + sha_file(a.block1))]
    for k in range(2, K_BLOCKS + 1):
        blocks.append(dict(allocate(analytes, f"{SALT_V04}{k}|"), k=k))
    out = {"protocol": "v0.4", "K": K_BLOCKS, "blocks": blocks}
    Path(a.out).write_text(json.dumps(out, indent=2) + "\n")
    for b in blocks:
        print(f"block {b['k']}: O={b['O']}\n         T={b['T']}\n         P={b['P']}")
    print("sha256", sha_file(a.out))


def _run_cycle(a, which):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    tag = {"dev": "FAS_G_development", "ext": "FAS_H_evaluation"}[which]
    done = out / f"{tag}_v04_result.json"
    if done.exists():
        raise SystemExit(f"REFUSED: {done} exists; the protocol allows one run per cycle (a crashed run leaves no result file)")
    if which == "ext":
        feas = out / "FAS_G_development_v04_result.json"
        if not feas.exists() or json.loads(feas.read_text())["feasibility"]["decision"] != "GO":
            raise SystemExit("REFUSED: evaluation cycle requires a recorded GO in FAS_G_development_v04_result.json")
    bl = json.loads(Path(a.blocks).read_text())
    blocks = bl["blocks"] if "blocks" in bl else [dict(bl, k=1)]
    blk = blocks[0]
    hashes = {"fas": sha_file(a.fas), "demo": sha_file(a.demo), "blocks": sha_file(a.blocks), "script": sha_file(__file__),
              "p6_n2": P6_SHA, "p5_gate3a": P5_SHA, "protocol": sha_file(a.protocol) if Path(a.protocol).exists() else None}
    # eligibility uses the union of all analytes over all blocks (every block uses the same 29 analytes)
    allvars = sorted({v for b in blocks for r in ("O", "P", "T") for v in b[r]})
    d = load_cycle(Path(a.fas), Path(a.demo), {"O": allvars, "P": [], "T": []})
    n_el = int(d["eligible"].sum())
    ids = d[d["eligible"]]["SEQN"].to_numpy().astype(np.int64).astype(str).tolist()
    print(f"{tag}: file rows {len(d)}, eligible {n_el}")
    if which == "ext":
        dev_ids = set(json.loads((out / "FAS_G_development_v04_result.json").read_text())["seqn_sha256_list"])
        if any(hashlib.sha256(s.encode()).hexdigest() in dev_ids for s in ids):
            raise SystemExit("REFUSED: SEQN overlap between cycles")
    t0 = time.time()
    seed = SEED_DEV if which == "dev" else SEED_EXT
    per_block, stds, pooled = [], [], None
    for b in blocks:
        O, X, Y, ids_b, std = arrays(d, b)
        assert ids_b == ids
        r = mirror_block(O, X, Y, ids, seed)                # same subjects + seed -> identical folds in every block
        per_block.append({"k": b["k"], "inference": {kk: v for kk, v in inference(d, r).items() if kk != "picks"} | {"picks": r["picks"]}})
        stds.append(std)
        if pooled is None:
            pooled = {l: {a: np.where(np.isnan(r[l][a]), 0.0, r[l][a]) for a in r[l]} for l in ("sq", "bd")}
            pooled.update(fold=r["fold"], picks=list(r["picks"]), claims=r["claims"], folds=r["folds"], cells=r["cells"].copy())
        else:
            for l in ("sq", "bd"):
                for a_ in r[l]:
                    pooled[l][a_] = pooled[l][a_] + np.where(np.isnan(r[l][a_]), 0.0, r[l][a_])
            pooled["picks"] += r["picks"]; pooled["claims"] += r["claims"]; pooled["folds"] += r["folds"]; pooled["cells"] += r["cells"]
        print(f"  block {b['k']}: claims {r['claims']}/{r['folds']}, picks {[p['NOVA-S'] for p in r['picks']]} vs G1 {[p['G1'] for p in r['picks']]}")
    inf = inference(d, pooled)                               # pooled participant-level sums over blocks (s1, v0.4)
    inf["per_block"] = per_block
    inf["total_departures"] = pooled["claims"]
    rec = {"script": "nhanes_fa_run.py", "version": VERSION, "cycle": tag, "inputs": hashes, "K": len(blocks),
           "standardisation": stds, "seed": seed, "runtime_s": time.time() - t0,
           "inference": inf, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
           "seqn_sha256_list": [hashlib.sha256(s.encode()).hexdigest() for s in ids]}
    res = pooled
    if which == "dev":
        sq = inf["sq"]; df = sq["df"]
        mde = (t_quantile(1 - ALPHA_ONE_SIDED, df) + t_quantile(POWER, df)) * sq["se"]
        nondegenerate = pooled["claims"] >= 1 and math.isfinite(sq["se"]) and sq["se"] > 0
        go = (n_el >= MIN_N) and nondegenerate and (mde <= MDE_CEILING)
        rec["feasibility"] = {"n_eligible": n_el, "min_n": MIN_N, "departures": pooled["claims"], "block_folds": pooled["folds"],
                              "se_dev": sq["se"], "df": df, "mde_dev": mde, "mde_ceiling": MDE_CEILING,
                              "nondegenerate": bool(nondegenerate), "decision": "GO" if go else "NO-GO",
                              "note": "Delta_dev plays no part in this decision; a zero-departure or zero-SE development result is NO-GO (v0.4)"}
        print("FEASIBILITY", json.dumps(rec["feasibility"], indent=1))
    else:
        sq = inf["sq"]
        c1 = sq["lower_97.5_g1_base"] > 0
        c2 = sq["lower_97.5_onesided"] > 0
        c3 = sq["delta_ext"] >= SUCCESS_FACTOR * sq["delta_g1_base"]
        rec["success"] = {"baseline_lower_bound_positive": bool(c1), "delta_ext_lower_bound_positive": bool(c2),
                          "delta_ext_ge_factor_times_baseline": bool(c3), "factor": SUCCESS_FACTOR,
                          "verdict": "SUPPORTIVE SUCCESS" if (c1 and c2 and c3) else "FAILURE"}
        print("SUCCESS", json.dumps(rec["success"], indent=1))
    # per-participant table stays local
    import csv
    with (out / f"{tag}_v04_participants.csv").open("w", newline="") as fh:
        wr = csv.writer(fh); wr.writerow(["seqn_sha256", "fold"] + [f"sq_{a}" for a in res["sq"]] + [f"bd_{a}" for a in res["bd"]])
        for i, s in enumerate(ids):
            wr.writerow([hashlib.sha256(s.encode()).hexdigest(), int(res["fold"][i])] + [res["sq"][a][i] for a in res["sq"]] + [res["bd"][a][i] for a in res["bd"]])
    done.write_text(json.dumps(rec, indent=2, default=float) + "\n")
    print(json.dumps({k: inf["sq"][k] for k in ("delta_ext", "se", "df", "lower_97.5_onesided", "delta_g1_base", "lower_97.5_g1_base")}, indent=1))
    print("result", done, sha_file(done))


# ---------------------------------------------------------------- selftest
def cmd_selftest(a):
    ok = True

    def check(c, msg):
        nonlocal ok
        print(("PASS " if c else "FAIL ") + msg); ok &= bool(c)
    # allocator
    for counts in [(10, 6, 13), (7, 7, 7), (1, 1, 27), (10, 6, 12), (11, 6, 13)]:
        an = ([{"var": f"S{i}", "class": "SFA"} for i in range(counts[0])] + [{"var": f"M{i}", "class": "MUFA"} for i in range(counts[1])]
              + [{"var": f"P{i}", "class": "PUFA"} for i in range(counts[2])])
        b = allocate(an); b2 = allocate(list(reversed(an)))
        spread = max(max(sum(b["class_of"][v] == c for v in b[r]) for r in ROLE_ORDER) - min(sum(b["class_of"][v] == c for v in b[r]) for r in ROLE_ORDER) for c in CLASS_ORDER)
        check(b == b2 and spread <= 1 and sum(b["sizes"].values()) == sum(counts), f"allocator {counts}: deterministic, sizes {b['sizes']}, class spread {spread}")
    # equivalence: mirror vs frozen run_block on synthetic blocks
    n2 = frozen(); g = n2.gate3a()
    rng = np.random.default_rng(7)
    for kind in ("null", "quadratic"):
        O, X, Y = n2.synth(rng, 150, kind)
        subj = [f"s{i}" for i in range(len(Y))]
        ref = n2.run_block(O, X, Y, subj, SEED_DEV, "rbf")
        mir = mirror_block(O, X, Y, subj, SEED_DEV)
        same_sums = all(abs(np.nansum(mir[l][a]) - ref[l][a]) < 1e-9 for l in ("sq", "bd") for a in ("G1",) + n2.VARIANTS)
        same_picks = mir["picks"] == ref["picks"] and mir["claims"] == ref["claims"] and mir["folds"] == ref["folds"]
        same_n = abs(mir["cells"].sum() - ref["n"]) < 1e-9
        check(same_sums and same_picks and same_n, f"equivalence with frozen run_block ({kind}): sums, picks, claims, cells identical")
        check(np.isfinite(mir["sq"]["BASE"]).sum() == len(Y), f"base arm computed for every row ({kind})")
    # t quantile against known values
    check(abs(t_quantile(0.975, 10) - 2.2281) < 1e-3 and abs(t_quantile(0.80, 15) - 0.8662) < 1e-3 and abs(t_quantile(0.975, 1000) - 1.9623) < 1e-3,
          "t quantiles (df 10, 15, 1000) match tables")
    # linearised variance vs bootstrap under iid single-stratum design (2 PSUs per stratum, many strata)
    rng = np.random.default_rng(3)
    n = 4000
    num = rng.normal(0.01, 0.5, n); den = rng.lognormal(0, 0.3, n) + 1
    w = np.ones(n); strata = np.repeat(np.arange(40), 100); psu = np.tile(np.repeat([1, 2], 50), 40)
    dl, se, df, ns, npsu = ratio_linearised(num, den, w, strata, psu)
    draws = []
    for b in range(2000):
        # resample PSUs within strata
        idx = []
        for h in range(40):
            for p in rng.choice([1, 2], 2):
                idx.append(np.flatnonzero((strata == h) & (psu == p)))
        ii = np.concatenate(idx); draws.append(num[ii].sum() / den[ii].sum())
    check(0.6 < se / np.std(draws) < 1.6 and df == 40, f"linearised SE {se:.5f} vs PSU-bootstrap SD {np.std(draws):.5f}, df {df}")
    # domain estimation: zeros outside domain do not change Delta, strata/PSU retained
    num2 = num.copy(); den2 = den.copy(); num2[:500] = 0; den2[:500] = 0
    dl2, se2, df2, ns2, npsu2 = ratio_linearised(num2, den2, w, strata, psu)
    check(df2 == df and ns2 == ns and abs(dl2 - num[500:].sum() / den[500:].sum()) < 1e-12, "domain estimation keeps all strata/PSUs")
    # loader logic through injected frames
    import pandas as pd
    blk = allocate([{"var": f"A{i}", "class": c} for i, c in enumerate(["SFA"] * 4 + ["MUFA"] * 2 + ["PUFA"] * 3)])
    an = blk["O"] + blk["P"] + blk["T"]
    fas = pd.DataFrame({"SEQN": [1, 2, 3, 4, 5], "WTFAS2YR": [1.0, 0.0, 2.0, 3.0, 1.0], **{v: [1.0, 1.0, 1.0, np.nan, 1.0] for v in an}})
    demo = pd.DataFrame({"SEQN": [1, 2, 3, 4, 5], "RIDAGEYR": [30, 40, 11, 50, 12], "SDMVSTRA": [1, 1, 2, 2, 2], "SDMVPSU": [1, 2, 1, 2, 2]})
    d = load_cycle(None, None, blk, _frames=(fas, demo))
    check(d["eligible"].tolist() == [True, False, False, False, True], "eligibility: age>=12, weight>0, complete analytes")
    try:
        load_cycle(None, None, blk, _frames=(pd.concat([fas, fas.iloc[[0]]]), demo)); check(False, "duplicate SEQN refused")
    except SystemExit as e:
        check("duplicate SEQN" in str(e), "duplicate SEQN refused")
    try:
        load_cycle(None, None, blk, _frames=(fas.drop(columns=[an[0]]), demo)); check(False, "missing variable refused")
    except SystemExit as e:
        check("variables missing" in str(e), "missing variable refused")
    # v0.4: five salts give five distinct allocations, deterministic, sizes exact
    an29 = ([{"var": f"S{i}", "class": "SFA"} for i in range(10)] + [{"var": f"M{i}", "class": "MUFA"} for i in range(6)]
            + [{"var": f"P{i}", "class": "PUFA"} for i in range(13)])
    b5 = [allocate(an29, SALT)] + [allocate(an29, f"{SALT_V04}{k}|") for k in range(2, K_BLOCKS + 1)]
    sets = [tuple(sorted(b["O"])) for b in b5]
    check(len(set(sets)) == K_BLOCKS and all(b["sizes"] == {"O": 10, "T": 10, "P": 9} for b in b5), "v0.4: 5 distinct deterministic allocations, sizes 10/10/9")
    # pooled inference: two identical blocks give the same Delta as one, with cells doubled
    rng = np.random.default_rng(11)
    O, X, Y = n2.synth(rng, 150, "quadratic"); subj = [f"s{i}" for i in range(len(Y))]
    r1 = mirror_block(O, X, Y, subj, SEED_DEV)
    pooled = {l: {a_: 2 * np.where(np.isnan(r1[l][a_]), 0.0, r1[l][a_]) for a_ in r1[l]} for l in ("sq", "bd")}
    pooled.update(fold=r1["fold"], picks=r1["picks"] * 2, claims=2 * r1["claims"], folds=2 * r1["folds"], cells=2 * r1["cells"])
    import pandas as pd
    dd = pd.DataFrame({"SEQN": range(150), "WTFAS2YR": 1.0, "SDMVSTRA": np.repeat(np.arange(15), 10), "SDMVPSU": np.tile([1, 2], 75), "eligible": True})
    i1, i2 = inference(dd, r1), inference(dd, pooled)
    check(abs(i1["sq"]["delta_ext"] - i2["sq"]["delta_ext"]) < 1e-12 and abs(i1["sq"]["se"] - i2["sq"]["se"]) < 1e-12, "pooling two identical blocks leaves Delta and SE unchanged")
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    al = sub.add_parser("alloc"); al.add_argument("--analytes", required=True); al.add_argument("--out", default="nhanes_fa_block.json")
    a5 = sub.add_parser("alloc5"); a5.add_argument("--analytes", required=True); a5.add_argument("--block1", required=True); a5.add_argument("--out", default="nhanes_fa_blocks_v0.4.json")
    for name in ("dev", "ext"):
        r = sub.add_parser(name)
        r.add_argument("--fas", required=True); r.add_argument("--demo", required=True); r.add_argument("--blocks", required=True)
        r.add_argument("--protocol", default=PROTOCOL); r.add_argument("--out", default="nhanes_fa_out")
    a = ap.parse_args()
    if a.cmd == "selftest": cmd_selftest(a)
    elif a.cmd == "alloc": cmd_alloc(a)
    elif a.cmd == "alloc5": cmd_alloc5(a)
    else: _run_cycle(a, a.cmd)


if __name__ == "__main__":
    main()
