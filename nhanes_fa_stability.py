#!/usr/bin/env python3
"""
nhanes_fa_stability.py  v0.1.0   Prospective plug-in prediction of G1 refit stability on NHANES 2013-2014
(NHANES_FA_STABILITY_PREDICTION_PROTOCOL_v0.4.md).

Every eligibility, fold and selection decision is made by the FROZEN code (nhanes_fa_run.py f85d5235..., which loads
nhanes_fa/frozen/p6_n2.py 681ce852... and p5_gate3a.py 4600e2f3...), verified by hash before use. This file adds:
`g1_scores` (the score vector behind the frozen `g1_pick`, tested for argmax equivalence including ties),
`novas_picks_only` (the selection logic of the frozen `run_block` with no held-out loss evaluation, tested for pick
equivalence), the plug-in models M_G and M_E fitted to 2011-2012, the exact-size prediction table, the in-sample
screens, the freeze, the guarded single external opening and the comparison. It adds no selection logic of its own.

  python3 -I nhanes_fa_stability.py selftest [--dev]            synthetic tests; --dev adds the development-data tests
  python3 -I nhanes_fa_stability.py fit --fas FAS_G.xpt --demo DEMO_G.xpt
  python3 -I nhanes_fa_stability.py bench                       time one cycle; extrapolate the full table
  python3 -I nhanes_fa_stability.py predict [--sizes a:b] [--workers k]   resumable; parts in nhanes_fa_stab_out/parts
  python3 -I nhanes_fa_stability.py assemble                    parts -> nhanes_fa_prediction_v0.4.json
  python3 -I nhanes_fa_stability.py screen                      protocol s3 (a)-(d); writes the assessable-size set
  python3 -I nhanes_fa_stability.py freeze                      hashes -> NHANES_FA_STABILITY_FREEZE_HASHES.txt
  python3 -I nhanes_fa_stability.py run-ext --fas FAS_H.xpt --demo DEMO_H.xpt   LOCKED: needs the authorisation file
  python3 -I nhanes_fa_stability.py compare
"""
from __future__ import annotations
import argparse, datetime, hashlib, importlib.util, json, math, os, sys, time
from pathlib import Path
import numpy as np

VERSION = "0.1.0"
PROTOCOL = "NHANES_FA_STABILITY_PREDICTION_PROTOCOL_v0.4.md"
SALT = "NHANES-FA-STAB-v0.4|"
HERE = Path(__file__).resolve().parent
WRAPPER, WRAPPER_SHA = HERE / "nhanes_fa_run.py", "f85d5235cd2fd609145be9ba2b444cd2c4b394a8b3ed6c761c942e765a730297"
P6_SHA = "681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e"
P5_SHA = "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"
BLOCKS, BLOCKS_SHA = HERE / "nhanes_fa_blocks_v0.4.json", "e3e1d8ac636afd9690e4b5b00b884e21ccf79dad77993e4a75d8776ef3a14338"
ANALYTES = HERE / "nhanes_fa_analytes.json"
DEV_RESULT = HERE / "nhanes_fa_out" / "FAS_G_development_v04_result.json"
SEED_EXT = 20261104                      # frozen 2013-2014 fold seed (v0.4 s4)
N_LO, N_HI, N_DEV = 1200, 2400, 1725     # verdict range and development size
B, B_S = 2000, 100                       # cycles per (model, size); NOVA-S picks-only cycles at N_DEV
U_CAP = 5                                # sharpness requirement: U <= 5 of 50 pairs
MODELS = ("MG", "ME")
OUT = HERE / "nhanes_fa_stab_out"
PARTS = OUT / "parts"
FIT = OUT / "dev_fit.json"                 # aggregate fit (mean, cov, picks, pairs): committable
ROWS = OUT / "dev_rows.npy"                # participant-level eligible matrix for M_E: stays local
IDS = OUT / "dev_ids.json"                 # SEQN identifiers of the eligible rows, same order: stays local
PRED = HERE / "nhanes_fa_prediction_v0.4.json"
SCREEN = OUT / "screen_v0.4.json"
FREEZE = HERE / "NHANES_FA_STABILITY_FREEZE_HASHES.txt"
AUTH = HERE / "NHANES_FA_STABILITY_AUTHORISATION.txt"
ATTEMPT = OUT / "FAS_H_opening_attempt.json"
EXT_RESULT = OUT / "FAS_H_stability_result.json"
ATTAINABLE = (0, 4, 6, 7, 8, 9, 10)


def sha_file(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def utc(): return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
def seed_for(model, n): return int(hashlib.sha256(f"{SALT}{model}|{n}".encode()).hexdigest()[:15], 16)


# ---------------------------------------------------------------- frozen code
_W = None


def wrapper():
    global _W
    if _W is None:
        for p, want in ((WRAPPER, WRAPPER_SHA), (HERE / "nhanes_fa" / "frozen" / "p6_n2.py", P6_SHA),
                        (HERE / "nhanes_fa" / "frozen" / "p5_gate3a.py", P5_SHA), (BLOCKS, BLOCKS_SHA)):
            if not p.exists() or sha_file(p) != want:
                raise SystemExit(f"REFUSED: {p.name} missing or sha256 != {want}")
        spec = importlib.util.spec_from_file_location("nhanes_fa_run_frozen", WRAPPER)
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
        _W = m
    return _W


def n2(): return wrapper().frozen()
def gate(): return n2().gate3a()
def blocks(): return json.loads(BLOCKS.read_text())["blocks"]
def allvars(): return sorted({v for b in blocks() for r in ("O", "P", "T") for v in b[r]})


# ---------------------------------------------------------------- G1 score vector (mirror of p5_gate3a.g1_pick)
def g1_scores(O, T, P):
    """Same standardisation, shrinkage, ridge and clip as the frozen g1_pick; returns the full score vector.
    argmax (lowest index) is tested equal to g1_pick in selftest (ii)."""
    Z = np.column_stack([O, T, P]).astype(float)
    mu, sd = Z.mean(axis=0), Z.std(axis=0, ddof=1)
    sd = np.where(sd > 1e-12, sd, 1.0)
    Z = (Z - mu) / sd
    n, d = Z.shape
    S = np.cov(Z, rowvar=False, ddof=1)
    alpha = min(1.0, max(0.05, d / (d + n)))
    S = (1 - alpha) * S + alpha * (np.trace(S) / d) * np.eye(d)
    o, t = O.shape[1], T.shape[1]
    io, it, ip = slice(0, o), slice(o, o + t), slice(o + t, d)
    Soo = S[io, io] + 1e-6 * np.eye(o)
    sol_p = np.linalg.solve(Soo, S[io, ip])
    C = S[it, ip] - S[it, io] @ sol_p
    v = np.diag(S[ip, ip] - S[ip, io] @ sol_p)
    return (C * C).sum(axis=0) / np.maximum(v, 1e-12)


def lowest_argmax(s): return int(np.flatnonzero(s == s.max())[0])


def runner_up(s, j):
    rest = np.delete(np.arange(len(s)), j)
    m = s[rest].max()
    return int(rest[np.flatnonzero(s[rest] == m)[0]])


# ---------------------------------------------------------------- NOVA-S picks only (mirror of p6_n2.run_block selection)
def novas_picks_only(O, X, Y, subjects, seed, family="rbf"):
    """Picks of G1, NOVA-V, NOVA-U, NOVA-S per outer fold exactly as run_block decides them; no held-out loss."""
    m = n2(); g = gate()
    fold = g.folds_for(subjects, seed)
    picks = []
    for k in range(m.K_OUT):
        tr, te = fold != k, fold == k
        if te.sum() == 0:
            continue
        jg = g.g1_pick(O[tr], Y[tr], X[tr])
        psi, best = m.inner_scores(family, O[tr], X[tr], Y[tr], m.SEED + k)
        ok = ~np.isnan(psi).any(axis=1)
        psi = psi[ok]
        Vh = psi.mean(axis=0)
        se = psi.std(axis=0, ddof=1) / np.sqrt(len(psi))
        jv = int(np.flatnonzero(Vh == Vh.max())[0])
        u = Vh - m.LAMBDA_U * se
        ju = int(np.flatnonzero(u == u.max())[0])
        claim = False
        if jv != jg:
            d = psi[:, jv] - psi[:, jg]
            sd = d.std(ddof=1) / np.sqrt(len(d))
            claim = bool(d.mean() > m.CLAIM_SE * sd)
        js = jv if claim else jg
        picks.append({"G1": jg, "NOVA-V": jv, "NOVA-U": ju, "NOVA-S": js})
    return picks


# ---------------------------------------------------------------- the statistic on one cycle
def block_cols(av):
    idx = {v: i for i, v in enumerate(av)}
    return [{r: np.array([idx[v] for v in b[r]]) for r in ("O", "P", "T")} for b in blocks()]


def disagreement_count(picks):
    return int(sum(1 for i in range(len(picks)) for j in range(i + 1, len(picks)) if picks[i] != picks[j]))


def cycle_stats(rows, cols, pairs, seed=SEED_EXT, with_ds=True, ids=None):
    """rows: (n, 29) in allvars order; identifiers 0..n-1 as strings unless given. Returns per-block picks, counts c_b,
    C, and signed differences on the fixed pairs."""
    g = gate()
    n = rows.shape[0]
    ids = [str(i) for i in range(n)] if ids is None else ids
    fold = g.folds_for(ids, seed)
    picks, cs, ds = [], [], []
    for b, c in enumerate(cols):
        pk, dd = [], []
        for k in range(n2().K_OUT):
            tr = fold != k
            O, P, T = rows[tr][:, c["O"]], rows[tr][:, c["P"]], rows[tr][:, c["T"]]
            s = g1_scores(O, T, P)
            pk.append(lowest_argmax(s))
            if with_ds:
                js, jr = pairs[b]
                dd.append(float(s[js] - s[jr]))
        picks.append(pk); cs.append(disagreement_count(pk)); ds.append(dd)
    return {"picks": picks, "c": cs, "C": int(sum(cs)), "ds": ds}


def frozen_cycle(d, seed, pairs=None):
    """Picks on a loaded NHANES frame exactly as the frozen wrapper makes them: `arrays` (v0.4 target scaling),
    `folds_for` on SEQN ids, `g1_pick(O, Y, X)` per outer fold. g1_scores is computed alongside and its lowest-index
    argmax is required to equal g1_pick on every training set; signed gaps use the fixed pairs when given."""
    w = wrapper(); g = gate(); picks, cs, ds = [], [], []
    for b, blk in enumerate(blocks()):
        O, X, Y, ids, _ = w.arrays(d, blk)
        fold = g.folds_for(ids, seed); pk, dd = [], []
        for k in range(n2().K_OUT):
            tr = fold != k
            jg = g.g1_pick(O[tr], Y[tr], X[tr]); s = g1_scores(O[tr], Y[tr], X[tr])
            if lowest_argmax(s) != jg:
                raise SystemExit(f"REFUSED: g1_scores argmax != frozen g1_pick (block {b+1}, fold {k})")
            pk.append(jg)
            if pairs is not None:
                js, jr = pairs[b]; dd.append(float(s[js] - s[jr]))
        picks.append(pk); cs.append(disagreement_count(pk)); ds.append(dd)
    return {"picks": picks, "c": cs, "C": int(sum(cs)), "ds": ds, "n_eligible": int(d["eligible"].sum())}


# ---------------------------------------------------------------- fit (2011-2012 only)
def cmd_fit(a):
    w = wrapper(); OUT.mkdir(exist_ok=True)
    av = allvars()
    d = w.load_cycle(Path(a.fas), Path(a.demo), {"O": av, "P": [], "T": []})
    e = d[d["eligible"]]
    rows = e[av].to_numpy(float)
    n = rows.shape[0]
    cols = block_cols(av)
    # full-fit picks, pairs and recorded-pick check (protocol s2.3, s3(d), selftest (i))
    g = gate(); pairs, full_scores, full_picks = [], [], []
    for c in cols:
        s = g1_scores(rows[:, c["O"]], rows[:, c["T"]], rows[:, c["P"]])
        j = lowest_argmax(s); jr = runner_up(s, j)
        if j != g.g1_pick(rows[:, c["O"]], rows[:, c["T"]], rows[:, c["P"]]):
            raise SystemExit("REFUSED: g1_scores argmax differs from frozen g1_pick on the full development fit")
        pairs.append([j, jr]); full_scores.append(s.tolist()); full_picks.append(j)
    dev = frozen_cycle(d, w.SEED_DEV, pairs)   # frozen path with the development seed must reproduce the recorded 25 picks
    rec = json.loads(DEV_RESULT.read_text())["inference"]["per_block"]
    recorded = [[p["G1"] for p in blk["inference"]["picks"]] for blk in rec]
    if dev["picks"] != recorded:
        raise SystemExit(f"REFUSED: development G1 picks {dev['picks']} != recorded {recorded}")
    fit = {"script": Path(__file__).name, "version": VERSION, "protocol": PROTOCOL, "created_utc": utc(),
           "inputs": {"fas": sha_file(a.fas), "demo": sha_file(a.demo), "blocks": BLOCKS_SHA},
           "n_eligible": int(n), "allvars": av, "mean": rows.mean(axis=0).tolist(), "cov_ddof1": np.cov(rows, rowvar=False, ddof=1).tolist(),
           "full_fit_picks": full_picks, "fixed_pairs": pairs, "full_fit_scores": full_scores,
           "full_fit_signed_gap": [float(full_scores[b][pairs[b][0]] - full_scores[b][pairs[b][1]]) for b in range(5)],
           "dev_cycle_picks": dev["picks"], "dev_c": dev["c"], "dev_C": dev["C"], "dev_ds": dev["ds"], "recorded_picks_match": True,
           "pick_names": [[blocks()[b]["P"][j] for j in pairs[b]] for b in range(5)]}
    FIT.write_text(json.dumps(fit, indent=1)); np.save(ROWS, rows)
    IDS.write_text(json.dumps(e["SEQN"].to_numpy().astype(np.int64).astype(str).tolist()))
    print(json.dumps({k: fit[k] for k in ("n_eligible", "full_fit_picks", "fixed_pairs", "pick_names", "full_fit_signed_gap", "dev_C")}, indent=1))
    print(f"fit written: {FIT} (aggregate, committable); {ROWS}, {IDS} (participant-level, stay local)")


def load_fit():
    if not FIT.exists() or not ROWS.exists():
        raise SystemExit("REFUSED: run `fit` first (dev_fit.json and dev_rows.npy required)")
    return json.loads(FIT.read_text()), np.load(ROWS)


# ---------------------------------------------------------------- prediction table
def draw(model, rng, n, fit, rows):
    if model == "MG":
        return rng.multivariate_normal(np.array(fit["mean"]), np.array(fit["cov_ddof1"]), size=n, method="cholesky")
    return rows[rng.integers(0, rows.shape[0], n)]


def quantile_set(hist):
    """Frozen discrete convention: L = min{c: P(C<=c) >= .025}, U = min{c: P(C<=c) >= .975}, Monte Carlo estimates."""
    p = np.asarray(hist, float) / float(np.sum(hist)); cdf = np.cumsum(p)
    L = int(np.flatnonzero(cdf >= 0.025)[0]); U = int(np.flatnonzero(cdf >= 0.975)[0])
    return L, U


def simulate_size(model, n, fit, rows, Bn, keep_ds=False):
    cols = block_cols(fit["allvars"]); pairs = fit["fixed_pairs"]
    rng = np.random.default_rng(seed_for(model, n))
    hist = np.zeros(51, int); pickc = np.zeros((5, 9), int); ds_neg = np.zeros(5, int); ds_n = 0
    ds_sum = np.zeros(5); ds_sq = np.zeros(5); ds_keep = [] if keep_ds else None
    t0 = time.time()
    for _ in range(Bn):
        st = cycle_stats(draw(model, rng, n, fit, rows), cols, pairs)
        hist[st["C"]] += 1
        for b in range(5):
            for j in st["picks"][b]: pickc[b, j] += 1
            dd = np.array(st["ds"][b]); ds_neg[b] += int((dd < 0).sum()); ds_sum[b] += dd.sum(); ds_sq[b] += (dd ** 2).sum()
        ds_n += 5
        if keep_ds: ds_keep.append(st["ds"])
    L, U = quantile_set(hist)
    out = {"model": model, "n": n, "B": Bn, "seed": seed_for(model, n), "hist_C": hist.tolist(), "L": L, "U": U,
           "P_C0": float(hist[0] / Bn), "pick_counts": pickc.tolist(), "modal_pick": [int(np.argmax(pickc[b])) for b in range(5)],
           "modal_prob": [float(pickc[b].max() / pickc[b].sum()) for b in range(5)],
           "ds_reversal_prob": (ds_neg / ds_n).tolist(), "ds_mean": (ds_sum / ds_n).tolist(),
           "ds_sd": np.sqrt(np.maximum(ds_sq / ds_n - (ds_sum / ds_n) ** 2, 0)).tolist(), "seconds": round(time.time() - t0, 1)}
    if keep_ds: out["ds_samples"] = ds_keep
    return out


def _part(args):
    model, n, Bn = args
    fit, rows = load_fit()
    PARTS.mkdir(parents=True, exist_ok=True)
    p = PARTS / f"{model}_{n}.json"
    if p.exists():
        return str(p)
    res = simulate_size(model, n, fit, rows, Bn, keep_ds=(n == N_DEV))
    tmp = p.with_suffix(".tmp"); tmp.write_text(json.dumps(res)); os.replace(tmp, p)
    return str(p)


def cmd_bench(a):
    fit, rows = load_fit()
    for model in MODELS:
        t0 = time.time(); simulate_size(model, N_DEV, fit, rows, 10); dt = (time.time() - t0) / 10
        total = dt * B * (N_HI - N_LO + 1)
        print(f"{model}: {dt*1000:.0f} ms per cycle at n={N_DEV}; full table for this model ~ {total/3600:.1f} core-hours")


def cmd_predict(a):
    lo, hi = (int(x) for x in a.sizes.split(":")) if a.sizes else (N_LO, N_HI)
    sizes = sorted(set(range(lo, hi + 1)) | ({N_DEV} if lo <= N_DEV <= hi or not a.sizes else set()))
    tasks = [(m, n, a.B) for m in MODELS for n in sizes if not (PARTS / f"{m}_{n}.json").exists()]
    print(f"{len(tasks)} (model, size) parts to compute with B={a.B}, workers={a.workers}")
    if a.workers > 1:
        import multiprocessing as mp
        with mp.Pool(a.workers) as pool:
            for i, p in enumerate(pool.imap_unordered(_part, tasks), 1):
                if i % 20 == 0: print(f"  {i}/{len(tasks)} done", flush=True)
    else:
        for i, t in enumerate(tasks, 1):
            _part(t)
            if i % 20 == 0: print(f"  {i}/{len(tasks)} done", flush=True)
    print("predict: parts complete; run `assemble`")


def cmd_assemble(a):
    fit, _ = load_fit()
    table = {m: {} for m in MODELS}; missing = []
    for m in MODELS:
        for n in list(range(N_LO, N_HI + 1)) + [N_DEV]:
            p = PARTS / f"{m}_{n}.json"
            if not p.exists(): missing.append((m, n)); continue
            r = json.loads(p.read_text())
            if r["B"] != B: raise SystemExit(f"REFUSED: part {p.name} has B={r['B']} != frozen B={B}")
            table[m][str(n)] = {k: v for k, v in r.items() if k != "ds_samples"}
            if n == N_DEV: table[m][str(n)]["ds_samples_sha256"] = hashlib.sha256(json.dumps(r.get("ds_samples")).encode()).hexdigest()
    complete = not missing
    pred = {"script": Path(__file__).name, "version": VERSION, "protocol": PROTOCOL, "salt": SALT, "created_utc": utc(),
            "fit_sha256": sha_file(FIT), "blocks_sha256": BLOCKS_SHA, "B": B, "verdict_range": [N_LO, N_HI], "n_dev": N_DEV,
            "quantile_convention": "L=min{c:P(C<=c)>=.025}, U=min{c:P(C<=c)>=.975}; Monte Carlo estimates from B cycles",
            "complete": complete, "missing_parts": missing[:50], "table": table}
    PRED.write_text(json.dumps(pred))
    print(f"assembled: complete={complete}, missing={len(missing)}; written {PRED}")
    if not complete:
        print("NOTE: an incomplete table is descriptive only and cannot support a freeze or an external opening (protocol s2)")


# ---------------------------------------------------------------- screens (protocol s3)
def cmd_screen(a):
    fit, _ = load_fit(); pred = json.loads(PRED.read_text())
    if not pred["complete"]: raise SystemExit("REFUSED: prediction table incomplete; screens require the full table")
    mg = pred["table"]["MG"]; d = mg[str(N_DEV)]
    res = {"created_utc": utc(), "prediction_sha256": sha_file(PRED)}
    res["a_dev_C_in_set"] = {"L": d["L"], "U": d["U"], "dev_C": fit["dev_C"], "pass": d["L"] <= fit["dev_C"] <= d["U"]}
    assessable = sorted(int(n) for n in mg if N_LO <= int(n) <= N_HI and mg[n]["U"] <= U_CAP)
    res["b_sharpness"] = {"U_cap": U_CAP, "assessable_sizes": assessable, "n_assessable": len(assessable),
                          "all_sizes_pass": len(assessable) == (N_HI - N_LO + 1), "any_size_passes": len(assessable) > 0,
                          "U_by_size_summary": {"min": min(mg[n]["U"] for n in mg), "max": max(mg[n]["U"] for n in mg)}}
    res["c_modal_picks"] = {"predicted": d["modal_pick"], "development": fit["full_fit_picks"], "pass": d["modal_pick"] == fit["full_fit_picks"]}
    res["d_full_fit_gap_positive"] = {"gaps": fit["full_fit_signed_gap"], "pass": all(g > 0 for g in fit["full_fit_signed_gap"])}
    res["ME_at_dev"] = {k: pred["table"]["ME"][str(N_DEV)][k] for k in ("L", "U", "P_C0", "modal_pick")}
    res["freeze_permitted"] = bool(res["a_dev_C_in_set"]["pass"] and res["b_sharpness"]["any_size_passes"] and res["c_modal_picks"]["pass"] and res["d_full_fit_gap_positive"]["pass"])
    OUT.mkdir(exist_ok=True); SCREEN.write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "b_sharpness"} | {"b_sharpness": {k: v for k, v in res["b_sharpness"].items() if k != "assessable_sizes"}}, indent=1))
    print("SCREEN", "PASS (freeze permitted; assessable sizes stored)" if res["freeze_permitted"] else "FAIL (protocol not to be frozen)")


# ---------------------------------------------------------------- freeze / authorisation / external run / compare
def frozen_files():
    return [HERE / PROTOCOL, Path(__file__).resolve(), PRED, FIT, SCREEN, WRAPPER, BLOCKS, ANALYTES,
            HERE / "nhanes_fa" / "frozen" / "p6_n2.py", HERE / "nhanes_fa" / "frozen" / "p5_gate3a.py"]


def cmd_freeze(a):
    if FREEZE.exists(): raise SystemExit("REFUSED: freeze record exists; it is never overwritten")
    sc = json.loads(SCREEN.read_text())
    if not sc["freeze_permitted"]: raise SystemExit("REFUSED: screen did not permit freezing")
    if sc["prediction_sha256"] != sha_file(PRED): raise SystemExit("REFUSED: prediction file changed since the screen")
    lines = [f"# {PROTOCOL} freeze, {utc()}; assessable sizes: {len(sc['b_sharpness']['assessable_sizes'])} integers in [{N_LO},{N_HI}]"]
    lines += [f"{sha_file(p)}  {p.relative_to(HERE)}" for p in frozen_files()]
    FREEZE.write_text("\n".join(lines) + "\n"); print(FREEZE.read_text())


def check_frozen():
    if not FREEZE.exists(): raise SystemExit("REFUSED: no freeze record")
    want = {l.split("  ", 1)[1]: l.split("  ", 1)[0] for l in FREEZE.read_text().splitlines() if l and not l.startswith("#")}
    for rel, h in want.items():
        if sha_file(HERE / rel) != h: raise SystemExit(f"REFUSED: {rel} changed since freeze")
    return want


def cmd_run_ext(a):
    """Single external opening. Order: authorisation -> attempt record -> read -> hash -> disjointness -> compute."""
    want = check_frozen()
    if not AUTH.exists(): raise SystemExit("REFUSED: no authorisation file")
    auth = AUTH.read_text()
    for rel in (PROTOCOL, Path(__file__).name, PRED.name):
        if want[rel] not in auth: raise SystemExit(f"REFUSED: authorisation does not cite the frozen hash of {rel}")
    if "AUTHORISED" not in auth: raise SystemExit("REFUSED: authorisation file lacks the word AUTHORISED")
    settings = {"fas": str(a.fas), "demo": str(a.demo), "seed": SEED_EXT, "protocol_sha256": want[PROTOCOL], "script_sha256": want[Path(__file__).name]}
    OUT.mkdir(exist_ok=True)
    att = json.loads(ATTEMPT.read_text()) if ATTEMPT.exists() else {"attempts": []}
    if any(x["status"] == "completed" for x in att["attempts"]): raise SystemExit("REFUSED: the single external opening has been completed")
    if att["attempts"] and att["attempts"][-1]["settings"] != settings: raise SystemExit("REFUSED: retry must use identical inputs and settings")
    if EXT_RESULT.exists(): raise SystemExit("REFUSED: external result exists")
    att["attempts"].append({"status": "started", "utc": utc(), "settings": settings}); ATTEMPT.write_text(json.dumps(att, indent=1))
    try:   # ---- from here on the external files are read; any failure is recorded on the attempt before re-raising
        w = wrapper(); fit, _ = load_fit(); av = fit["allvars"]
        hashes = {"fas": sha_file(a.fas), "demo": sha_file(a.demo)}
        att["attempts"][-1]["input_sha256"] = hashes; ATTEMPT.write_text(json.dumps(att, indent=1))
        if len(att["attempts"]) > 1 and att["attempts"][-2].get("input_sha256") not in (None, hashes):
            raise SystemExit("REFUSED: retry inputs differ from the first attempt")
        d = w.load_cycle(Path(a.fas), Path(a.demo), {"O": av, "P": [], "T": []})
        e = d[d["eligible"]]; n = int(len(e))
        ids = e["SEQN"].to_numpy().astype(np.int64).astype(str).tolist()
        dev_ids = set(json.loads(DEV_RESULT.read_text())["seqn_sha256_list"])
        if any(hashlib.sha256(s.encode()).hexdigest() in dev_ids for s in ids): raise SystemExit("REFUSED: SEQN overlap between cycles")
        fc = frozen_cycle(d, SEED_EXT, fit["fixed_pairs"]); picks, cs, ds = fc["picks"], fc["c"], fc["ds"]
        # S3 picks-only NOVA path (descriptive), with the v0.4 target scaling as in the frozen wrapper
        s3 = []
        for b in blocks():
            O, X, Y, subj, _ = w.arrays(d, b)
            s3.append(novas_picks_only(O, X, Y, subj, SEED_EXT))
        res = {"script": Path(__file__).name, "version": VERSION, "protocol": PROTOCOL, "created_utc": utc(), "inputs": hashes,
               "n_file": int(len(d)), "n_eligible": n, "assessable": N_LO <= n <= N_HI, "picks": picks, "c": cs, "C": int(sum(cs)),
               "ds": ds, "ds_reversals": int(sum(1 for b in ds for x in b if x < 0)), "g1_pick_argmax_check": True,
               "S3_picks": s3, "S3_agreement_NOVA_S_G1": float(np.mean([p["NOVA-S"] == p["G1"] for blk in s3 for p in blk]))}
        EXT_RESULT.write_text(json.dumps(res, indent=1))
    except BaseException as ex:
        att["attempts"][-1]["status"] = "failed"; att["attempts"][-1]["error"] = f"{type(ex).__name__}: {ex}"[:500]
        att["attempts"][-1]["failed_utc"] = utc(); ATTEMPT.write_text(json.dumps(att, indent=1)); raise
    att["attempts"][-1]["status"] = "completed"; att["attempts"][-1]["result_sha256"] = sha_file(EXT_RESULT); ATTEMPT.write_text(json.dumps(att, indent=1))
    print(json.dumps({k: res[k] for k in ("n_eligible", "assessable", "c", "C", "ds_reversals", "S3_agreement_NOVA_S_G1")}, indent=1))


def cmd_compare(a):
    check_frozen(); pred = json.loads(PRED.read_text()); res = json.loads(EXT_RESULT.read_text()); sc = json.loads(SCREEN.read_text())
    n = res["n_eligible"]; out = {"n_eligible": n, "C": res["C"], "c": res["c"]}
    if not (N_LO <= n <= N_HI): out["verdict"] = "not assessable (outside verdict range)"
    elif n not in sc["b_sharpness"]["assessable_sizes"]: out["verdict"] = "not assessable (sharpness requirement failed at this size)"
    else:
        for m in MODELS:
            t = pred["table"][m][str(n)]; hist = np.array(t["hist_C"]); Bn = hist.sum()
            out[m] = {"L": t["L"], "U": t["U"], "P_C_ge_obs": float(hist[res["C"]:].sum() / Bn), "P_C_le_obs": float(hist[:res["C"] + 1].sum() / Bn),
                      "verdict": "compatible" if t["L"] <= res["C"] <= t["U"] else ("above" if res["C"] > t["U"] else "below")}
        out["verdict"] = f"MG: {out['MG']['verdict']}; ME (secondary): {out['ME']['verdict']}"
    (OUT / "comparison.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out, indent=1))


# ---------------------------------------------------------------- self-tests
def cmd_selftest(a):
    fails = []
    def check(c, msg):
        print(("PASS " if c else "FAIL ") + msg)
        if not c: fails.append(msg)
    g = gate(); m = n2(); rng = np.random.default_rng(11)
    # (ii) argmax equivalence incl. exact ties
    bad = 0
    for trial in range(30):
        n = int(rng.integers(60, 400)); O = rng.normal(size=(n, 10)); P = rng.normal(size=(n, 9))
        T = O[:, :5] @ rng.normal(size=(5, 10)) + 0.3 * P[:, [trial % 9]] + rng.normal(size=(n, 10))
        bad += lowest_argmax(g1_scores(O, T, P)) != g.g1_pick(O, T, P)
    check(bad == 0, f"(ii) argmax(g1_scores) == g1_pick on 30 random cases ({bad} mismatches)")
    n = 300; O = rng.normal(size=(n, 10)); z = rng.normal(size=n)
    P = np.column_stack([rng.normal(size=(n, 3)), z, z, rng.normal(size=(n, 4))]); T = np.column_stack([z + 0.1 * rng.normal(size=n)] * 10)
    s = g1_scores(O, T, P)
    check(abs(s[3] - s[4]) < 1e-12 and lowest_argmax(s) == 3 == g.g1_pick(O, T, P), "(ii) exact duplicate candidates: tie resolved to the lowest index by both")
    T2 = T * 1e-3 + 7.0; check(lowest_argmax(g1_scores(O, T2, P)) == g.g1_pick(O, T2, P) == 3, "(ii) argmax equivalence preserved under target rescaling")
    # (iii) fold pattern
    ids = [str(i) for i in range(1725)]; f = g.folds_for(ids, SEED_EXT)
    check(sorted(np.bincount(f, minlength=5).tolist()) == [345, 345, 345, 345, 345], f"(iii) folds_for on synthetic ids: sizes {np.bincount(f).tolist()}")
    # (iv) attainable per-block counts, exhaustive over 9^5 patterns
    import itertools
    seen = sorted({disagreement_count(p) for p in itertools.product(range(9), repeat=5)})
    check(seen == list(ATTAINABLE), f"(iv) attainable c_b = {seen}")
    # (v) planted tie vs separated: synthetic Gaussian 29-column cycles
    av = allvars(); cols = block_cols(av); d = len(av)
    def synth_fit(gap):
        mean = np.zeros(d); A = rng.normal(size=(d, d)) * 0.15; cov = A @ A.T + np.eye(d)
        return {"allvars": av, "mean": mean.tolist(), "cov_ddof1": cov.tolist(), "fixed_pairs": [[0, 1]] * 5}
    fit_sep = synth_fit(1.0)
    # separated: make two P columns of block 1 strongly/weakly related to T -> stable pick; tie: duplicate column structure
    r = simulate_size("MG", 400, fit_sep, None, 20)
    check(sum(r["hist_C"]) == 20 and all(0 <= x <= 50 for x in range(51)), "(v) simulate_size runs; histogram sums to B")
    # a planted exact population tie between two candidates must produce disagreement in a substantial share of cycles
    cov = np.array(fit_sep["cov_ddof1"]); c0 = cols[0]; i, j = c0["P"][0], c0["P"][1]
    cov[j, :] = cov[i, :]; cov[:, j] = cov[:, i]; cov[j, j] = cov[i, i] + 1e-9  # candidate j duplicates candidate i (population tie)
    for t_idx in c0["T"]:
        cov[t_idx, i] = cov[i, t_idx] = 0.8 * math.sqrt(cov[i, i] * cov[t_idx, t_idx]); cov[t_idx, j] = cov[j, t_idx] = cov[t_idx, i]
    w_, V = np.linalg.eigh(cov); cov = (V * np.maximum(w_, 1e-6)) @ V.T
    fit_tie = dict(fit_sep, cov_ddof1=cov.tolist())
    r_tie = simulate_size("MG", 400, fit_tie, None, 30)
    share = 1 - r_tie["hist_C"][0] / 30
    check(share > 0.5, f"(v) planted population tie in block 1 -> cycles with C>0: {share:.2f} (expected large)")
    # (vii) picks-only path == frozen run_block picks
    rng2 = np.random.default_rng(7)
    for kind in ("null", "quadratic"):
        O, X, Y = m.synth(rng2, 150, kind); subj = [f"s{i}" for i in range(len(Y))]
        ref = m.run_block(O, X, Y, subj, 20261103, "rbf")
        check(novas_picks_only(O, X, Y, subj, 20261103) == ref["picks"], f"(vii) novas_picks_only picks == frozen run_block picks ({kind})")
    # quantile convention
    h = np.zeros(51, int); h[0] = 1900; h[4] = 60; h[6] = 40
    check(quantile_set(h) == (0, 4), f"(conv) quantile set on a test histogram = {quantile_set(h)} (expected (0, 4))")
    if a.dev:
        fit, rows = load_fit(); w = wrapper()
        check(fit["recorded_picks_match"] and fit["dev_C"] == 0, "(i) development picks reproduce the recorded 25 picks; C=0")
        # (i) the simulation path (raw rows, standardised inside g1_scores) and the frozen path (v0.4-scaled Y) pick identically
        same = True
        for c in cols:
            fold = g.folds_for([str(i) for i in range(rows.shape[0])], 20261103)
            for k in range(5):
                tr = fold != k; T = rows[tr][:, c["T"]]; sd = T.std(axis=0, ddof=1); Ys = T / np.maximum(sd, m.FLOOR * np.median(sd))
                same &= lowest_argmax(g1_scores(rows[tr][:, c["O"]], T, rows[tr][:, c["P"]])) == g.g1_pick(rows[tr][:, c["O"]], Ys, rows[tr][:, c["P"]])
        check(same, "(i) raw-row g1_scores argmax == frozen g1_pick on v0.4-scaled targets, 25 training sets")
        check(all(x > 0 for x in fit["full_fit_signed_gap"]), f"(viii) full-fit signed gaps positive: {['%.3g' % x for x in fit['full_fit_signed_gap']]}")
        # (vi) M_E identity resample reproduces the development picks with the development seed
        st = cycle_stats(rows, cols, fit["fixed_pairs"], seed=w.SEED_DEV, ids=json.loads(IDS.read_text()))
        check(st["picks"] == fit["dev_cycle_picks"] and st["C"] == 0, "(vi) simulation path on the identity resample with SEQN ids and the development seed reproduces the 25 picks, C=0")
    print("SELFTEST", "PASS" if not fails else f"FAIL ({len(fails)})"); return 0 if not fails else 1


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("selftest"); s.add_argument("--dev", action="store_true")
    s = sub.add_parser("fit"); s.add_argument("--fas", required=True); s.add_argument("--demo", required=True)
    sub.add_parser("bench")
    s = sub.add_parser("predict"); s.add_argument("--sizes", default=None); s.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1)); s.add_argument("--B", type=int, default=B)
    sub.add_parser("assemble"); sub.add_parser("screen"); sub.add_parser("freeze")
    s = sub.add_parser("run-ext"); s.add_argument("--fas", required=True); s.add_argument("--demo", required=True)
    sub.add_parser("compare")
    a = ap.parse_args()
    return {"selftest": cmd_selftest, "fit": cmd_fit, "bench": cmd_bench, "predict": cmd_predict, "assemble": cmd_assemble,
            "screen": cmd_screen, "freeze": cmd_freeze, "run-ext": cmd_run_ext, "compare": cmd_compare}[a.cmd](a) or 0


if __name__ == "__main__":
    sys.exit(main())
