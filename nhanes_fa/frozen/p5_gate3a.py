#!/usr/bin/env python3
"""Gate 3A: does better conditional prediction reduce acquisition regret against G1? (development data only)

Implements GATE3A_PROTOCOL_v1.0.md. Entry points:
    python3 -I p5_gate3a.py selftest
        the four synthetic self-tests of protocol section 6.3
    python3 -I p5_gate3a.py run --bench <uapa_run/paper4_predictor_benchmark> --protocol <file>
                                --protocol-sha256 <hash> [--out gate3a_out] [--workers 4]
        the real run; refuses unless the protocol hash matches and its section 9 is filled; design rows only,
        through the archived F1 loader (audit.load_record, SHA-256 checked), which asserts each record's
        design and evaluation ID-set hashes
    python3 -I p5_gate3a.py secondary --bench <...> [--out gate3a_out]
        the evaluation-half secondary of protocol section 4; refuses before gate3a_summary.json exists

Estimation (protocol section 3), per block and outer fold k (5 folds over design subjects, seed 20261019):
  on D_-k only: inner 5-fold CV (seed 20261019 + k) chooses each family's hyperparameters per input set and
  gives V_f(j) = L_inner(O) - L_inner(O u {j}); j_f = argmax_j V_f(j) (ties: lowest index); j_G1 is
  re-estimated on D_-k with the frozen Paper IV shrinkage rule. Every model is refitted on D_-k and scored
  once on D_k. Loss per cell: min(e^2, 25). Records: losses pooled over blocks, folds and cells; reported
  numbers are equal-record means with a record bootstrap (10,000 resamples, seed 20261019).
"""
import argparse
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

VERSION = "0.1.0"
SEED = 20261019
CAP = 25.0
LAMBDAS = np.array([0.01, 0.1, 1.0, 10.0])
BANDWIDTHS = (0.5, 1.0, 2.0)
FAMILIES = ("P0", "P1", "P2", "P3")
TESTED = ("P1", "P2", "P3")
MDE = 0.0029
ALPHA_BONF = 0.05 / 3
K_OUT, K_IN = 5, 5
B_REPS = 10_000


# ---------------------------------------------------------------- transforms and design matrices
class LogRule:
    """Protocol section 2: per feature, fitted on training rows only."""

    def fit(self, X):
        self.neg = (X < 0).any(axis=0)
        pos = np.where(X > 0, X, np.inf).min(axis=0)
        self.drop = (~self.neg) & ~np.isfinite(pos)
        self.c = np.where(np.isfinite(pos), 0.5 * pos, 1.0)
        return self

    def transform(self, X):
        out = np.array(X, dtype=float, copy=True)
        lg = ~self.neg & ~self.drop
        out[:, lg] = np.log(np.maximum(out[:, lg], 0.0) + self.c[lg])
        return out[:, ~self.drop]


def _standardize(tr, te):
    mu = tr.mean(axis=0)
    sd = tr.std(axis=0, ddof=1) if tr.shape[0] > 1 else np.ones(tr.shape[1])
    keep = sd > 1e-12
    return (tr[:, keep] - mu[keep]) / sd[keep], (te[:, keep] - mu[keep]) / sd[keep]


def features(family, Otr, Ote, xtr=None, xte=None):
    """Returns (train, test) inputs for family f with base panel O and optional candidate column x."""
    Atr = Otr if xtr is None else np.column_stack([Otr, xtr])
    Ate = Ote if xte is None else np.column_stack([Ote, xte])
    if family == "P0":
        return _standardize(Atr, Ate)
    rule = LogRule().fit(Atr)
    Ltr, Lte = rule.transform(Atr), rule.transform(Ate)
    Ltr, Lte = _standardize(Ltr, Lte)
    if family in ("P1", "P3"):
        return Ltr, Lte
    # P2: [u, u^2] for the base panel; with a candidate m: [u, u^2, m, m^2, u*m]
    n_base = Otr.shape[1] - int(rule.drop[:Otr.shape[1]].sum())
    has_m = xtr is not None and not rule.drop[-1] and Ltr.shape[1] > n_base
    if has_m:
        u_tr, m_tr = Ltr[:, :-1], Ltr[:, -1:]
        u_te, m_te = Lte[:, :-1], Lte[:, -1:]
        Ftr = np.column_stack([u_tr, u_tr ** 2, m_tr, m_tr ** 2, u_tr * m_tr])
        Fte = np.column_stack([u_te, u_te ** 2, m_te, m_te ** 2, u_te * m_te])
    else:
        Ftr = np.column_stack([Ltr, Ltr ** 2])
        Fte = np.column_stack([Lte, Lte ** 2])
    return _standardize(Ftr, Fte)


def fit_predict(family, Xtr, Ytr, Xte):
    """Predictions on Xte for every hyperparameter of the family: array (n_hyper, n_te, n_targets)."""
    my = Ytr.mean(axis=0)
    Yc = Ytr - my
    if family != "P3":
        if Xtr.shape[1] == 0:
            return np.repeat(my[None, None, :], len(LAMBDAS), axis=0).repeat(Xte.shape[0], axis=1)
        U, s, Vt = np.linalg.svd(Xtr, full_matrices=False)
        UtY = U.T @ Yc
        out = []
        for lam in LAMBDAS:
            B = Vt.T @ ((s / (s * s + lam))[:, None] * UtY)
            out.append(Xte @ B + my)
        return np.stack(out)
    d2 = ((Xtr[:, None, :] - Xtr[None, :, :]) ** 2).sum(-1)
    med = np.sqrt(np.median(d2[np.triu_indices(len(Xtr), 1)])) if len(Xtr) > 1 else 1.0
    med = med if med > 1e-12 else 1.0
    d2te = ((Xte[:, None, :] - Xtr[None, :, :]) ** 2).sum(-1)
    out = []
    for bw in BANDWIDTHS:
        s2 = 2.0 * (bw * med) ** 2
        K = np.exp(-d2 / s2)
        Kte = np.exp(-d2te / s2)
        w, Q = np.linalg.eigh(K)
        QtY = Q.T @ Yc
        for lam in LAMBDAS:
            alpha = Q @ (QtY / (np.maximum(w, 0.0) + lam)[:, None])
            out.append(Kte @ alpha + my)
    return np.stack(out)


def n_hyper(family):
    return len(LAMBDAS) * (len(BANDWIDTHS) if family == "P3" else 1)


def bounded_sum(P, Y):
    """P: (n_hyper, n, t) predictions; returns summed bounded loss per hyper."""
    E = P - Y[None]
    return np.minimum(E * E, CAP).sum(axis=(1, 2))


# ---------------------------------------------------------------- G1 (frozen Paper IV rule)
def g1_pick(O, T, P):
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
    sol_t = np.linalg.solve(Soo, S[io, it])
    sol_p = np.linalg.solve(Soo, S[io, ip])
    C = S[it, ip] - S[it, io] @ sol_p
    v = np.diag(S[ip, ip] - S[ip, io] @ sol_p)
    score = (C * C).sum(axis=0) / np.maximum(v, 1e-12)
    return int(np.argmax(score))


# ---------------------------------------------------------------- folds
def folds_for(subjects, seed):
    uniq = sorted(set(subjects))
    perm = np.random.default_rng(seed).permutation(len(uniq))
    f = {s: int(perm[i] % K_OUT) for i, s in enumerate(uniq)}
    return np.array([f[s] for s in subjects])


def inner_select(family, O, X, Y, seed):
    """Inner CV on the given rows: chosen hyper and inner loss for O and for each candidate."""
    n = len(Y)
    fold = np.random.default_rng(seed).permutation(n) % K_IN
    m = X.shape[1]
    tot = np.zeros((m + 1, n_hyper(family)))
    for f in range(K_IN):
        tr, te = fold != f, fold == f
        if te.sum() == 0 or tr.sum() < 3:
            continue
        for j in range(-1, m):
            xtr = None if j < 0 else X[tr, j]
            xte = None if j < 0 else X[te, j]
            Ftr, Fte = features(family, O[tr], O[te], xtr, xte)
            tot[j + 1] += bounded_sum(fit_predict(family, Ftr, Y[tr], Fte), Y[te])
    best = tot.argmin(axis=1)
    loss = tot[np.arange(m + 1), best] / (n * Y.shape[1])
    return best, loss


# ---------------------------------------------------------------- one block
def run_block(O, X, Y, subjects, record_seed, leaky=False):
    """Returns per-family outer-fold loss sums and selections for one block."""
    fold = folds_for(subjects, record_seed)
    m = X.shape[1]
    res = {f: {"S_O": np.zeros(K_OUT), "S_j": np.zeros((K_OUT, m)), "jhat": [], "n_cells": np.zeros(K_OUT)}
           for f in FAMILIES}
    g1 = []
    for k in range(K_OUT):
        tr, te = fold != k, fold == k
        if te.sum() == 0:
            g1.append(-1)
            continue
        jg = g1_pick(O[tr], Y[tr], X[tr])
        g1.append(jg)
        for f in FAMILIES:
            best, loss = inner_select(f, O[tr], X[tr], Y[tr], SEED + k)
            V = loss[0] - loss[1:]
            jhat = int(np.flatnonzero(V == V.max())[0])
            r = res[f]
            r["jhat"].append(jhat)
            r["n_cells"][k] = te.sum() * Y.shape[1]
            Ftr, Fte = features(f, O[tr], O[te])
            r["S_O"][k] = bounded_sum(fit_predict(f, Ftr, Y[tr], Fte)[[best[0]]], Y[te])[0]
            for j in range(m):
                Ftr, Fte = features(f, O[tr], O[te], X[tr, j], X[te, j])
                r["S_j"][k, j] = bounded_sum(fit_predict(f, Ftr, Y[tr], Fte)[[best[j + 1]]], Y[te])[0]
    return res, g1


def summarize_block(res, g1):
    """Per family: record-level sums for each arm plus diagnostics inputs."""
    out = {}
    ks = [k for k in range(K_OUT) if g1[k] >= 0]
    for f, r in res.items():
        jh = r["jhat"]
        s_hat = sum(r["S_j"][k, jh[i]] for i, k in enumerate(ks))
        s_g1 = sum(r["S_j"][k, g1[k]] for k in ks)
        s_leak = sum(r["S_j"][k].min() for k in ks)
        s_o = r["S_O"][ks].sum()
        n = r["n_cells"][ks].sum()
        Sj = r["S_j"][ks].sum(axis=0)
        odd = [k for k in ks if k % 2 == 0]
        even = [k for k in ks if k % 2 == 1]
        va = r["S_O"][odd].sum() - r["S_j"][odd].sum(axis=0)
        vb = r["S_O"][even].sum() - r["S_j"][even].sum(axis=0)
        out[f] = {"hat": s_hat, "g1": s_g1, "leak": s_leak, "O": s_o, "n": n,
                  "hind_gain": s_o - Sj.min(), "agree": sum(jh[i] == g1[k] for i, k in enumerate(ks)),
                  "nfold": len(ks), "split_rho": spearman(va, vb)}
    return out


def spearman(a, b):
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    ra, rb = np.argsort(np.argsort(a)), np.argsort(np.argsort(b))
    return float(np.corrcoef(ra, rb)[0, 1])


# ---------------------------------------------------------------- aggregation and inference
def boot_lower(vals, q, seed=SEED, reps=B_REPS):
    v = np.asarray(vals, dtype=float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(v), size=(reps, len(v)))
    means = v[idx].mean(axis=1)
    return float(np.quantile(means, q))


def aggregate(records):
    """records: {record: {family: summed dict over blocks}} -> summary dict (aggregates only)."""
    fam = {f: {"acq": [], "pred": [], "dep": [], "leak": [], "A": [], "C": [], "agree": 0, "nfold": 0}
           for f in FAMILIES}
    rhos = {f: [] for f in FAMILIES}
    for rec, d in records.items():
        for f in FAMILIES:
            x, p1, p0 = d[f], d["P1"], d["P0"]
            R = lambda s, n=x["n"]: s / n
            fam[f]["acq"].append(1 - R(x["hat"]) / R(x["g1"]))
            fam[f]["pred"].append(1 - R(x["O"]) / (p1["O"] / p1["n"]))
            fam[f]["dep"].append(1 - R(x["hat"]) / (p0["g1"] / p0["n"]))
            fam[f]["leak"].append(1 - R(x["leak"]) / R(x["g1"]))
            fam[f]["A"].append(x["hind_gain"] / x["O"])
            w_f = R(x["O"]) - R(x["hat"])
            w_1 = p1["O"] / p1["n"] - p1["hat"] / p1["n"]
            fam[f]["C"].append((w_f - w_1) / (p1["O"] / p1["n"]))
            fam[f]["agree"] += x["agree"]
            fam[f]["nfold"] += x["nfold"]
            rhos[f].extend(d[f]["rhos"])
    out = {}
    for f, v in fam.items():
        acq = np.array(v["acq"])
        row = {"delta_acq": float(acq.mean()), "delta_acq_lower_98.33": boot_lower(acq, ALPHA_BONF),
               "delta_pred_vs_P1": float(np.mean(v["pred"])),
               "delta_pred_upper95": boot_lower(v["pred"], 0.95),
               "delta_dep": float(np.mean(v["dep"])), "delta_acq_leaky_selection": float(np.mean(v["leak"])),
               "diag_A_hindsight": float(np.mean(v["A"])), "diag_C_vs_P1": float(np.mean(v["C"])),
               "diag_C_lower_98.33": boot_lower(v["C"], ALPHA_BONF) if f in ("P2", "P3") else None,
               "agreement_with_G1": v["agree"] / max(v["nfold"], 1),
               "diag_B_split_half_spearman_median": float(np.nanmedian(rhos[f])) if rhos[f] else None}
        row["qualifies"] = bool(f in TESTED and row["delta_acq"] >= MDE and row["delta_acq_lower_98.33"] > 0)
        out[f] = row
    return out


def accumulate(records, rec, blk_summary):
    d = records.setdefault(rec, {f: {"hat": 0.0, "g1": 0.0, "leak": 0.0, "O": 0.0, "n": 0.0, "hind_gain": 0.0,
                                     "agree": 0, "nfold": 0, "rhos": []} for f in FAMILIES})
    for f, x in blk_summary.items():
        for k in ("hat", "g1", "leak", "O", "n", "hind_gain", "agree", "nfold"):
            d[f][k] += x[k]
        d[f]["rhos"].append(x["split_rho"])


# ---------------------------------------------------------------- synthetic data and self-tests
def synth_block(rng, n, kind):
    """kind: 'planted' (one candidate with value only through z^2, one linear decoy) or 'null'."""
    O = rng.lognormal(0.0, 0.6, size=(n, 20))
    lo = np.log(O)
    W = rng.normal(0, 0.25, size=(20, 20))
    T = lo @ W
    X = rng.lognormal(0.0, 0.6, size=(n, 20))
    if kind == "planted":
        z = rng.normal(size=n)
        X[:, 7] = z                                   # raw values with negatives: identity branch of the log rule
        T[:, :10] += 1.2 * (z[:, None] ** 2 - 1.0)    # value only through z^2
        w = np.log(X[:, 3])
        T[:, :10] += 0.35 * w[:, None]                # linear decoy that G1 and the linear families prefer
    T += rng.normal(0, 1.0, size=T.shape)
    Y = (T - T.mean(0)) / T.std(0, ddof=1)
    return O, X, Y


def selftest(args):
    t0 = time.time()
    report = {}
    # 1. log-rule domain
    Xt = np.array([[-1.0, 0.0, 2.0, 0.0], [3.0, 0.0, 0.5, 1.0], [0.5, 0.0, 0.0, 4.0]])
    lr = LogRule().fit(Xt)
    Xo = lr.transform(np.array([[-5.0, 1.0, -3.0, -2.0]]))
    ok1 = (bool(lr.neg[0]) and bool(lr.drop[1]) and not lr.neg[2] and Xo.shape == (1, 3)
           and np.isfinite(Xo).all() and abs(Xo[0, 0] + 5.0) < 1e-12
           and abs(Xo[0, 1] - math.log(0.25)) < 1e-12)
    ok1 = bool(ok1)
    report["log_rule_domain"] = ok1
    # 2. planted nonlinear candidate
    rng = np.random.default_rng(1)
    picks = {f: [] for f in FAMILIES}
    g1s = []
    for b in range(args.planted_blocks):
        O, X, Y = synth_block(rng, 160, "planted")
        res, g1 = run_block(O, X, Y, list(range(160)), SEED + b)
        g1s += g1
        for f in FAMILIES:
            picks[f] += res[f]["jhat"]
    rate = {f: float(np.mean(np.array(picks[f]) == 7)) for f in FAMILIES}
    report["planted_selection_rate"] = rate
    report["planted_G1_rate"] = float(np.mean(np.array(g1s) == 7))
    ok2 = bool(rate["P2"] >= 0.8 and rate["P3"] >= 0.5 and rate["P1"] <= 0.2 and rate["P0"] <= 0.2)
    report["planted_ok"] = ok2
    # 3 and 4. null replicates: qualifications and leakage
    quals, honest, leaky = 0, [], []
    for r in range(args.null_reps):
        rng = np.random.default_rng(100 + r)
        records = {}
        for rec in range(args.null_records):
            O, X, Y = synth_block(rng, 50, "null")
            res, g1 = run_block(O, X, Y, list(range(50)), SEED + rec)
            accumulate(records, f"R{rec}", summarize_block(res, g1))
        s = aggregate(records)
        quals += sum(s[f]["qualifies"] for f in TESTED)
        honest.append(np.mean([s[f]["delta_acq"] for f in TESTED]))
        leaky.append(np.mean([s[f]["delta_acq_leaky_selection"] for f in TESTED]))
    se = lambda v: float(np.std(v, ddof=1) / math.sqrt(len(v)))
    report["null_qualifications"] = [int(quals), 3 * args.null_reps]
    report["null_honest_mean_delta_acq"] = [float(np.mean(honest)), se(honest)]
    report["null_leaky_mean_delta_acq"] = [float(np.mean(leaky)), se(leaky)]
    ok3 = bool(quals <= 3)
    ok4 = bool((np.mean(honest) <= 2 * se(honest)) and (np.mean(leaky) > 2 * se(leaky)))
    report["null_ok"], report["leakage_ok"] = bool(ok3), bool(ok4)
    report["all_pass"] = bool(ok1 and ok2 and ok3 and ok4)
    report["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(report, indent=1))
    return 0 if report["all_pass"] else 1


# ---------------------------------------------------------------- archive binding (development data only)
AUDIT_SHA = "2d0128a599720293a852be4c0766fab8268db5aa31fd5ffb65fc3a8eb95642bd"
PREDICTOR_SHA = "9a591bf7be15733cf7801de1a548aab5544cfa72d0af20ef26066cbc0447dfe2"
_AUDIT = None


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _bind(bench):
    """Import the archived F1 loader after checking it is byte-identical to the one F1 ran."""
    global _AUDIT
    if _AUDIT is None:
        for f, want in (("audit.py", AUDIT_SHA), ("predictor.py", PREDICTOR_SHA)):
            got = sha_file(os.path.join(bench, f))
            if got != want:
                raise SystemExit(f"{f} in {bench} has SHA-256 {got}, expected {want}")
        sys.path.insert(0, bench)
        import audit  # noqa: E402  (archived loader; asserts the design and evaluation ID-set hashes)
        _AUDIT = audit
    return _AUDIT


def block_arrays(b, frame):
    O = frame[b["O"]].to_numpy(float)
    X = frame[b["P"]].to_numpy(float)
    T = frame[b["T"]].to_numpy(float)
    sd = np.array([b["target_reference"]["sd_ddof1_observed_design"][t] for t in b["T"]], float)
    sd = np.where(sd > 0, sd, 1.0)
    return O, X, T / sd


def log_branches(A):
    neg = (A < 0).any(axis=0)
    nopos = ~(A > 0).any(axis=0)
    return {"log": int((~neg & ~nopos).sum()), "identity_negative": int(neg.sum()),
            "dropped_no_positive": int((~neg & nopos).sum())}


def run_record(task):
    bench, study = task
    audit = _bind(bench)
    key, m, _df, dev, _test = audit.load_record(study)
    del _df, _test  # evaluation rows are never used in the primary run
    out = {"key": key, "blocks": [], "unusable": [], "branches": {"log": 0, "identity_negative": 0,
                                                                  "dropped_no_positive": 0}}
    ids = dev.index.to_numpy().astype(str)
    for b in m["blocks"]:
        O, X, Y = block_arrays(b, dev)
        if not (np.isfinite(O).all() and np.isfinite(X).all() and np.isfinite(Y).all()):
            out["unusable"].append({"block": b["block_id"], "reason": "non-finite design value"})
            continue
        for k, v in log_branches(np.column_stack([O, X])).items():
            out["branches"][k] += v
        try:
            res, g1 = run_block(O, X, Y, list(ids), SEED)
        except (ValueError, np.linalg.LinAlgError) as e:
            out["unusable"].append({"block": b["block_id"], "reason": str(e)[:120]})
            continue
        sb = summarize_block(res, g1)
        g1_all = g1_pick(O, Y, X)
        sb["_g1_full_design_agreement"] = sum(int(j == g1_all) for j in g1 if j >= 0)
        sb["_g1_folds"] = sum(int(j >= 0) for j in g1)
        out["blocks"].append(sb)
    return out


def protocol_ready(path, want_sha):
    got = sha_file(path)
    if got != want_sha:
        raise SystemExit(f"protocol SHA-256 is {got}; the command named {want_sha}. Refusing to run.")
    text = open(path, encoding="utf-8").read()
    if "______" in text:
        raise SystemExit("protocol section 9 still has blank fields (OSF registration). Refusing to run.")
    return got


def cmd_run(a):
    from concurrent.futures import ProcessPoolExecutor
    out = os.path.abspath(a.out)
    if os.path.exists(os.path.join(out, "gate3a_summary.json")):
        raise SystemExit("gate3a_summary.json exists: the real run has already been made.")
    proto = protocol_ready(a.protocol, a.protocol_sha256)
    bench = os.path.abspath(a.bench)
    audit = _bind(bench)
    freeze = json.loads((audit.ROOT / "paper4_design_models/p4d_design_model_freeze_manifest.json").read_text())
    studies = freeze["studies"]
    os.makedirs(out, exist_ok=True)
    manifest = {"script": f"p5_gate3a.py v{VERSION}", "script_sha256": sha_file(os.path.abspath(__file__)),
                "protocol_sha256": proto, "audit_sha256": AUDIT_SHA, "predictor_sha256": PREDICTOR_SHA,
                "archive_root": str(audit.ROOT), "records": len(studies), "seed": SEED, "outer_k": K_OUT,
                "inner_k": K_IN, "lambdas": LAMBDAS.tolist(), "bandwidths": list(BANDWIDTHS), "cap": CAP,
                "evaluation_rows_read": False, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    with open(os.path.join(out, "gate3a_run_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1)
    records, unusable, branches, g1agree = {}, [], {"log": 0, "identity_negative": 0, "dropped_no_positive": 0}, [0, 0]
    with ProcessPoolExecutor(max_workers=a.workers) as pool, open(os.path.join(out, "gate3a_blocks.jsonl"), "w") as fh:
        for i, r in enumerate(pool.map(run_record, [(bench, s) for s in studies])):
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            for k in branches:
                branches[k] += r["branches"][k]
            unusable += [dict(record=r["key"], **u) for u in r["unusable"]]
            for sb in r["blocks"]:
                g1agree[0] += sb.pop("_g1_full_design_agreement")
                g1agree[1] += sb.pop("_g1_folds")
                accumulate(records, r["key"], sb)
            print(f"[{i + 1}/{len(studies)}] {r['key']} blocks {len(r['blocks'])} unusable {len(r['unusable'])}",
                  flush=True)
    summary = {"script": manifest["script"], "script_sha256": manifest["script_sha256"],
               "protocol_sha256": proto, "records_with_usable_blocks": len(records),
               "unusable_blocks": len(unusable), "log_rule_branches": branches,
               "g1_fold_vs_full_design_agreement": g1agree, "families": aggregate(records),
               "MDE": MDE, "bonferroni_alpha": ALPHA_BONF, "evaluation_rows_read": False}
    summary["qualifying_families"] = [f for f in TESTED if summary["families"][f]["qualifies"]]
    with open(os.path.join(out, "gate3a_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    with open(os.path.join(out, "gate3a_unusable.json"), "w") as fh:
        json.dump(unusable, fh, indent=1)
    print(json.dumps({k: v for k, v in summary.items() if k != "families"}, indent=1))
    for f, row in summary["families"].items():
        print(f, json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in row.items()}))
    print("gate3a_summary.json", sha_file(os.path.join(out, "gate3a_summary.json")))


def secondary_record(task):
    """Selections and fits on all design rows; scored once on the evaluation half (protocol section 4)."""
    bench, study = task
    audit = _bind(bench)
    key, m, _df, dev, test = audit.load_record(study)
    del _df
    out = {"key": key, "blocks": []}
    for b in m["blocks"]:
        O, X, Y = block_arrays(b, dev)
        Ot, Xt, Yt = block_arrays(b, test)
        if not all(np.isfinite(a).all() for a in (O, X, Y)):
            continue
        ok = np.isfinite(Ot).all(axis=1) & np.isfinite(Xt).all(axis=1) & np.isfinite(Yt).all(axis=1)
        Ot, Xt, Yt = Ot[ok], Xt[ok], Yt[ok]
        if len(Yt) == 0:
            continue
        jg = g1_pick(O, Y, X)
        row = {}
        for f in FAMILIES:
            best, loss = inner_select(f, O, X, Y, SEED)
            V = loss[0] - loss[1:]
            jh = int(np.flatnonzero(V == V.max())[0])
            sums = {}
            for arm, j in (("hat", jh), ("g1", jg)):
                Ftr, Fte = features(f, O, Ot, X[:, j], Xt[:, j])
                sums[arm] = float(bounded_sum(fit_predict(f, Ftr, Y, Fte)[[best[j + 1]]], Yt)[0])
            row[f] = {"hat": sums["hat"], "g1": sums["g1"], "n": float(Yt.size)}
        out["blocks"].append(row)
    return out


def cmd_secondary(a):
    from concurrent.futures import ProcessPoolExecutor
    out = os.path.abspath(a.out)
    if not os.path.exists(os.path.join(out, "gate3a_summary.json")):
        raise SystemExit("refusing: the secondary runs only after gate3a_summary.json exists")
    if os.path.exists(os.path.join(out, "gate3a_secondary.json")):
        raise SystemExit("gate3a_secondary.json exists: the secondary has already been run.")
    bench = os.path.abspath(a.bench)
    audit = _bind(bench)
    studies = json.loads((audit.ROOT / "paper4_design_models/p4d_design_model_freeze_manifest.json").read_text())["studies"]
    per = {f: [] for f in FAMILIES}
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for i, r in enumerate(pool.map(secondary_record, [(bench, s) for s in studies])):
            if not r["blocks"]:
                continue
            for f in FAMILIES:
                h = sum(b[f]["hat"] for b in r["blocks"])
                g = sum(b[f]["g1"] for b in r["blocks"])
                per[f].append(1 - h / g if g > 0 else 0.0)
            print(f"[{i + 1}/{len(studies)}] {r['key']}", flush=True)
    res = {"note": "secondary only; evaluation half previously opened by Paper III, Paper IV and the 9 Oct benchmark",
           "families": {f: {"delta_acq": float(np.mean(v)), "lower_98.33": boot_lower(v, ALPHA_BONF),
                            "records": len(v)} for f, v in per.items() if v}}
    with open(os.path.join(out, "gate3a_secondary.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    st = sub.add_parser("selftest")
    st.add_argument("--planted-blocks", type=int, default=4)
    st.add_argument("--null-reps", type=int, default=20)
    st.add_argument("--null-records", type=int, default=8)
    for name in ("run", "secondary"):
        rn = sub.add_parser(name)
        rn.add_argument("--bench", required=True)
        rn.add_argument("--out", default="gate3a_out")
        rn.add_argument("--workers", type=int, default=4)
        if name == "run":
            rn.add_argument("--protocol", required=True)
            rn.add_argument("--protocol-sha256", required=True)
    a = ap.parse_args()
    if a.cmd == "selftest":
        sys.exit(selftest(a))
    (cmd_run if a.cmd == "run" else cmd_secondary)(a)


if __name__ == "__main__":
    main()
