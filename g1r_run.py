#!/usr/bin/env python3
"""G1-R executable implementation v0.3.0 (protocol G1_REGULARISATION_PROTOCOL_v0.5.md). DRAFT until `freeze`.

Gaussian-only study of covariance regularisation and finite-sample acquisition decisions. Comparator A0 is the frozen
p5_gate3a.g1_pick (sha256 4600e2f3...), hash-verified and called unchanged. Every other arm is defined here in closed
form on the standardised sample covariance. Primary outcome: EXACT population regret of each pick,
(G_{j*} - G_jhat) / R_empty. No evaluation sample, no fitted predictor, no cross-fitting, no outer folds.

v0.2.0 (owner's code review, 7 Nov 2026): A1 gate on scaled C, v and scores with defined norms and the scaled clip
threshold; replicate-level gates are ASSERTIONS that stop the run; freeze binds the passing self-test report to the
current file hashes, records its checksum, and refuses to overwrite; complete shrinkage (delta == 1) distinguished
from exact score ties; secondary contrasts A2-S vs A2 and A2' vs A2 with sparse handling; GCV objective checked
direct-vs-SVD across the whole search set with the selected penalty verified; sampler moment check on six fixed
conditions; replicate-level integrity checker (g1r_integrity.py) required before analyze.

v0.3.0 (owner's second review): explicit finiteness and domain gates at runtime (NaN comparisons are false, so every
quantity is checked with isfinite before any tolerance test); g1r_integrity.py is a frozen dependency and its report is
bound to the replicate-log, freeze and checker checksums, all re-verified by analyze; analyze defines its result when no
contrast is assessable.

Commands: selftest | freeze | run [--cond a:b] | analyze.  Salt G1R-V02| fixed (owner decision).
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, sys, time
from pathlib import Path
import numpy as np

VERSION = "0.3.0"
PROTOCOL = "G1_REGULARISATION_PROTOCOL_v0.5.md"
SALT = "G1R-V02|"
HERE = Path(__file__).resolve().parent
FROZEN_DIR = HERE / "nhanes_fa" / "frozen"
P5_SHA = "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"
OUT = HERE / "g1r_out"; FREEZE = OUT / "FREEZE.json"; RUNLOG = OUT / "replicates.jsonl"; SELFTEST = OUT / "selftest.json"
INTEGRITY = OUT / "integrity.json"; CHECKER = HERE / "g1r_integrity.py"

D_O, D_T, P = 20, 20, 20; D = D_O + D_T + P
EPS, CLIP, PSI = 1e-6, 1e-12, 0.25
GATE_TOL = 1e-12                                 # A1 == A0 relative tolerance on C, v, scores
N_TRAIN = (102, 410, 820, 1640, 3280)
PHIS = (0.0, 0.5, 0.8, 0.95)
GAMMAS = (0.02, 0.05, 0.20)
A_GRID = (0.05, 0.15, 0.20)
HETERO = (0, 1)
R_REPS = 1000
LAMBDA_GRID = np.concatenate([[0.0], np.logspace(-6, 2, 81)])
LAMBDA_SET = frozenset(float(x) for x in LAMBDA_GRID)   # exact membership test for a selected penalty
BETAS = (0.0, 0.5, 1.0)
B_BOOT = 2000
SPARSE_MIN = 10
MARGIN_H1, MARGIN_H2 = 0.9, 1.1
S2_EQUAL = 0.5
S2_HETERO = 0.25 * 3.0 ** (np.arange(P) / (P - 1))        # 0.2500 ... 0.7500, mean 0.457484
BASE_SHARE_FACTOR = 2.0
MOMENT_N, MOMENT_BATCHES, MOMENT_Z, MOMENT_TOL = 400_000, 20, 4.0, 0.02

sys.path.insert(0, str(HERE))
import nova_d_calib as CAL   # noqa: E402

io, it, ip = slice(0, D_O), slice(D_O, D_O + D_T), slice(D_O + D_T, D)
ARMS = ("A0", "A2", "A2S", "A2p", "A3", "A4", "A2_b0.5", "A2_b1")


class GateError(RuntimeError):
    """A replicate-level gate failed. The run stops; nothing is recorded for the failing replicate."""


def sha_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def sha_text(s: str) -> str:
    return hashlib.sha256(s.encode()).hexdigest()


# ---------------------------------------------------------------- frozen comparator
_G = None


def frozen_gate3a():
    global _G
    if _G is None:
        p5 = FROZEN_DIR / "p5_gate3a.py"
        if not p5.exists() or sha_file(p5) != P5_SHA:
            raise SystemExit(f"REFUSED: {p5} missing or sha256 != {P5_SHA}")
        spec = importlib.util.spec_from_file_location("p5_gate3a_frozen", p5)
        m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); _G = m
    return _G


def alpha_frozen(n_train: int) -> float:
    return min(1.0, max(0.05, D / (D + n_train)))


def alpha_nofloor(n_train: int) -> float:
    return min(1.0, D / (D + n_train))


# ---------------------------------------------------------------- conditions and seeds (explicit)
def cond_string(c) -> str:
    return f"n_train={c['n_train']}|phi={c['phi']}|gamma={c['gamma']}|A={c['A']}|hetero={c['hetero']}"


def seed_for(*parts) -> int:
    return int(hashlib.sha256((SALT + "|".join(str(x) for x in parts)).encode()).hexdigest()[:15], 16)


def conditions():
    """Enumeration order outer->inner: n_train, phi, gamma, A, hetero. idx 0..359."""
    out, i = [], 0
    for n in N_TRAIN:
        for phi in PHIS:
            for g in GAMMAS:
                for A in A_GRID:
                    for h in HETERO:
                        out.append(dict(idx=i, n_train=n, phi=phi, gamma=g, A=A, hetero=h)); i += 1
    return out


def find_cond(conds, **spec):
    return next(c for c in conds if all(c[k] == v for k, v in spec.items()))


# ---------------------------------------------------------------- structure, population covariance, sampler
def structure(c):
    cell = CAL.build("M1", 1, c["gamma"], c["A"], 0.0, PSI)
    if not cell.feasible:
        raise SystemExit(f"REFUSED: infeasible {cond_string(c)}")
    rng = np.random.default_rng(seed_for("structure", cond_string(c)))
    r, s_eps2 = cell.r, cell.sigma2
    G = np.array([CAL.G_at(k, r, c["gamma"]) for k in range(1, P + 1)])
    Sig = c["phi"] ** np.abs(np.subtract.outer(np.arange(D_O), np.arange(D_O)))
    L = np.linalg.cholesky(Sig) if c["phi"] > 0 else np.eye(D_O)
    s2 = np.full(P, S2_EQUAL) if not c["hetero"] else S2_HETERO[rng.permutation(P)]
    def scaled_dirs(k, var):
        w = rng.normal(size=(k, D_O)); q = np.einsum("ij,jk,ik->i", w, Sig, w); return w * np.sqrt(var / q)[:, None]
    a = scaled_dirs(P, 1.0 - s2)
    w = rng.normal(size=(P, D_T)); w /= np.linalg.norm(w, axis=1, keepdims=True)
    lam = (w * np.sqrt(G / s2)[:, None]).T
    condvar = (lam ** 2 * s2[None, :]).sum(axis=1) + s_eps2
    cc = BASE_SHARE_FACTOR * condvar.max()
    b = np.vstack([scaled_dirs(1, cc - condvar[k])[0] for k in range(D_T)])
    R0 = G.sum() + D_T * s_eps2
    return dict(cond=c, r=r, s_eps2=s_eps2, Sig=Sig, L=L, s2=s2, a=a, lam=lam, b=b, G=G, jstar=int(np.argmax(G)),
                R0=R0, A_or=G.max() / R0, c_marg=cc, condvar=condvar)


def population_cov(st):
    Sig, a, b, lam, s2, se2 = st["Sig"], st["a"], st["b"], st["lam"], st["s2"], st["s_eps2"]
    Sop = Sig @ a.T; Sot = Sig @ b.T
    Spp = a @ Sig @ a.T + np.diag(s2); Stp = b @ Sig @ a.T + lam * s2[None, :]
    Stt = b @ Sig @ b.T + (lam * s2[None, :]) @ lam.T + se2 * np.eye(D_T)
    return np.block([[Sig, Sot, Sop], [Sot.T, Stt, Stp], [Sop.T, Stp.T, Spp]])


def sample(st, n, seed):
    rng = np.random.default_rng(seed)
    O = rng.normal(size=(n, D_O)) @ st["L"].T
    u = rng.normal(size=(n, P)) * np.sqrt(st["s2"])[None, :]
    X = O @ st["a"].T + u
    T = O @ st["b"].T + u @ st["lam"].T + rng.normal(scale=math.sqrt(st["s_eps2"]), size=(n, D_T))
    return O, T, X, u


def regret(st, j):
    return float((st["G"][st["jstar"]] - st["G"][j]) / st["R0"])


# ---------------------------------------------------------------- scoring forms
def standardise_cov(Sigma):
    sd = np.sqrt(np.diag(Sigma)); return Sigma / np.outer(sd, sd)


def schur(S, lam):
    sol = np.linalg.solve(S[io, io] + lam * np.eye(D_O), S[io, ip])
    return S[it, ip] - S[it, io] @ sol, np.diag(S[ip, ip] - S[ip, io] @ sol)


def residual(S, lam):
    Sinv = np.linalg.inv(S[io, io] + lam * np.eye(D_O)); BT, BP = Sinv @ S[io, it], Sinv @ S[io, ip]
    C = S[it, ip] - S[it, io] @ BP - BT.T @ S[io, ip] + BT.T @ S[io, io] @ BP
    v = np.diag(S[ip, ip] - S[ip, io] @ BP - BP.T @ S[io, ip] + BP.T @ S[io, io] @ BP)
    return C, v


def residual_via_identity(S, lam):
    C, v = schur(S, lam); Sinv = np.linalg.inv(S[io, io] + lam * np.eye(D_O)); BT, BP = Sinv @ S[io, it], Sinv @ S[io, ip]
    return C - lam * BT.T @ BP, v - lam * (BP ** 2).sum(axis=0)


def score(C, v, beta=0.0, floor=0.0, scale=1.0):
    vb = (1 - beta) * v + beta * v.mean() + floor
    return scale * (C * C).sum(axis=0) / np.maximum(vb, CLIP), vb


def frozen_algebra(S_std, alpha):
    """Frozen g1_pick algebra on a standardised covariance; returns (C, v, scores) in the shrunk scale."""
    Sa = (1 - alpha) * S_std + alpha * (np.trace(S_std) / D) * np.eye(D)
    C, v = schur(Sa, EPS); s, _ = score(C, v); return C, v, s


def a1_scaled(S_std, alpha):
    """Equivalence control: ridge lam, floor, scaled by (1 - alpha). Returns (C, v, scores) comparable to frozen_algebra."""
    lam, floor = (alpha + EPS) / (1 - alpha), alpha / (1 - alpha)
    C, v = schur(S_std, lam); s, vb = score(C, v, floor=floor, scale=(1 - alpha))
    return (1 - alpha) * C, (1 - alpha) * vb, s


def relmax(x, y):
    return float(np.max(np.abs(x - y)) / max(np.max(np.abs(y)), 1e-300))


def a1_gate(S_std, alpha):
    """Returns dict(C, v, s) relative max-abs discrepancies, argmax equality, and clip status at the scaled threshold."""
    Cm, vm, sm = frozen_algebra(S_std, alpha); C1, v1, s1 = a1_scaled(S_std, alpha)
    return dict(C=relmax(C1, Cm), v=relmax(v1, vm), s=relmax(s1, sm), argmax_equal=bool(np.argmax(s1) == np.argmax(sm)),
                clip_active=bool(min(vm.min(), v1.min()) <= CLIP), clip_threshold_unscaled=CLIP / (1 - alpha))


def pick(scores):
    return int(np.argmax(scores))


def top_tie_count(scores):
    return int(np.sum(scores == scores.max()))


# ---------------------------------------------------------------- GCV and Ledoit-Wolf
def gcv_objective_svd(O, Y, n):
    W, sv, _ = np.linalg.svd(O, full_matrices=False); e = sv ** 2
    Yt = W.T @ Y; resid0 = np.sum(Y ** 2) - np.sum(Yt ** 2); rowsq = np.sum(Yt ** 2, axis=1)
    vals = []
    for lam in LAMBDA_GRID:
        g = (n - 1) * lam; shrink = g / (e + g)
        rss = resid0 + np.sum(shrink ** 2 * rowsq); trH = np.sum(e / (e + g))
        vals.append((rss / n) / (1.0 - trH / n) ** 2)
    return np.array(vals)


def gcv_objective_direct(O, Y, n):
    vals = []
    for lam in LAMBDA_GRID:
        g = (n - 1) * lam; H = O @ np.linalg.solve(O.T @ O + g * np.eye(D_O), O.T)
        R = Y - H @ Y; vals.append((np.sum(R ** 2) / n) / (1.0 - np.trace(H) / n) ** 2)
    return np.array(vals)


def gcv_select(vals):
    best, bestv = None, np.inf
    for lam, val in zip(LAMBDA_GRID, vals):
        if val < bestv * (1 - 1e-12): best, bestv = lam, val
    return float(best), bool(best == LAMBDA_GRID[-1])


def ledoit_wolf_delta(Zc):
    n = Zc.shape[0]; Sn = Zc.T @ Zc / n; mu = np.trace(Sn) / D
    X2 = Zc ** 2; beta_ = np.sum(X2.T @ X2) / n - np.sum(Sn ** 2); delta_ = np.sum((Sn - mu * np.eye(D)) ** 2)
    beta_ = min(beta_ / n, delta_); return float(beta_ / delta_) if delta_ > 0 else 1.0


def ledoit_wolf_delta_direct(Zc):
    n = Zc.shape[0]; Sn = Zc.T @ Zc / n; mu = np.trace(Sn) / D
    beta_ = sum(np.sum((np.outer(z, z) - Sn) ** 2) for z in Zc) / n ** 2; gam = np.sum((Sn - mu * np.eye(D)) ** 2)
    return float(min(1.0, beta_ / gam)) if gam > 0 else 1.0


# ---------------------------------------------------------------- all arms on one sample (gates are assertions)
def finite(*arrs):
    return all(np.all(np.isfinite(np.asarray(a, dtype=float))) for a in arrs)


def scored(C, v, **kw):
    """score() with a finiteness gate on inputs and outputs (NaN > tol is False, so finiteness is checked explicitly)."""
    if not finite(C, v): raise GateError("non-finite conditional quantities")
    s, vb = score(C, v, **kw)
    if not finite(s, vb): raise GateError("non-finite scores")
    return s


def run_arms(O, T, X, n_train):
    g = frozen_gate3a()
    if not finite(O, T, X): raise GateError("non-finite sample")
    Z = np.column_stack([O, T, X]); Zc = Z - Z.mean(axis=0); sd = Zc.std(axis=0, ddof=1)
    if not (finite(sd) and np.all(sd > 0)): raise GateError("zero or non-finite column SD")
    Zs = Zc / sd; S = np.cov(Zs, rowvar=False, ddof=1)
    if not finite(S): raise GateError("non-finite covariance")
    a_n, a3 = alpha_frozen(n_train), alpha_nofloor(n_train)
    if not (a_n < 1.0 and a3 < 1.0): raise GateError("alpha >= 1")
    out, info = {}, {}
    out["A0"] = int(g.g1_pick(O, T, X))
    Cm, vm, sm = frozen_algebra(S, a_n)
    if not finite(Cm, vm, sm): raise GateError("non-finite frozen-mirror quantities")
    if pick(sm) != out["A0"]: raise GateError(f"mirror argmax {pick(sm)} != frozen pick {out['A0']}")
    gt = a1_gate(S, a_n)
    for k in ("C", "v", "s"):
        if not (math.isfinite(gt[k]) and 0.0 <= gt[k] <= GATE_TOL): raise GateError(f"A1 gate failed on {k}: {gt[k]!r}")
    if not gt["argmax_equal"]: raise GateError("A1 argmax != frozen argmax")
    if gt["clip_active"]: raise GateError("denominator clip active")
    info.update(A1_C=gt["C"], A1_v=gt["v"], A1_s=gt["s"])
    Os, Ys = Zs[:, io], Zs[:, D_O:]
    vals = gcv_objective_svd(Os, Ys, n_train)
    if not finite(vals): raise GateError("non-finite GCV objective")
    lam_c, ep = gcv_select(vals)
    if lam_c not in LAMBDA_SET: raise GateError(f"lambda_c {lam_c!r} not in the search set")
    CR, vR = residual(S, lam_c); CS, vS = schur(S, lam_c)
    out["A2"] = pick(scored(CR, vR)); out["A2S"] = pick(scored(CS, vS))
    out["A2_b0.5"] = pick(scored(CR, vR, beta=0.5)); out["A2_b1"] = pick(scored(CR, vR, beta=1.0))
    C0, v0 = schur(S, 0.0); out["A2p"] = pick(scored(C0, v0))
    C3, v3, s3 = frozen_algebra(S, a3)
    if not finite(C3, v3, s3): raise GateError("non-finite A3 quantities")
    out["A3"] = pick(s3)
    delta = ledoit_wolf_delta(Zs)
    if not (math.isfinite(delta) and 0.0 <= delta <= 1.0): raise GateError(f"Ledoit-Wolf intensity out of [0,1]: {delta!r}")
    S_lw = (1 - delta) * S + delta * (np.trace(S) / D) * np.eye(D)
    C4, v4 = schur(S_lw, EPS); s4 = scored(C4, v4); out["A4"] = pick(s4)
    kappa_hat = float(np.linalg.cond(S[io, io]))
    if not (math.isfinite(kappa_hat) and kappa_hat >= 1.0): raise GateError(f"invalid condition number {kappa_hat!r}")
    info.update(lambda_c=float(lam_c), endpoint_hit=bool(ep), delta_lw=float(delta), complete_shrinkage=bool(delta >= 1.0),
                A4_exact_top_tie=bool(top_tie_count(s4) > 1), kappa_hat=kappa_hat)
    for k, v in out.items():
        if not (isinstance(v, int) and 0 <= v < P): raise GateError(f"invalid pick {k}={v!r}")
    return out, info


def replicate(c, st, rep):
    O, T, X, _ = sample(st, c["n_train"], seed_for("rep", cond_string(c), rep))
    picks, info = run_arms(O, T, X, c["n_train"])
    return dict(idx=c["idx"], rep=rep, picks=picks, regret={k: regret(st, j) for k, j in picks.items()}, **info)


# ---------------------------------------------------------------- Objective A, designation, moment check
def pop_argmax_tables(c, st):
    S = standardise_cov(population_cov(st)); n = c["n_train"]
    res = dict(gaussian=int(st["jstar"]), A0=pick(frozen_algebra(S, alpha_frozen(n))[2]), A3=pick(frozen_algebra(S, alpha_nofloor(n))[2]),
               A2p=pick(score(*schur(S, 0.0))[0]))
    lam1 = (alpha_frozen(n) + EPS) / (1 - alpha_frozen(n))
    res["A2_matched_lambda_residual"] = pick(score(*residual(S, lam1))[0]); res["A2S_matched_lambda_schur"] = pick(score(*schur(S, lam1))[0])
    SOO, SOY = S[io, io], S[io, D_O:]; B0 = np.linalg.solve(SOO, SOY)
    pe = [float(np.trace((np.linalg.solve(SOO + l * np.eye(D_O), SOY) - B0).T @ SOO @ (np.linalg.solve(SOO + l * np.eye(D_O), SOY) - B0))) for l in LAMBDA_GRID]
    res["gcv_infinite_data_benchmark_lambda"] = float(LAMBDA_GRID[int(np.argmin(pe))])
    mu = np.trace(S) / D; pi_pop = np.sum(np.outer(np.diag(S), np.diag(S)) + S ** 2); gam_pop = np.sum((S - mu * np.eye(D)) ** 2)
    d_pop = min(1.0, pi_pop / (n * gam_pop)); res["lw_proxy_delta"] = float(d_pop)
    res["A4_at_lw_proxy"] = pick(score(*schur((1 - d_pop) * S + d_pop * mu * np.eye(D), EPS))[0])
    res["A0_flips"] = bool(res["A0"] != res["gaussian"]); return res


D1_LIST = [dict(gamma=0.02, A=A, hetero=0, n_train=410, phi=phi) for A in (0.15, 0.20, 0.05) for phi in (0.5, 0.8, 0.0, 0.95)]
D2_LIST = [dict(gamma=0.20, A=A, hetero=0, n_train=410, phi=phi) for A in (0.15, 0.20, 0.05) for phi in (0.5, 0.8, 0.0, 0.95)]


def designate(conds):
    d1 = d2 = None; log = {"D1": [], "D2": []}
    for pos, spec in enumerate(D1_LIST):
        c = find_cond(conds, **spec); t = pop_argmax_tables(c, structure(c)); log["D1"].append((c["idx"], t["A0_flips"]))
        if t["A0_flips"] and d1 is None: d1 = (c["idx"], pos)
    for pos, spec in enumerate(D2_LIST):
        c = find_cond(conds, **spec); t = pop_argmax_tables(c, structure(c)); log["D2"].append((c["idx"], not t["A0_flips"]))
        if not t["A0_flips"] and d2 is None: d2 = (c["idx"], pos)
    return d1, d2, log


def moment_conditions(conds, d1, d2):
    """Six fixed conditions: D1, D2, and four covering both variance constructions and both conditioning extremes."""
    return [d1[0], d2[0],
            find_cond(conds, n_train=102, phi=0.0, gamma=0.02, A=0.05, hetero=1)["idx"],
            find_cond(conds, n_train=3280, phi=0.95, gamma=0.20, A=0.20, hetero=1)["idx"],
            find_cond(conds, n_train=820, phi=0.95, gamma=0.05, A=0.05, hetero=0)["idx"],
            find_cond(conds, n_train=1640, phi=0.0, gamma=0.05, A=0.15, hetero=1)["idx"]]


def moment_check(st, seed):
    """Sampler instantiates the calibrated distribution: sampled G_j, R_empty, A vs analytic; batch-means z and absolute caps."""
    n = MOMENT_N; O, T, X, u = sample(st, n, seed); Rt = T - O @ st["b"].T
    Ghat = np.array([((Rt * u[:, [j]]).mean(axis=0) ** 2).sum() / u[:, j].var(ddof=1) for j in range(P)])
    R0row = (Rt ** 2).sum(axis=1); bsz = n // MOMENT_BATCHES
    R0b = np.array([R0row[b * bsz:(b + 1) * bsz].mean() for b in range(MOMENT_BATCHES)])
    Ab = np.array([(((Rt[b * bsz:(b + 1) * bsz] * u[b * bsz:(b + 1) * bsz, [st["jstar"]]]).mean(axis=0) ** 2).sum()
                    / u[b * bsz:(b + 1) * bsz, st["jstar"]].var(ddof=1)) / R0row[b * bsz:(b + 1) * bsz].mean() for b in range(MOMENT_BATCHES)])
    z = lambda vals, target: float(abs(vals.mean() - target) / (vals.std(ddof=1) / math.sqrt(len(vals))))
    rep = dict(G_abs_err_max_over_G1=float(np.max(np.abs(Ghat - st["G"]))), R0_rel_err=float(abs(R0b.mean() / st["R0"] - 1)), R0_z=z(R0b, st["R0"]),
               A_est=float(Ab.mean()), A_target=float(st["cond"]["A"]), A_abs_err=float(abs(Ab.mean() - st["cond"]["A"])), A_z=z(Ab, st["cond"]["A"]),
               marginal_var_rel_spread=float(np.max(np.abs(T.var(axis=0, ddof=1) / st["c_marg"] - 1))), argmaxG_sampled_is_jstar=bool(np.argmax(Ghat) == st["jstar"]))
    rep["ok"] = (rep["G_abs_err_max_over_G1"] < MOMENT_TOL and rep["R0_rel_err"] < MOMENT_TOL and rep["R0_z"] < MOMENT_Z
                 and rep["A_abs_err"] < MOMENT_TOL and rep["A_z"] < MOMENT_Z and rep["marginal_var_rel_spread"] < MOMENT_TOL)
    return rep


# ---------------------------------------------------------------- self-test
def cmd_selftest(a):
    t0 = time.time(); g = frozen_gate3a(); conds = conditions(); rng = np.random.default_rng(11)
    rep = {"version": VERSION, "hashes": current_hashes(), "n_conditions": len(conds)}   # same set the freeze binds to, incl. the checker
    # 1. mirror == frozen pick
    mism = 0
    for _ in range(200):
        n = int(rng.integers(60, 600)); O, T, X = rng.normal(size=(n, D_O)), rng.normal(size=(n, D_T)), rng.normal(size=(n, P)); T += 0.5 * O + 0.3 * X * rng.normal(size=P)
        Z = np.column_stack([O, T, X]); S = np.cov((Z - Z.mean(0)) / Z.std(0, ddof=1), rowvar=False, ddof=1)
        mism += int(pick(frozen_algebra(S, alpha_frozen(n))[2]) != g.g1_pick(O, T, X))
    rep["mirror_mismatches_of_200"] = mism
    # 2. A1 gate on C, v, scores
    alphas = sorted({alpha_frozen(n) for n in N_TRAIN} | {D / (D + m) for m in (102, 103, 409, 410, 819, 820)})
    mats = [standardise_cov(population_cov(structure(c))) for c in conds]
    for _ in range(1000):
        Araw = rng.normal(size=(D, D + 5)); mats.append(standardise_cov(Araw @ Araw.T / (D + 5) + 1e-3 * np.eye(D)))
    worst = {"C": 0.0, "v": 0.0, "s": 0.0}; amism = clipact = 0
    for S in mats:
        for al in alphas:
            assert al < 1.0; gt = a1_gate(S, al)
            for k in worst: worst[k] = max(worst[k], gt[k])
            amism += int(not gt["argmax_equal"]); clipact += int(gt["clip_active"])
    rep["A1_gate"] = dict(matrices=len(mats), alphas=len(alphas), max_rel_discrepancy=worst, argmax_mismatches=amism, clip_active_cases=clipact,
                          clip_threshold_unscaled_range=[CLIP / (1 - max(alphas)), CLIP / (1 - min(alphas))], norm="max-abs relative to max-abs of the frozen quantity")
    # 3. designation and Objective A
    d1, d2, dlog = designate(conds); rep["designation"] = dict(D1=d1, D2=d2, log=dlog)
    rep["objective_A"] = {c["idx"]: pop_argmax_tables(c, structure(c)) for c in conds}
    # 4. calibration
    rep["calibration_max_abs_A_err"] = float(max(abs(structure(c)["A_or"] - c["A"]) for c in conds))
    # 5. sampler moment check on six fixed conditions
    mc = moment_conditions(conds, d1, d2) if (d1 and d2) else []
    rep["moment_checks"] = {str(i): moment_check(structure(conds[i]), seed_for("moment", i)) for i in mc}
    # 6. exact-regret identity via covariance route
    worst_r = 0.0
    for key in [k for k in (d1, d2) if k]:
        c = conds[key[0]]; st = structure(c); Sig = population_cov(st)
        for j in rng.integers(0, P, 10):
            idx = list(range(D_O)) + [D_O + D_T + int(j)]; A_ = Sig[np.ix_(idx, idx)]; Bm = Sig[np.ix_(idx, range(D_O, D_O + D_T))]
            Rj = np.trace(Sig[it, it] - Bm.T @ np.linalg.solve(A_, Bm)); R0 = np.trace(Sig[it, it] - Sig[it, io] @ np.linalg.solve(Sig[io, io], Sig[io, it]))
            worst_r = max(worst_r, abs((R0 - Rj) / R0 - st["G"][j] / st["R0"]))
    rep["exact_regret_identity_max_abs_err"] = float(worst_r)
    # 7. GCV: objective direct vs SVD over the whole search set on 20 random problems; selected penalty identical; penalty-scaling identity
    gcv_worst, gcv_sel_mism = 0.0, 0
    for _ in range(20):
        n = int(rng.integers(102, 900)); O = rng.normal(size=(n, D_O)) @ np.linalg.cholesky(0.7 ** np.abs(np.subtract.outer(np.arange(D_O), np.arange(D_O))))
        Y = O @ rng.normal(size=(D_O, D_T + P)) * 0.3 + rng.normal(size=(n, D_T + P)); O -= O.mean(0); Y -= Y.mean(0)
        v1, v2 = gcv_objective_svd(O, Y, n), gcv_objective_direct(O, Y, n); gcv_worst = max(gcv_worst, float(np.max(np.abs(v1 - v2) / np.abs(v2))))
        gcv_sel_mism += int(gcv_select(v1)[0] != gcv_select(v2)[0])
    rep["gcv_objective_direct_vs_svd_max_rel_err"] = gcv_worst; rep["gcv_selected_penalty_mismatches_of_20"] = gcv_sel_mism
    O = rng.normal(size=(300, D_O)); Y = rng.normal(size=(300, D_T + P)); Z = np.column_stack([O, Y]); Zs = (Z - Z.mean(0)) / Z.std(0, ddof=1)
    S = np.cov(Zs, rowvar=False, ddof=1); lam = 0.3; gam = 299 * lam
    lhs = np.linalg.solve(S[io, io] + lam * np.eye(D_O), S[io, D_O:]); rhs = np.linalg.solve(Zs[:, io].T @ Zs[:, io] + gam * np.eye(D_O), Zs[:, io].T @ Zs[:, D_O:])
    rep["gcv_penalty_identity_max_abs_err"] = float(np.max(np.abs(lhs - rhs)))
    # 8. Ledoit-Wolf independent formula; residual identity; forms
    rep["lw_vs_direct_formula_abs_err"] = float(abs(ledoit_wolf_delta_direct(Zs) - ledoit_wolf_delta(Zs)))
    try:
        from sklearn.covariance import ledoit_wolf as sk_lw
        rep["lw_vs_sklearn_abs_err"] = float(abs(sk_lw(Zs, assume_centered=True)[1] - ledoit_wolf_delta(Zs)))
    except Exception as e:  # noqa: BLE001
        rep["lw_vs_sklearn_abs_err"] = f"sklearn unavailable here: {type(e).__name__} (direct-formula check is the gate)"
    CR, vR = residual(S, lam); CRi, vRi = residual_via_identity(S, lam); C0a, v0a = residual(S, 0.0); C0b, v0b = schur(S, 0.0)
    rep["residual_identity_max_abs_err"] = float(max(np.max(np.abs(CR - CRi)), np.max(np.abs(vR - vRi))))
    rep["forms_coincide_at_lambda0_max_abs_err"] = float(max(np.max(np.abs(C0a - C0b)), np.max(np.abs(v0a - v0b))))
    rep["forms_differ_at_lambda_pos"] = bool(np.max(np.abs(CR - schur(S, lam)[0])) > 1e-6)
    rep["cp_upper_k0_n1000"] = clopper_pearson_upper(0, 1000); rep["cp_upper_k0_n1000_expected"] = 1 - 0.05 ** (1 / 1000)
    # 9. replicate gates exercised: 5 replicates at the two designated conditions must pass as assertions
    for key in [k for k in (d1, d2) if k]:
        c = conds[key[0]]; st = structure(c); [replicate(c, st, r) for r in range(5)]
    rep["replicate_gate_smoke"] = "5 replicates at D1 and D2 passed all assertions" if (d1 and d2) else "skipped"
    c = conds[0]; st = structure(c); t1 = time.perf_counter(); [replicate(c, st, r) for r in range(5)]; rep["sec_per_replicate_n102"] = (time.perf_counter() - t1) / 5
    c = conds[-1]; st = structure(c); t1 = time.perf_counter(); [replicate(c, st, r) for r in range(5)]; rep["sec_per_replicate_n3280"] = (time.perf_counter() - t1) / 5
    rep["ok"] = (mism == 0 and all(w < GATE_TOL for w in worst.values()) and amism == 0 and clipact == 0 and d1 is not None and d2 is not None
                 and rep["calibration_max_abs_A_err"] < 1e-9 and len(rep["moment_checks"]) == 6 and all(m["ok"] for m in rep["moment_checks"].values())
                 and worst_r < 1e-12 and gcv_worst < 1e-10 and gcv_sel_mism == 0 and rep["gcv_penalty_identity_max_abs_err"] < 1e-12
                 and rep["lw_vs_direct_formula_abs_err"] < 1e-10 and (not isinstance(rep["lw_vs_sklearn_abs_err"], float) or rep["lw_vs_sklearn_abs_err"] < 1e-10)
                 and rep["residual_identity_max_abs_err"] < 1e-12 and rep["forms_coincide_at_lambda0_max_abs_err"] < 1e-12 and rep["forms_differ_at_lambda_pos"]
                 and abs(rep["cp_upper_k0_n1000"] - rep["cp_upper_k0_n1000_expected"]) < 1e-12)
    rep["elapsed_s"] = time.time() - t0
    OUT.mkdir(exist_ok=True); SELFTEST.write_text(json.dumps(rep, indent=1, default=str))
    print(json.dumps({k: v for k, v in rep.items() if k not in ("objective_A",)}, indent=1, default=str)); print("SELFTEST", "PASS" if rep["ok"] else "FAIL")
    return 0 if rep["ok"] else 1


# ---------------------------------------------------------------- freeze / run / analyze
def current_hashes():
    return {"p5_gate3a.py": sha_file(FROZEN_DIR / "p5_gate3a.py"), "g1r_run.py": sha_file(__file__), "nova_d_calib.py": sha_file(HERE / "nova_d_calib.py"),
            "g1r_integrity.py": sha_file(CHECKER)}


def cmd_freeze(a):
    OUT.mkdir(exist_ok=True)
    if FREEZE.exists(): raise SystemExit("REFUSED: FREEZE.json exists; a freeze is never overwritten")
    if RUNLOG.exists(): raise SystemExit("REFUSED: outputs exist")
    if not (HERE / PROTOCOL).exists(): raise SystemExit(f"REFUSED: {PROTOCOL} not found")
    if not SELFTEST.exists(): raise SystemExit("REFUSED: no self-test report")
    st = json.loads(SELFTEST.read_text())
    if not st.get("ok"): raise SystemExit("REFUSED: self-test report is not a pass")
    if st.get("version") != VERSION or st.get("hashes") != current_hashes(): raise SystemExit("REFUSED: self-test report is stale (version or file hashes differ from current files)")
    rec = dict(protocol=PROTOCOL, protocol_sha256=sha_file(HERE / PROTOCOL), **current_hashes(), selftest_sha256=sha_file(SELFTEST), version=VERSION,
               salt=SALT, n_train=N_TRAIN, phis=PHIS, gammas=GAMMAS, A_grid=A_GRID, hetero=HETERO, reps=R_REPS, lambda_grid_len=len(LAMBDA_GRID),
               betas=BETAS, margins=(MARGIN_H1, MARGIN_H2), sparse_min=SPARSE_MIN, B_boot=B_BOOT, gate_tol=GATE_TOL,
               D1=st["designation"]["D1"], D2=st["designation"]["D2"], moment_conditions=list(st["moment_checks"].keys()),
               n_conditions=len(conditions()), frozen_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    FREEZE.write_text(json.dumps(rec, indent=1)); print(json.dumps(rec, indent=1)); return 0


def check_frozen():
    if not FREEZE.exists(): raise SystemExit("REFUSED: not frozen")
    f = json.loads(FREEZE.read_text())
    if sha_file(HERE / PROTOCOL) != f["protocol_sha256"]: raise SystemExit("REFUSED: protocol changed since freeze")
    for k, v in current_hashes().items():
        if f[k] != v: raise SystemExit(f"REFUSED: {k} changed since freeze")
    if not SELFTEST.exists() or sha_file(SELFTEST) != f["selftest_sha256"]: raise SystemExit("REFUSED: self-test report changed since freeze")
    return f


def cmd_run(a):
    check_frozen(); conds = conditions(); lo, hi = (0, len(conds)) if a.cond is None else tuple(int(x) for x in a.cond.split(":"))
    done = set()
    if RUNLOG.exists():
        for line in RUNLOG.read_text().splitlines():
            if line.strip(): d = json.loads(line); done.add((d["idx"], d["rep"]))
    t0 = time.time(); k = 0
    with RUNLOG.open("a") as fh:
        for c in conds[lo:hi]:
            st = structure(c)
            for r in range(R_REPS):
                if (c["idx"], r) in done: continue
                fh.write(json.dumps(replicate(c, st, r)) + "\n"); k += 1          # GateError propagates and stops the run
            fh.flush(); print(f"[{time.strftime('%H:%M:%S')}] cond {c['idx']} done; total {k}; {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


def clopper_pearson_upper(k, n, conf=0.95):
    """One-sided 95 % Clopper-Pearson upper bound on P(jhat != j*): p with P(Bin(n, p) <= k) = 1 - conf, by bisection on the exact log-space binomial CDF."""
    if k >= n: return 1.0
    lc = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) for i in range(k + 1)]
    def cdf(p):
        if p <= 0: return 1.0
        if p >= 1: return 0.0
        return sum(math.exp(lc[i] + i * math.log(p) + (n - i) * math.log1p(-p)) for i in range(k + 1))
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if cdf(mid) > 1 - conf: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def cmd_analyze(a):
    f = check_frozen(); conds = conditions()
    if not INTEGRITY.exists(): raise SystemExit("REFUSED: no integrity report; run g1r_integrity.py first")
    ig = json.loads(INTEGRITY.read_text())
    if not ig.get("ok"): raise SystemExit("REFUSED: integrity report is not a PASS")
    if ig.get("replicates_sha256") != sha_file(RUNLOG): raise SystemExit("REFUSED: replicate log changed since the integrity report")
    if ig.get("freeze_sha256") != sha_file(FREEZE): raise SystemExit("REFUSED: freeze record changed since the integrity report")
    if ig.get("checker_sha256") != sha_file(CHECKER) or ig.get("checker_sha256") != f["g1r_integrity.py"]: raise SystemExit("REFUSED: integrity checker differs from the frozen one or from the report")
    recs = [json.loads(l) for l in RUNLOG.read_text().splitlines() if l.strip()]; by = {}
    for r in recs: by.setdefault(r["idx"], []).append(r)
    missing = [c["idx"] for c in conds if len(by.get(c["idx"], [])) != R_REPS]
    if missing: raise SystemExit(f"REFUSED: {len(missing)} conditions incomplete")
    rng = np.random.default_rng(seed_for("analyze")); out = {"freeze": f, "conditions": []}
    def contrast(ra, rb, factor):          # bootstrap of mean(factor * ra - rb); LB = 5th percentile (one-sided 95 %)
        q = factor * ra - rb; bs = np.array([q[rng.integers(0, len(q), len(q))].mean() for _ in range(B_BOOT)])
        return dict(mean=float(q.mean()), lcb95=float(np.percentile(bs, 5)), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
    def pair(R, sparse, x, y, with_margins):
        if sparse[x] or sparse[y]: return {"assessable": False, "reason": f"sparse regret in {x} or {y}; Clopper-Pearson bounds reported in arms"}
        d = dict(assessable=True, plain=contrast(R[x], R[y], 1.0))
        if with_margins: d["Q1"] = contrast(R[x], R[y], MARGIN_H1); d["Q2"] = contrast(R[x], R[y], MARGIN_H2)
        return d
    for c in conds:
        rs = by[c["idx"]]; R = {arm: np.array([r["regret"][arm] for r in rs]) for arm in ARMS}; st = structure(c)
        row = dict(c, jstar=st["jstar"], arms={}, sparse={}, contrasts_vs_A0={}, secondary={})
        for arm in ARMS:
            nz = int((R[arm] > 0).sum()); mis = int(sum(r["picks"][arm] != st["jstar"] for r in rs))
            row["arms"][arm] = dict(mean_regret=float(R[arm].mean()), p_optimal=1 - mis / len(rs), nonzero=nz, cp95_upper_misselection=clopper_pearson_upper(mis, len(rs)))
            row["sparse"][arm] = nz < SPARSE_MIN
        for arm in ARMS[1:]: row["contrasts_vs_A0"][arm] = pair(R, row["sparse"], "A0", arm, True)
        row["secondary"]["A2S_vs_A2_form_effect"] = pair(R, row["sparse"], "A2", "A2S", False)      # mean(r_A2 - r_A2S)
        row["secondary"]["A2p_vs_A2_tuning_effect"] = pair(R, row["sparse"], "A2", "A2p", False)
        row["A1_max_rel"] = {k: float(max(r[f"A1_{k}"] for r in rs)) for k in ("C", "v", "s")}
        row["endpoint_hit_rate"] = float(np.mean([r["endpoint_hit"] for r in rs])); row["complete_shrinkage_rate"] = float(np.mean([r["complete_shrinkage"] for r in rs]))
        row["A4_exact_top_tie_rate"] = float(np.mean([r["A4_exact_top_tie"] for r in rs]))
        row["lambda_c_median"] = float(np.median([r["lambda_c"] for r in rs])); row["delta_lw_median"] = float(np.median([r["delta_lw"] for r in rs]))
        row["kappa_hat_median"] = float(np.median([r["kappa_hat"] for r in rs]))
        out["conditions"].append(row)
    d1, d2 = f["D1"][0], f["D2"][0]; get = lambda i: next(r for r in out["conditions"] if r["idx"] == i)
    c1, c2 = get(d1)["contrasts_vs_A0"]["A2"], get(d2)["contrasts_vs_A0"]["A2"]
    assessable = [(r["contrasts_vs_A0"]["A2"]["plain"]["mean"], r["idx"]) for r in out["conditions"] if r["contrasts_vs_A0"]["A2"]["assessable"]]
    out["decisions"] = dict(H1_D1=("not assessable" if not c1["assessable"] else bool(c1["Q1"]["lcb95"] > 0)),
                            H2_D2=("not assessable" if not c2["assessable"] else bool(c2["Q2"]["lcb95"] >= 0)),
                            n_conditions_with_assessable_A2_contrast=len(assessable),
                            grid_max_plain_A2_gain_upward_biased=(max(assessable) if assessable else {"defined": False, "reason": "no condition has an assessable A2 contrast (all sparse)"}))
    (OUT / "g1r_result.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out["decisions"], indent=1)); return 0


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest); sub.add_parser("freeze").set_defaults(fn=cmd_freeze)
    r = sub.add_parser("run"); r.add_argument("--cond", default=None); r.set_defaults(fn=cmd_run)
    sub.add_parser("analyze").set_defaults(fn=cmd_analyze)
    a = ap.parse_args(); sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
