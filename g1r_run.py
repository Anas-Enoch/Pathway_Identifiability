#!/usr/bin/env python3
"""G1-R executable implementation (protocol G1_REGULARISATION_PROTOCOL_v0.3.md). DRAFT until `freeze`.

Gaussian-only study of covariance regularisation and finite-sample acquisition decisions. The comparator A0 is the
frozen p5_gate3a.g1_pick (sha256 4600e2f3...), hash-verified and called unchanged. Every other arm is defined
here in closed form on the standardised sample covariance. The primary outcome is the EXACT population regret of
each pick, (G_{j*} - G_jhat) / R_empty, from the condition's population covariance; no evaluation sample, no fitted
predictor, no cross-fitting, no outer folds.

Commands: selftest | freeze | run [--cond a:b] | analyze.  `run` refuses unless frozen and hashes match.
Salt G1R-V02| is fixed (owner decision 7 Nov 2026) and is not re-salted when the document version changes.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, sys, time
from pathlib import Path
import numpy as np

VERSION = "0.1.0"
PROTOCOL = "G1_REGULARISATION_PROTOCOL_v0.3.md"
SALT = "G1R-V02|"
HERE = Path(__file__).resolve().parent
FROZEN_DIR = HERE / "nhanes_fa" / "frozen"
P5_SHA = "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"
OUT = HERE / "g1r_out"; FREEZE = OUT / "FREEZE.json"; RUNLOG = OUT / "replicates.jsonl"

# ---- design constants (protocol v0.3 s4-s7) -----------------------------------------------------------------
D_O, D_T, P = 20, 20, 20; D = D_O + D_T + P
EPS, CLIP, PSI = 1e-6, 1e-12, 0.25
N_TRAIN = (102, 410, 820, 1640, 3280)
PHIS = (0.0, 0.5, 0.8, 0.95)
GAMMAS = (0.02, 0.05, 0.20)
A_GRID = (0.05, 0.15, 0.20)
HETERO = (0, 1)
R_REPS = 1000
LAMBDA_GRID = np.concatenate([[0.0], np.logspace(-6, 2, 81)])
BETAS = (0.0, 0.5, 1.0)
B_BOOT = 2000
SPARSE_MIN = 10
MARGIN_H1, MARGIN_H2 = 0.9, 1.1
S2_EQUAL = 0.5
S2_HETERO = 0.25 * 3.0 ** (np.arange(P) / (P - 1))        # 0.2500 ... 0.7500, mean 0.457484
BASE_SHARE_FACTOR = 2.0

sys.path.insert(0, str(HERE))
import nova_d_calib as CAL   # noqa: E402  (r, sigma_eps^2 for target A; H == 0)

io, it, ip = slice(0, D_O), slice(D_O, D_O + D_T), slice(D_O + D_T, D)


def sha_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


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
    """Enumeration order: n_train, phi, gamma, A, hetero (outer to inner). idx 0..359."""
    out, i = [], 0
    for n in N_TRAIN:
        for phi in PHIS:
            for g in GAMMAS:
                for A in A_GRID:
                    for h in HETERO:
                        out.append(dict(idx=i, n_train=n, phi=phi, gamma=g, A=A, hetero=h)); i += 1
    return out


# ---------------------------------------------------------------- structure, population covariance, sampler
def structure(c):
    cell = CAL.build("M1", 1, c["gamma"], c["A"], 0.0, PSI)       # H == 0; rank/mech irrelevant beyond gamma
    if not cell.feasible:
        raise SystemExit(f"REFUSED: infeasible {cond_string(c)}")
    rng = np.random.default_rng(seed_for("structure", cond_string(c)))
    r, s_eps2 = cell.r, cell.sigma2
    G = np.array([CAL.G_at(k, r, c["gamma"]) for k in range(1, P + 1)])     # candidate index == Gaussian rank
    Sig = c["phi"] ** np.abs(np.subtract.outer(np.arange(D_O), np.arange(D_O)))
    L = np.linalg.cholesky(Sig) if c["phi"] > 0 else np.eye(D_O)
    s2 = np.full(P, S2_EQUAL) if not c["hetero"] else S2_HETERO[rng.permutation(P)]
    def scaled_dirs(k, var):
        w = rng.normal(size=(k, D_O)); q = np.einsum("ij,jk,ik->i", w, Sig, w); return w * np.sqrt(var / q)[:, None]
    a = scaled_dirs(P, 1.0 - s2)                                              # Var(a_j'O) = 1 - s2_j, per candidate
    w = rng.normal(size=(P, D_T)); w /= np.linalg.norm(w, axis=1, keepdims=True)
    lam = (w * np.sqrt(G / s2)[:, None]).T                                  # (D_T, P): s2_j sum_k lam^2 = G_j
    condvar = (lam ** 2 * s2[None, :]).sum(axis=1) + s_eps2                   # Var(T_k | O)
    cc = BASE_SHARE_FACTOR * condvar.max()
    b = np.vstack([scaled_dirs(1, cc - condvar[k])[0] for k in range(D_T)])
    R0 = G.sum() + D_T * s_eps2
    return dict(cond=c, r=r, s_eps2=s_eps2, Sig=Sig, L=L, s2=s2, a=a, lam=lam, b=b, G=G, jstar=int(np.argmax(G)),
                R0=R0, A_or=G.max() / R0, c_marg=cc, condvar=condvar)


def population_cov(st):
    Sig, a, b, lam, s2, se2 = st["Sig"], st["a"], st["b"], st["lam"], st["s2"], st["s_eps2"]
    Soo = Sig; Sop = Sig @ a.T; Sot = Sig @ b.T
    Spp = a @ Sig @ a.T + np.diag(s2)
    Stp = b @ Sig @ a.T + lam * s2[None, :]
    Stt = b @ Sig @ b.T + (lam * s2[None, :]) @ lam.T + se2 * np.eye(D_T)
    return np.block([[Soo, Sot, Sop], [Sot.T, Stt, Stp], [Sop.T, Stp.T, Spp]])


def sample(st, n, seed):
    rng = np.random.default_rng(seed)
    O = rng.normal(size=(n, D_O)) @ st["L"].T
    u = rng.normal(size=(n, P)) * np.sqrt(st["s2"])[None, :]
    X = O @ st["a"].T + u
    T = O @ st["b"].T + u @ st["lam"].T + rng.normal(scale=math.sqrt(st["s_eps2"]), size=(n, D_T))
    return O, T, X


def regret(st, j):
    return float((st["G"][st["jstar"]] - st["G"][j]) / st["R0"])


# ---------------------------------------------------------------- scoring forms on a (standardised) covariance
def standardise_cov(Sigma):
    sd = np.sqrt(np.diag(Sigma)); return Sigma / np.outer(sd, sd)


def schur(S, lam):
    sol = np.linalg.solve(S[io, io] + lam * np.eye(D_O), S[io, ip])
    C = S[it, ip] - S[it, io] @ sol; v = np.diag(S[ip, ip] - S[ip, io] @ sol)
    return C, v


def residual(S, lam):
    Sinv = np.linalg.inv(S[io, io] + lam * np.eye(D_O))
    BT, BP = Sinv @ S[io, it], Sinv @ S[io, ip]
    C = S[it, ip] - S[it, io] @ BP - BT.T @ S[io, ip] + BT.T @ S[io, io] @ BP
    v = np.diag(S[ip, ip] - S[ip, io] @ BP - BP.T @ S[io, ip] + BP.T @ S[io, io] @ BP)
    return C, v


def residual_via_identity(S, lam):
    """C^R = C^S - lam B_T' B_P ;  v^R_j = v^S_j - lam ||B_P,j||^2  (owner's identity; used as a gate)."""
    C, v = schur(S, lam)
    Sinv = np.linalg.inv(S[io, io] + lam * np.eye(D_O)); BT, BP = Sinv @ S[io, it], Sinv @ S[io, ip]
    return C - lam * BT.T @ BP, v - lam * (BP ** 2).sum(axis=0)


def score(C, v, beta=0.0, floor=0.0, scale=1.0):
    vb = (1 - beta) * v + beta * v.mean() + floor
    return scale * (C * C).sum(axis=0) / np.maximum(vb, CLIP), bool(np.min(vb) <= CLIP)


def frozen_algebra(S_std, alpha):
    """The frozen g1_pick algebra on a standardised covariance: shrink, condition with ridge EPS, score."""
    Sa = (1 - alpha) * S_std + alpha * (np.trace(S_std) / D) * np.eye(D)
    C, v = schur(Sa, EPS)
    return score(C, v)[0]


def pick(scores):
    return int(np.argmax(scores))      # ties -> lowest index (NumPy)


# ---------------------------------------------------------------- GCV (s4.2a) and Ledoit-Wolf (s4.6)
def gcv_lambda(O, Y, n):
    """Shared-penalty multivariate ridge GCV over LAMBDA_GRID; gamma = (n - 1) lambda_c. Returns (lambda_c, endpoint_hit)."""
    W, sv, _ = np.linalg.svd(O, full_matrices=False); e = sv ** 2            # O = W diag(sv) V'
    Yt = W.T @ Y; resid0 = np.sum(Y ** 2) - np.sum(Yt ** 2)                   # ||Y - W W' Y||^2
    rowsq = np.sum(Yt ** 2, axis=1)
    best, bestv = None, np.inf
    for lam in LAMBDA_GRID:
        g = (n - 1) * lam
        shrink = g / (e + g)                                                   # 1 - e/(e+g)
        rss = resid0 + np.sum(shrink ** 2 * rowsq); trH = np.sum(e / (e + g))
        val = (rss / n) / (1.0 - trH / n) ** 2
        if val < bestv * (1 - 1e-12):                                           # strict improvement; ties -> smallest lam
            best, bestv = lam, val
    return float(best), bool(best == LAMBDA_GRID[-1])


def ledoit_wolf_delta(Zc):
    """Ledoit & Wolf (2004) intensity on centred data with the 1/n convention (scikit-learn ledoit_wolf, assume_centered=True)."""
    n = Zc.shape[0]; Sn = Zc.T @ Zc / n; mu = np.trace(Sn) / D
    X2 = Zc ** 2
    beta_ = np.sum(X2.T @ X2) / n - np.sum(Sn ** 2)                            # sum_ij (1/n) sum_t (z_ti z_tj)^2 - S_ij^2
    delta_ = np.sum((Sn - mu * np.eye(D)) ** 2)
    beta_ = min(beta_ / n, delta_)                                            # as in sklearn: shrinkage = min(beta/delta, 1)
    return float(beta_ / delta_) if delta_ > 0 else 1.0


# ---------------------------------------------------------------- all arms on one sample
ARMS = ("A0", "A2", "A2S", "A2p", "A3", "A4", "A2_b0.5", "A2_b1")


def run_arms(O, T, X, st, n_train):
    g = frozen_gate3a()
    Z = np.column_stack([O, T, X]); Zc = Z - Z.mean(axis=0); sd = Zc.std(axis=0, ddof=1); Zs = Zc / sd
    S = np.cov(Zs, rowvar=False, ddof=1)
    a_n = alpha_frozen(n_train); a3 = alpha_nofloor(n_train); assert a_n < 1.0 and a3 < 1.0
    out, flags = {}, {}
    # A0 frozen call; mirror scores for the A1 gate
    out["A0"] = int(g.g1_pick(O, T, X))
    s_mirror = frozen_algebra(S, a_n); flags["mirror_eq_frozen"] = bool(pick(s_mirror) == out["A0"])
    # A1 equivalence control (gate quantity)
    lam1, floor1 = (a_n + EPS) / (1 - a_n), a_n / (1 - a_n)
    C1, v1 = schur(S, lam1); s1, clip1 = score(C1, v1, floor=floor1, scale=(1 - a_n))
    flags["A1_rel_discrepancy"] = float(np.max(np.abs(s1 - s_mirror)) / max(np.max(np.abs(s_mirror)), 1e-300))
    flags["A1_clip_active"] = clip1
    # A2 decoupled (residual form, GCV), A2S (Schur form, same lambda), beta variants
    Os, Ys = Zs[:, io], Zs[:, D_O:]
    lam_c, ep = gcv_lambda(Os, Ys, n_train)
    CR, vR = residual(S, lam_c); CS, vS = schur(S, lam_c)
    out["A2"] = pick(score(CR, vR)[0]); out["A2S"] = pick(score(CS, vS)[0])
    out["A2_b0.5"] = pick(score(CR, vR, beta=0.5)[0]); out["A2_b1"] = pick(score(CR, vR, beta=1.0)[0])
    # A2' lambda_c = 0 (forms coincide)
    C0, v0 = schur(S, 0.0); out["A2p"] = pick(score(C0, v0)[0])
    # A3 vanishing floor (frozen algebra at alpha without floor)
    out["A3"] = pick(frozen_algebra(S, a3))
    # A4 Ledoit-Wolf: delta from 1/n formula on centred standardised data, applied to the ddof=1 covariance
    delta = ledoit_wolf_delta(Zs)
    S_lw = (1 - delta) * S + delta * (np.trace(S) / D) * np.eye(D)
    C4, v4 = schur(S_lw, EPS); s4, _ = score(C4, v4); out["A4"] = pick(s4)
    flags["A4_degenerate_tie"] = bool(np.allclose(s4, s4[0]))
    kappa_hat = float(np.linalg.cond(S[io, io]))
    return out, dict(lambda_c=lam_c, endpoint_hit=ep, delta_lw=delta, kappa_hat=kappa_hat, **flags)


def replicate(c, st, rep):
    O, T, X = sample(st, c["n_train"], seed_for("rep", cond_string(c), rep))
    picks, info = run_arms(O, T, X, st, c["n_train"])
    return dict(idx=c["idx"], rep=rep, picks=picks, regret={k: regret(st, j) for k, j in picks.items()}, **info)


# ---------------------------------------------------------------- Objective A and designation (exact)
def pop_argmax_tables(c, st):
    S = standardise_cov(population_cov(st)); n = c["n_train"]
    res = dict(gaussian=int(st["jstar"]))
    res["A0"] = pick(frozen_algebra(S, alpha_frozen(n))); res["A3"] = pick(frozen_algebra(S, alpha_nofloor(n)))
    res["A2p"] = pick(score(*schur(S, 0.0))[0])
    lam1 = (alpha_frozen(n) + EPS) / (1 - alpha_frozen(n))
    res["A2_matched_lambda_residual"] = pick(score(*residual(S, lam1))[0]); res["A2S_matched_lambda_schur"] = pick(score(*schur(S, lam1))[0])
    # infinite-data population benchmark for GCV: PE(lambda) in population-standardised coordinates (NOT finite-n expected GCV)
    SOO, SOY = S[io, io], S[io, D_O:]; B0 = np.linalg.solve(SOO, SOY)
    pe = [float(np.trace((np.linalg.solve(SOO + l * np.eye(D_O), SOY) - B0).T @ SOO @ (np.linalg.solve(SOO + l * np.eye(D_O), SOY) - B0))) for l in LAMBDA_GRID]
    res["gcv_infinite_data_benchmark_lambda"] = float(LAMBDA_GRID[int(np.argmin(pe))])
    # Ledoit-Wolf deterministic diagnostic proxy (Gaussian fourth moments; not the expectation of the estimated clipped intensity)
    mu = np.trace(S) / D; pi_pop = np.sum(np.outer(np.diag(S), np.diag(S)) + S ** 2); gam_pop = np.sum((S - mu * np.eye(D)) ** 2)
    d_pop = min(1.0, pi_pop / (n * gam_pop)); res["lw_proxy_delta"] = float(d_pop)
    res["A4_at_lw_proxy"] = pick(score(*schur((1 - d_pop) * S + d_pop * mu * np.eye(D), EPS))[0])
    res["A0_flips"] = bool(res["A0"] != res["gaussian"])
    return res


D1_LIST = [dict(gamma=0.02, A=A, hetero=0, n_train=410, phi=phi) for A in (0.15, 0.20, 0.05) for phi in (0.5, 0.8, 0.0, 0.95)]
D2_LIST = [dict(gamma=0.20, A=A, hetero=0, n_train=410, phi=phi) for A in (0.15, 0.20, 0.05) for phi in (0.5, 0.8, 0.0, 0.95)]


def designate(conds):
    def find(spec):
        return next(c for c in conds if all(c[k] == v for k, v in spec.items()))
    d1 = d2 = None; log = {"D1": [], "D2": []}
    for pos, spec in enumerate(D1_LIST):
        c = find(spec); t = pop_argmax_tables(c, structure(c)); log["D1"].append((c["idx"], t["A0_flips"]))
        if t["A0_flips"] and d1 is None: d1 = (c["idx"], pos)
    for pos, spec in enumerate(D2_LIST):
        c = find(spec); t = pop_argmax_tables(c, structure(c)); log["D2"].append((c["idx"], not t["A0_flips"]))
        if not t["A0_flips"] and d2 is None: d2 = (c["idx"], pos)
    return d1, d2, log


# ---------------------------------------------------------------- self-test
def cmd_selftest(a):
    t0 = time.time(); g = frozen_gate3a(); rep = {"version": VERSION, "hashes": {"p5_gate3a.py": sha_file(FROZEN_DIR / "p5_gate3a.py"),
                                                                                  "g1r_run.py": sha_file(__file__), "nova_d_calib.py": sha_file(HERE / "nova_d_calib.py")}}
    conds = conditions(); rep["n_conditions"] = len(conds); rng = np.random.default_rng(11)
    # 1/2. mirror == frozen pick on random samples; A1 == mirror to 1e-12 on exact covariances (all grid alphas + inherited fold sizes) and on sample covariances
    mism = 0
    for _ in range(200):
        n = int(rng.integers(60, 600)); O, T, X = rng.normal(size=(n, D_O)), rng.normal(size=(n, D_T)), rng.normal(size=(n, P)); T += 0.5 * O + 0.3 * X * rng.normal(size=P)
        Z = np.column_stack([O, T, X]); S = np.cov((Z - Z.mean(0)) / Z.std(0, ddof=1), rowvar=False, ddof=1)
        mism += int(pick(frozen_algebra(S, alpha_frozen(n))) != g.g1_pick(O, T, X))
    rep["mirror_mismatches_of_200"] = mism
    alphas = sorted({alpha_frozen(n) for n in N_TRAIN} | {D / (D + m) for m in (102, 103, 409, 410, 819, 820)})
    worst, amism, clipact = 0.0, 0, 0
    mats = [standardise_cov(population_cov(structure(c))) for c in conds]
    for _ in range(1000):
        Araw = rng.normal(size=(D, D + 5)); mats.append(standardise_cov(Araw @ Araw.T / (D + 5) + 1e-3 * np.eye(D)))
    for S in mats:
        for al in alphas:
            assert al < 1.0
            sm = frozen_algebra(S, al); C1, v1 = schur(S, (al + EPS) / (1 - al)); s1, cl = score(C1, v1, floor=al / (1 - al), scale=1 - al)
            worst = max(worst, np.max(np.abs(s1 - sm)) / np.max(np.abs(sm))); amism += int(pick(s1) != pick(sm)); clipact += int(cl)
    rep["A1_gate"] = dict(matrices=len(mats), alphas=len(alphas), max_rel_discrepancy=float(worst), argmax_mismatches=amism, clip_active_cases=clipact)
    # 3. Objective A tables and designation
    d1, d2, dlog = designate(conds); rep["designation"] = dict(D1=d1, D2=d2, log=dlog)
    rep["objective_A"] = {c["idx"]: pop_argmax_tables(c, structure(c)) for c in conds}
    # 4. calibration
    rep["calibration_max_abs_A_err"] = float(max(abs(structure(c)["A_or"] - c["A"]) for c in conds))
    # 6. exact-regret identity via the covariance route on designated conditions
    worst_r = 0.0
    for key in [k for k in (d1, d2) if k]:
        c = conds[key[0]]; st = structure(c); Sig = population_cov(st)
        for j in rng.integers(0, P, 10):
            idx = list(range(D_O)) + [D_O + D_T + int(j)]; A_ = Sig[np.ix_(idx, idx)]; Bm = Sig[np.ix_(idx, range(D_O, D_O + D_T))]
            Rj = np.trace(Sig[it, it] - Bm.T @ np.linalg.solve(A_, Bm)); R0 = np.trace(Sig[it, it] - Sig[it, io] @ np.linalg.solve(Sig[io, io], Sig[io, it]))
            worst_r = max(worst_r, abs((R0 - Rj) / R0 - st["G"][j] / st["R0"]))
    rep["exact_regret_identity_max_abs_err"] = float(worst_r)
    # 7. GCV matched-penalty identity and LW vs independent formula; 8. residual identity
    O = rng.normal(size=(300, D_O)); Y = rng.normal(size=(300, D_T + P)); Z = np.column_stack([O, Y]); Zs = (Z - Z.mean(0)) / Z.std(0, ddof=1)
    S = np.cov(Zs, rowvar=False, ddof=1); lam = 0.3; gam = 299 * lam
    lhs = np.linalg.solve(S[io, io] + lam * np.eye(D_O), S[io, D_O:]); rhs = np.linalg.solve(Zs[:, io].T @ Zs[:, io] + gam * np.eye(D_O), Zs[:, io].T @ Zs[:, D_O:])
    rep["gcv_penalty_identity_max_abs_err"] = float(np.max(np.abs(lhs - rhs)))
    rep["lw_vs_direct_formula_abs_err"] = float(abs(ledoit_wolf_delta_direct(Zs) - ledoit_wolf_delta(Zs)))
    try:
        from sklearn.covariance import ledoit_wolf as sk_lw
        _, sk_delta = sk_lw(Zs, assume_centered=True); rep["lw_vs_sklearn_abs_err"] = float(abs(sk_delta - ledoit_wolf_delta(Zs)))
    except Exception as e:  # noqa: BLE001
        rep["lw_vs_sklearn_abs_err"] = f"sklearn unavailable here: {type(e).__name__} (direct-formula check is the gate)"
    # Clopper-Pearson sanity: k = 0 of 1000 -> 1 - 0.05**(1/1000); k = 1000 -> 1
    rep["cp_upper_k0_n1000"] = clopper_pearson_upper(0, 1000); rep["cp_upper_k0_n1000_expected"] = 1 - 0.05 ** (1 / 1000)
    CR, vR = residual(S, lam); CRi, vRi = residual_via_identity(S, lam); C0a, v0a = residual(S, 0.0); C0b, v0b = schur(S, 0.0)
    rep["residual_identity_max_abs_err"] = float(max(np.max(np.abs(CR - CRi)), np.max(np.abs(vR - vRi))))
    rep["forms_coincide_at_lambda0_max_abs_err"] = float(max(np.max(np.abs(C0a - C0b)), np.max(np.abs(v0a - v0b))))
    rep["forms_differ_at_lambda_pos"] = bool(np.max(np.abs(CR - schur(S, lam)[0])) > 1e-6)
    # timing
    c = conds[0]; st = structure(c); t1 = time.perf_counter(); [replicate(c, st, r) for r in range(5)]; rep["sec_per_replicate_n102"] = (time.perf_counter() - t1) / 5
    c = conds[-1]; st = structure(c); t1 = time.perf_counter(); [replicate(c, st, r) for r in range(5)]; rep["sec_per_replicate_n3280"] = (time.perf_counter() - t1) / 5
    rep["ok"] = (mism == 0 and worst < 1e-12 and amism == 0 and clipact == 0 and d1 is not None and d2 is not None
                 and rep["calibration_max_abs_A_err"] < 1e-9 and worst_r < 1e-12 and rep["gcv_penalty_identity_max_abs_err"] < 1e-12
                 and rep["residual_identity_max_abs_err"] < 1e-12 and rep["forms_coincide_at_lambda0_max_abs_err"] < 1e-12 and rep["forms_differ_at_lambda_pos"]
                 and rep["lw_vs_direct_formula_abs_err"] < 1e-10
                 and (not isinstance(rep["lw_vs_sklearn_abs_err"], float) or rep["lw_vs_sklearn_abs_err"] < 1e-10)
                 and abs(rep["cp_upper_k0_n1000"] - rep["cp_upper_k0_n1000_expected"]) < 1e-12)
    rep["elapsed_s"] = time.time() - t0
    OUT.mkdir(exist_ok=True); (OUT / "selftest.json").write_text(json.dumps(rep, indent=1, default=str))
    short = {k: v for k, v in rep.items() if k != "objective_A"}; print(json.dumps(short, indent=1, default=str)); print("SELFTEST", "PASS" if rep["ok"] else "FAIL")
    return 0 if rep["ok"] else 1


# ---------------------------------------------------------------- freeze / run / analyze
def cmd_freeze(a):
    OUT.mkdir(exist_ok=True)
    if RUNLOG.exists(): raise SystemExit("REFUSED: outputs exist")
    if not (HERE / PROTOCOL).exists(): raise SystemExit(f"REFUSED: {PROTOCOL} not found")
    st = OUT / "selftest.json"
    if not st.exists() or not json.loads(st.read_text()).get("ok"): raise SystemExit("REFUSED: selftest has not passed")
    sel = json.loads(st.read_text())["designation"]
    rec = dict(protocol=PROTOCOL, protocol_sha256=sha_file(HERE / PROTOCOL), g1r_run=sha_file(__file__), nova_d_calib=sha_file(HERE / "nova_d_calib.py"),
               p5_gate3a=P5_SHA, salt=SALT, n_train=N_TRAIN, phis=PHIS, gammas=GAMMAS, A_grid=A_GRID, hetero=HETERO, reps=R_REPS,
               lambda_grid_len=len(LAMBDA_GRID), betas=BETAS, margins=(MARGIN_H1, MARGIN_H2), sparse_min=SPARSE_MIN, B_boot=B_BOOT,
               D1=sel["D1"], D2=sel["D2"], n_conditions=len(conditions()), frozen_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    FREEZE.write_text(json.dumps(rec, indent=1)); print(json.dumps(rec, indent=1)); return 0


def check_frozen():
    if not FREEZE.exists(): raise SystemExit("REFUSED: not frozen")
    f = json.loads(FREEZE.read_text())
    for k, p in (("protocol_sha256", HERE / PROTOCOL), ("g1r_run", Path(__file__)), ("nova_d_calib", HERE / "nova_d_calib.py")):
        if sha_file(p) != f[k]: raise SystemExit(f"REFUSED: {p.name} changed since freeze")
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
                fh.write(json.dumps(replicate(c, st, r)) + "\n"); k += 1
            fh.flush(); print(f"[{time.strftime('%H:%M:%S')}] cond {c['idx']} done; total {k}; {(time.time()-t0)/60:.1f} min", flush=True)
    return 0


def clopper_pearson_upper(k, n, conf=0.95):
    """One-sided 95 % Clopper-Pearson upper bound on P(jhat != j*): the p solving P(Bin(n, p) <= k) = 1 - conf,
    found by bisection on the exact binomial CDF (log-space; no scipy)."""
    if k >= n:
        return 1.0
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


def ledoit_wolf_delta_direct(Zc):
    """Independent formula for the same intensity: beta = (1/n^2) sum_t ||z_t z_t' - S_n||_F^2, delta = min(1, beta / ||S_n - mu I||_F^2)."""
    n = Zc.shape[0]; Sn = Zc.T @ Zc / n; mu = np.trace(Sn) / D
    beta_ = sum(np.sum((np.outer(z, z) - Sn) ** 2) for z in Zc) / n ** 2
    gam = np.sum((Sn - mu * np.eye(D)) ** 2)
    return float(min(1.0, beta_ / gam)) if gam > 0 else 1.0


def cmd_analyze(a):
    f = check_frozen(); conds = conditions()
    recs = [json.loads(l) for l in RUNLOG.read_text().splitlines() if l.strip()]; by = {}
    for r in recs: by.setdefault(r["idx"], []).append(r)
    missing = [c["idx"] for c in conds if len(by.get(c["idx"], [])) != R_REPS]
    if missing and not a.partial: raise SystemExit(f"REFUSED: {len(missing)} conditions incomplete")
    rng = np.random.default_rng(seed_for("analyze")); out = {"freeze": f, "conditions": []}
    def contrast(r0, rk, factor):
        q = factor * r0 - rk; m = float(q.mean()); bs = np.array([q[rng.integers(0, len(q), len(q))].mean() for _ in range(B_BOOT)])
        return dict(mean=m, lcb95=float(np.percentile(bs, 5)), ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))])
    for c in conds:
        rs = by.get(c["idx"], []);
        if not rs: continue
        R = {arm: np.array([r["regret"][arm] for r in rs]) for arm in ARMS}; st = structure(c)
        row = dict(c, jstar=st["jstar"], arms={}, contrasts={}, sparse={})
        for arm in ARMS:
            nz = int((R[arm] > 0).sum()); mis = int(sum(r["picks"][arm] != st["jstar"] for r in rs))
            row["arms"][arm] = dict(mean_regret=float(R[arm].mean()), p_optimal=1 - mis / len(rs), nonzero=nz,
                                    cp95_upper_misselection=clopper_pearson_upper(mis, len(rs)))
            row["sparse"][arm] = nz < SPARSE_MIN
        for arm in ARMS[1:]:
            if row["sparse"]["A0"] or row["sparse"][arm]:
                row["contrasts"][arm] = {"assessable": False, "reason": "sparse regret in A0 or arm; Clopper-Pearson bounds reported in arms"}
            else:
                row["contrasts"][arm] = dict(assessable=True, plain=contrast(R["A0"], R[arm], 1.0), Q1=contrast(R["A0"], R[arm], MARGIN_H1), Q2=contrast(R["A0"], R[arm], MARGIN_H2))
        row["A1_max_rel_discrepancy"] = float(max(r["A1_rel_discrepancy"] for r in rs)); row["mirror_all_equal_frozen"] = all(r["mirror_eq_frozen"] for r in rs)
        row["endpoint_hit_rate"] = float(np.mean([r["endpoint_hit"] for r in rs])); row["A4_degenerate_tie_rate"] = float(np.mean([r["A4_degenerate_tie"] for r in rs]))
        row["lambda_c_median"] = float(np.median([r["lambda_c"] for r in rs])); row["delta_lw_median"] = float(np.median([r["delta_lw"] for r in rs]))
        row["kappa_hat_median"] = float(np.median([r["kappa_hat"] for r in rs]))
        out["conditions"].append(row)
    d1, d2 = f["D1"][0], f["D2"][0]; get = lambda i: next(r for r in out["conditions"] if r["idx"] == i)
    c1, c2 = get(d1)["contrasts"]["A2"], get(d2)["contrasts"]["A2"]
    out["decisions"] = dict(
        H1_D1=("not assessable" if not c1["assessable"] else bool(c1["Q1"]["lcb95"] > 0)),
        H2_D2=("not assessable" if not c2["assessable"] else bool(c2["Q2"]["lcb95"] >= 0)),
        grid_max_plain_A2_gain_upward_biased=max((r["contrasts"]["A2"]["plain"]["mean"], r["idx"]) for r in out["conditions"] if r["contrasts"]["A2"]["assessable"]))
    (OUT / "g1r_result.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out["decisions"], indent=1)); return 0


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest); sub.add_parser("freeze").set_defaults(fn=cmd_freeze)
    r = sub.add_parser("run"); r.add_argument("--cond", default=None); r.set_defaults(fn=cmd_run)
    z = sub.add_parser("analyze"); z.add_argument("--partial", action="store_true"); z.set_defaults(fn=cmd_analyze)
    a = ap.parse_args(); sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
