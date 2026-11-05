#!/usr/bin/env python3
"""NOVA-D executable implementation (protocol NOVA_D_PROTOCOL_v0.3.md).

Synthetic test of whether the frozen NOVA-S estimator recovers a known nonlinear acquisition advantage, and
whether near-tie gating improves or obstructs that recovery.

SELECTION CODE IS NOT HERE. Every acquisition decision is made by the frozen modules
  p6_n2.py    v0.1.0  sha256 681ce852...   (NOVA-V / NOVA-U / NOVA-S, cross-fitted orthogonal scores)
  p5_gate3a.py        sha256 4600e2f3...   (G1 plug-in rule, log rule, features, folds)
which are hash-verified before import and called exactly as the NHANES wrapper (nhanes_fa_run.py v0.2.0) called
them. This file only: builds the calibrated synthetic distributions, samples from them, feeds the frozen procedure,
evaluates the resulting picks on an independent evaluation sample, applies the near-tie POLICY (outside the
selector), and aggregates.

Commands
  selftest   hash checks; G1-score mirror == g1_pick exactly; fold-loop mirror == frozen run_block to 1e-9;
             calibration feasible + verified for all conditions; sampler moment checks vs closed forms;
             log-rule branch assertion; timing. Produces NO NOVA-D outcome.
  freeze     records hashes of the protocol, this file and nova_d_calib.py; refuses if any output exists.
  run        refuses unless frozen and hashes match; one run per freeze; resumable and partitionable
             (--cond a:b) because every (condition, replicate) is deterministic and independent.
  analyze    pooled ratios per condition and gamma; kappa; paired bootstrap; prespecified decisions.
"""
from __future__ import annotations
import argparse, hashlib, importlib.util, json, math, os, sys, time
from pathlib import Path
import numpy as np

VERSION = "0.1.0"
PROTOCOL = "NOVA_D_PROTOCOL_v0.3.md"
SALT = "NOVA-D-SYNTH-V03|"
HERE = Path(__file__).resolve().parent
FROZEN_DIR = HERE / "nhanes_fa" / "frozen"
P6_SHA = "681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e"
P5_SHA = "4600e2f3e3bb652affadcc6f6935380d8b9a77a4237e7cb306f056f849478130"
OUT = HERE / "nova_d_out"
FREEZE = OUT / "FREEZE.json"
RUNLOG = OUT / "replicates.jsonl"

# ---- design constants (protocol v0.3 sections 4-6) --------------------------------------------------------------
D_O, D_T, P = 20, 20, 20
PHI, S2, PSI = 0.5, 0.5, 0.25
N_EVAL, R_REPS = 5000, 100
N_GRID = (128, 512, 1024)                 # 2048 replaced prospectively on measured cost (protocol s13)
GAMMAS = (0.02, 0.05, 0.10, 0.25, math.inf)
BASE_SHARE_FACTOR = 2.0                   # c = 2 * max_k condvar_k: equal marginal target variances (s4)
PRIMARY = ("M2", 0.10, 512, 0.15)
POSITIVE_CONTROL = ("M2", 0.40, 1024, 0.15)
B_BOOT = 2000
MOMENT_N, MOMENT_TOL = 400_000, 0.02

sys.path.insert(0, str(HERE))
import nova_d_calib as CAL                 # noqa: E402  (deterministic calibration; hashed with this file)


def sha_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------- frozen module
_N2 = None


def frozen():
    """Load p6_n2.py after verifying both hashes; p6_n2 re-verifies p5_gate3a itself."""
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


# ---------------------------------------------------------------- seeds and conditions
def seed_for(*parts) -> int:
    h = hashlib.sha256((SALT + "|".join(str(p) for p in parts)).encode()).hexdigest()
    return int(h[:15], 16)                                  # < 2**60


def conditions():
    """The 72 conditions in a fixed order: (index, mech, rank, gamma_rel, rho_pop, n, A)."""
    out, i = [], 0
    for mech, rank, gam, rhos in CAL.MECHS:
        for rho in rhos:
            for n in N_GRID:
                for A in CAL.A_GRID:
                    out.append(dict(idx=i, mech=mech, rank=rank, gam=gam, rho=rho, n=n, A=A)); i += 1
    return out


# ---------------------------------------------------------------- structure (fixed per condition) and sampler
def unit_rows(rng, k, d):
    v = rng.normal(size=(k, d)); return v / np.linalg.norm(v, axis=1, keepdims=True)


def structure(c):
    """Calibrated, fully specified distribution for one condition. Deterministic in the condition's structure seed."""
    cell = CAL.build(c["mech"], c["rank"], c["gam"], c["A"], c["rho"], PSI)
    if not cell.feasible:
        raise SystemExit(f"REFUSED: condition {c} infeasible")
    rng = np.random.default_rng(seed_for(c["mech"], c["rho"], c["n"], c["A"], "structure"))
    r, H, s_eps2 = cell.r, cell.H, cell.sigma2
    G = np.array([CAL.G_at(k, r, c["gam"]) for k in range(1, P + 1)])         # G-rank k == candidate index k-1
    nl = c["rank"] - 1 if c["rho"] > 0 else -1                                  # designated nonlinear candidate
    Hv = np.zeros(P);  Hv[nl] = H if nl >= 0 else 0.0
    Sig = PHI ** np.abs(np.subtract.outer(np.arange(D_O), np.arange(D_O)))
    L = np.linalg.cholesky(Sig)
    def scaled_dirs(k, var):                       # rows w with w' Sig w = var
        w = rng.normal(size=(k, D_O)); q = np.einsum("ij,jk,ik->i", w, Sig, w)
        return w * np.sqrt(var / q)[:, None]
    a = scaled_dirs(P, 1.0 - S2)                                                  # Var(a_j'O) = 1 - s^2
    lam = (unit_rows(rng, P, D_T) * np.sqrt(G / S2)[:, None]).T                   # (D_T, P): sum_k lam^2 = G_j / s^2
    nu = np.zeros((D_T, P))
    if nl >= 0:
        nu[:, nl] = unit_rows(rng, 1, D_T)[0] * np.sqrt(H)                        # sum_k nu^2 = H
    condvar = (lam ** 2).sum(axis=1) * S2 + (nu ** 2).sum(axis=1) + s_eps2       # Var(T_k | O), per k
    cc = BASE_SHARE_FACTOR * condvar.max()
    b = np.vstack([scaled_dirs(1, cc - condvar[k])[0] for k in range(D_T)])      # Var(b_k'O) = c - condvar_k
    V = G + Hv
    jG, jstar = int(np.argmax(G)), int(np.argmax(V))
    R0 = G.sum() + Hv.sum() + D_T * s_eps2
    A_or = V[jG] / R0; rho_or = (V[jstar] / V[jG] - 1.0) / (1.0 - A_or)
    return dict(cond=c, r=r, H=H, s_eps2=s_eps2, L=L, Sig=Sig, a=a, lam=lam, nu=nu, b=b, nl=nl, G=G, Hv=Hv, V=V,
                jG=jG, jstar=jstar, R0=R0, A_or=A_or, rho_or=rho_or, c_marg=cc, condvar=condvar)


def hermite2(z):
    return (z * z - 1.0) / math.sqrt(2.0)


def sample(st, n, seed):
    rng = np.random.default_rng(seed)
    O = rng.normal(size=(n, D_O)) @ st["L"].T
    u = rng.normal(scale=math.sqrt(S2), size=(n, P))
    X = O @ st["a"].T + u
    h = hermite2(u / math.sqrt(S2))
    T = O @ st["b"].T + u @ st["lam"].T + h @ st["nu"].T + rng.normal(scale=math.sqrt(st["s_eps2"]), size=(n, D_T))
    return O, X, T, u


def moment_check(st, n=MOMENT_N, seed=0, tol=MOMENT_TOL, batches=20, zmax=4.0):
    """Owner's safeguard (v0.3 s11.6): the sampler must instantiate the calibrated distribution.
    Independently estimated population moments are compared with the closed forms by (i) a z-test from batch
    means (Monte Carlo SE) and (ii) an absolute cap at the experiment's own scale (G_(1) = 1 for G and H; tol for
    A and rho). ASSERT-ONLY: no calibrated parameter is changed in response."""
    O, X, T, u = sample(st, n, seed)
    Rt = T - O @ st["b"].T                                   # T - E[T|O], exact residual (b known)
    h = hermite2(u / math.sqrt(S2))
    rep = {}
    orth = (u * h).mean(axis=0); orth_se = (u * h).std(axis=0, ddof=1) / math.sqrt(n)
    rep["orthogonality_max_abs_z"] = float(np.max(np.abs(orth / orth_se)))
    Ghat = np.array([((Rt * u[:, [j]]).mean(axis=0) ** 2).sum() / u[:, j].var(ddof=1) for j in range(P)])
    Vrow = np.stack([((u[:, [j]] * st["lam"][:, j] + h[:, [j]] * st["nu"][:, j]) ** 2).sum(axis=1) for j in range(P)], axis=1)
    Vhat = Vrow.mean(axis=0)
    R0row = (Rt ** 2).sum(axis=1)
    rep["G_abs_err_max_over_G1"] = float(np.max(np.abs(Ghat - st["G"])))          # G_(1) == 1 by scale
    rep["H_abs_err_max_over_G1"] = float(np.max(np.abs((Vhat - Ghat) - st["Hv"])))
    # batch means for A and rho (Monte Carlo SE of the ratio estimates)
    bsz = n // batches; As, rhos, R0s = [], [], []
    for b in range(batches):
        sl = slice(b * bsz, (b + 1) * bsz)
        Vb = Vrow[sl].mean(axis=0); R0b = R0row[sl].mean()
        Ab = Vb[st["jG"]] / R0b; As.append(Ab); R0s.append(R0b)
        rhos.append((Vb[st["jstar"]] / Vb[st["jG"]] - 1.0) / (1.0 - Ab))
    As, rhos, R0s = map(np.array, (As, rhos, R0s))
    def ztest(vals, target):
        m, se = vals.mean(), vals.std(ddof=1) / math.sqrt(len(vals))
        return float(m), float(se), float(abs(m - target) / se) if se > 0 else 0.0
    rep["A_est"], rep["A_se"], rep["A_z"] = ztest(As, st["cond"]["A"]); rep["A_target"] = float(st["cond"]["A"])
    rep["rho_est"], rep["rho_se"], rep["rho_z"] = ztest(rhos, st["cond"]["rho"]); rep["rho_target"] = float(st["cond"]["rho"])
    rep["R0_est"], rep["R0_se"], rep["R0_z"] = ztest(R0s, st["R0"]); rep["R0_target"] = float(st["R0"])
    rep["marginal_var_rel_spread"] = float(np.max(np.abs(T.var(axis=0, ddof=1) / st["c_marg"] - 1.0)))
    rep["argmaxV_is_designated"] = bool(st["jstar"] == (st["nl"] if st["nl"] >= 0 else st["jG"]))
    rep["negatives_in_every_feature"] = bool((np.column_stack([O, X]) < 0).any(axis=0).all())
    ok = (rep["orthogonality_max_abs_z"] < zmax
          and rep["G_abs_err_max_over_G1"] < tol and rep["H_abs_err_max_over_G1"] < tol
          and rep["A_z"] < zmax and abs(rep["A_est"] - st["cond"]["A"]) < tol
          and rep["rho_z"] < zmax and abs(rep["rho_est"] - st["cond"]["rho"]) < tol
          and rep["R0_z"] < zmax and abs(rep["R0_est"] / st["R0"] - 1.0) < tol
          and rep["marginal_var_rel_spread"] < tol and rep["argmaxV_is_designated"]
          and rep["negatives_in_every_feature"])
    rep["ok"] = bool(ok)
    return rep


# ---------------------------------------------------------------- G1 score vector (diagnostic; the pick is still g1_pick)
def g1_scores(O, T, Pm):
    """Verbatim body of the frozen p5_gate3a.g1_pick, returning the score vector instead of its argmax.
    Self-test asserts argmax(g1_scores) == g1_pick exactly. Used only for Gamma_hat_obs (the near-tie policy)."""
    Z = np.column_stack([O, T, Pm]).astype(float)
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


def gamma_obs(score):
    s = np.sort(score)[::-1]
    return float((s[0] - s[1]) / s[0]) if s[0] > 0 else math.inf     # s[0] == 0: no near-tie; policy uses G1


# ---------------------------------------------------------------- the fold loop: frozen decisions, external losses
def standardise(Ttr, Tev):
    n2 = frozen()
    sd = Ttr.std(axis=0, ddof=1); sd = np.where(np.isfinite(sd) & (sd > 0), sd, 0.0)
    floor = n2.FLOOR * np.median(sd[sd > 0]) if (sd > 0).any() else 1.0
    den = np.maximum(sd, floor)
    return Ttr / den, Tev / den


def fold_loop(O, X, Y, subjects, seed, Oev=None, Xev=None, Yev=None, st=None, require_identity=True):
    """Same calls and same decisions as p6_n2.run_block(family='rbf') / the NHANES mirror. Returns, per outer fold:
    the picks, claim, t, Gamma_hat_obs, within-sample held-out losses (for the equivalence test) and, when an
    evaluation sample is given, external pooled losses per distinct candidate and for the base panel."""
    n2 = frozen(); g = n2.gate3a()
    fold = g.folds_for(subjects, seed)
    n, m = len(Y), X.shape[1]
    arms = ("G1",) + n2.VARIANTS + ("BASE",)
    sq = {a: np.full(n, np.nan) for a in arms}
    rows = []
    for k in range(n2.K_OUT):
        tr, te = fold != k, fold == k
        if te.sum() == 0:
            continue
        if require_identity:                                      # NOVA-D data only: every feature has negatives
            lr = g.LogRule().fit(np.column_stack([O[tr], X[tr]]))
            if lr.neg.sum() != O.shape[1] + X.shape[1] or lr.drop.any():
                raise SystemExit("REFUSED: log rule left the identity branch on synthetic data (protocol s4)")
        jg = g.g1_pick(O[tr], Y[tr], X[tr])
        gam = gamma_obs(g1_scores(O[tr], Y[tr], X[tr]))
        psi, best = n2.inner_scores("rbf", O[tr], X[tr], Y[tr], n2.SEED + k)
        ok = ~np.isnan(psi).any(axis=1); psi = psi[ok]
        Vh = psi.mean(axis=0); se = psi.std(axis=0, ddof=1) / np.sqrt(len(psi))
        jv = int(np.flatnonzero(Vh == Vh.max())[0])
        u = Vh - n2.LAMBDA_U * se; ju = int(np.flatnonzero(u == u.max())[0])
        claim, tstat = False, 0.0
        if jv != jg:
            dd = psi[:, jv] - psi[:, jg]; sdd = dd.std(ddof=1) / np.sqrt(len(dd))
            tstat = float(dd.mean() / sdd) if sdd > 0 else math.inf
            claim = bool(dd.mean() > n2.CLAIM_SE * sdd)
        js = jv if claim else jg
        pk = {"G1": jg, "NOVA-V": jv, "NOVA-U": ju, "NOVA-S": js}
        cache = {}
        for a, j in pk.items():                                   # within-sample held-out (equivalence test)
            if j not in cache:
                Pp = n2.predict_grid("rbf", O[tr], O[te], Y[tr], X[tr, j], X[te, j])[best[j + 1]]
                cache[j] = ((Pp - Y[te]) ** 2).sum(axis=1)
            sq[a][te] = cache[j]
        P0 = n2.predict_grid("rbf", O[tr], O[te], Y[tr])[best[0]]
        sq["BASE"][te] = ((P0 - Y[te]) ** 2).sum(axis=1)
        row = dict(fold=k, n_tr=int(tr.sum()), picks=pk, claim=claim, t=tstat, gamma_obs=gam)
        if Oev is not None:                                       # external evaluation: same hypers, new rows
            ext = {}
            for j in set(pk.values()):
                Pp = n2.predict_grid("rbf", O[tr], Oev, Y[tr], X[tr, j], Xev[:, j])[best[j + 1]]
                ext[int(j)] = float(((Pp - Yev) ** 2).sum())
            P0e = n2.predict_grid("rbf", O[tr], Oev, Y[tr])[best[0]]
            row["ext_loss"] = ext; row["ext_base"] = float(((P0e - Yev) ** 2).sum()); row["n_ev"] = int(len(Yev))
            if st is not None:
                row["oracle"] = {a: float(st["V"][st["jstar"]] - st["V"][j]) / st["R0"] for a, j in pk.items()}
        rows.append(row)
    return rows, sq


def policy_pick(row, gamma):
    return row["picks"]["NOVA-S"] if row["gamma_obs"] <= gamma else row["picks"]["G1"]


def replicate(c, st, rep):
    O, X, T, _ = sample(st, c["n"], seed_for(c["mech"], c["rho"], c["n"], c["A"], rep, "design"))
    Oev, Xev, Tev, _ = sample(st, N_EVAL, seed_for(c["mech"], c["rho"], c["n"], c["A"], rep, "evaluation"))
    Y, Yev = standardise(T, Tev)
    subjects = [str(i) for i in range(c["n"])]
    rows, _ = fold_loop(O, X, Y, subjects, seed_for(c["mech"], c["rho"], c["n"], c["A"], rep, "folds"),
                        Oev, Xev, Yev, st)
    # pooled external loss sums over folds, per arm and per gamma policy (paired: same evaluation rows)
    sums = {"BASE": 0.0, "G1": 0.0, "NOVA-V": 0.0, "NOVA-U": 0.0, "NOVA-S": 0.0}
    pol = {str(gm): 0.0 for gm in GAMMAS}
    for r in rows:
        sums["BASE"] += r["ext_base"]
        for a in ("G1", "NOVA-V", "NOVA-U", "NOVA-S"):
            sums[a] += r["ext_loss"][r["picks"][a]]
        for gm in GAMMAS:
            pol[str(gm)] += r["ext_loss"][policy_pick(r, gm)]
    dep = sum(int(r["picks"]["NOVA-S"] != r["picks"]["G1"]) for r in rows)
    ben = sum(int(r["ext_loss"][r["picks"]["NOVA-S"]] < r["ext_loss"][r["picks"]["G1"]]) for r in rows
              if r["picks"]["NOVA-S"] != r["picks"]["G1"])
    opp = sum(int(st["jstar"] != r["picks"]["G1"]) for r in rows)              # a beneficial departure existed
    excl = {str(gm): sum(int(st["jstar"] != r["picks"]["G1"] and r["gamma_obs"] > gm) for r in rows)
            for gm in GAMMAS}
    return dict(idx=c["idx"], rep=rep, folds=len(rows), sums=sums, policy=pol, departures=dep, beneficial=ben,
                opportunities=opp, screen_excluded=excl,
                picks=[r["picks"] for r in rows], t=[r["t"] for r in rows], gamma_obs=[r["gamma_obs"] for r in rows],
                oracle_regret=[r["oracle"] for r in rows], jstar=st["jstar"], jG=st["jG"])


# ---------------------------------------------------------------- self-test (no NOVA-D outcome)
def cmd_selftest(a):
    t0 = time.time(); n2 = frozen(); g = n2.gate3a(); rep = {"version": VERSION}
    rep["hashes"] = {"p6_n2.py": sha_file(FROZEN_DIR / "p6_n2.py"), "p5_gate3a.py": sha_file(FROZEN_DIR / "p5_gate3a.py"),
                     "nova_d_run.py": sha_file(__file__), "nova_d_calib.py": sha_file(HERE / "nova_d_calib.py")}
    rng = np.random.default_rng(1)
    # 1. G1 score mirror: argmax equals the frozen pick exactly, on random blocks of several sizes
    mism = 0
    for i in range(200):
        n = int(rng.integers(40, 400))
        O, Tt, Pm = rng.normal(size=(n, D_O)), rng.normal(size=(n, D_T)), rng.normal(size=(n, P))
        Tt += O[:, :D_T] * 0.5 + Pm * rng.normal(size=P) * 0.3
        mism += int(int(np.argmax(g1_scores(O, Tt, Pm))) != g.g1_pick(O, Tt, Pm))
    rep["g1_score_mirror_mismatches_of_200"] = mism
    # 2. fold-loop mirror == frozen run_block on p6_n2's own synthetic blocks (within-sample losses)
    maxdiff, pickdiff = 0.0, 0
    for i, kind in enumerate(("null", "quadratic", "null", "quadratic")):
        O, X, Y = n2.synth(np.random.default_rng(100 + i), 160, kind)
        subj = [str(k) for k in range(len(Y))]
        ref = n2.run_block(O, X, Y, subj, 777 + i, "rbf")
        rows, sq = fold_loop(O, X, Y, subj, 777 + i, require_identity=False)   # frozen lognormal blocks: log branch expected
        for arm in ("G1", "NOVA-S", "NOVA-V", "NOVA-U"):
            ours = float(np.nansum(sq[arm])); theirs = float(ref["sum"][arm]) if "sum" in ref else float(np.nansum(ref["sq"][arm]))
            maxdiff = max(maxdiff, abs(ours - theirs) / max(abs(theirs), 1e-9))
        if "picks" in ref:
            pickdiff += sum(int(rp != {k: v for k, v in r["picks"].items()}) for rp, r in zip(ref["picks"], rows))
    rep["fold_loop_equivalence_max_rel_diff"] = maxdiff; rep["fold_loop_pick_mismatches"] = pickdiff
    # 3. calibration: all conditions feasible and verified; moment checks on four designated conditions
    conds = conditions(); rep["n_conditions"] = len(conds)
    infeas = [c["idx"] for c in conds if not CAL.build(c["mech"], c["rank"], c["gam"], c["A"], c["rho"], PSI).feasible]
    rep["infeasible_conditions"] = infeas
    rep["moment_checks"] = {}
    for lab, key in (("primary", PRIMARY), ("positive_control", POSITIVE_CONTROL),
                     ("aim_control_M4", ("M4", 0.10, 512, 0.15)), ("instability_control_M3", ("M3", 0.0, 512, 0.15)),
                     ("null_M1", ("M1", 0.0, 128, 0.05)), ("M2_rho05_A20", ("M2", 0.05, 128, 0.20))):
        c = next(c for c in conds if (c["mech"], c["rho"], c["n"], c["A"]) == key)
        rep["moment_checks"][lab] = moment_check(structure(c), seed=seed_for("moment", lab))
    # 4. timing of the frozen solve (no NOVA-D distribution involved)
    rep["timing_s"] = {}
    for n in N_GRID:
        ntr = int(0.8 * 0.8 * n)
        O = rng.normal(size=(ntr, D_O)); Y = rng.normal(size=(ntr, D_T)); x = rng.normal(size=ntr)
        t = time.perf_counter(); n2.predict_grid("rbf", O, O[: ntr // 4], Y, x, x[: ntr // 4]); dt = time.perf_counter() - t
        rep["timing_s"][str(n)] = {"inner_predict_grid": dt, "projected_per_replicate": 25 * 21 * dt * 1.3}
    rep["ok"] = (mism == 0 and maxdiff < 1e-9 and pickdiff == 0 and not infeas
                 and all(m["ok"] for m in rep["moment_checks"].values()))
    rep["elapsed_s"] = time.time() - t0
    OUT.mkdir(exist_ok=True); (OUT / "selftest.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1)); print("SELFTEST", "PASS" if rep["ok"] else "FAIL")
    return 0 if rep["ok"] else 1


# ---------------------------------------------------------------- freeze / run / analyze
def cmd_freeze(a):
    OUT.mkdir(exist_ok=True)
    if RUNLOG.exists():
        raise SystemExit("REFUSED: outputs exist; a freeze must precede the first outcome")
    prot = HERE / PROTOCOL
    if not prot.exists():
        raise SystemExit(f"REFUSED: {PROTOCOL} not found")
    st = OUT / "selftest.json"
    if not st.exists() or not json.loads(st.read_text()).get("ok"):
        raise SystemExit("REFUSED: selftest has not passed")
    rec = {"protocol": PROTOCOL, "protocol_sha256": sha_file(prot), "nova_d_run.py": sha_file(__file__),
           "nova_d_calib.py": sha_file(HERE / "nova_d_calib.py"), "p6_n2.py": P6_SHA, "p5_gate3a.py": P5_SHA,
           "salt": SALT, "n_grid": N_GRID, "A_grid": CAL.A_GRID, "psi": PSI, "reps": R_REPS, "n_eval": N_EVAL,
           "gammas": [str(x) for x in GAMMAS], "primary": PRIMARY, "positive_control": POSITIVE_CONTROL,
           "n_conditions": len(conditions()), "frozen_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    FREEZE.write_text(json.dumps(rec, indent=1)); print(json.dumps(rec, indent=1)); return 0


def check_frozen():
    if not FREEZE.exists():
        raise SystemExit("REFUSED: not frozen")
    f = json.loads(FREEZE.read_text())
    for k, p in (("protocol_sha256", HERE / PROTOCOL), ("nova_d_run.py", Path(__file__)), ("nova_d_calib.py", HERE / "nova_d_calib.py")):
        if sha_file(p) != f[k]:
            raise SystemExit(f"REFUSED: {p.name} changed since freeze")
    return f


def done_pairs():
    s = set()
    if RUNLOG.exists():
        for line in RUNLOG.read_text().splitlines():
            if line.strip():
                d = json.loads(line); s.add((d["idx"], d["rep"]))
    return s


def cmd_run(a):
    check_frozen(); conds = conditions()
    lo, hi = (0, len(conds)) if a.cond is None else tuple(int(x) for x in a.cond.split(":"))
    done = done_pairs(); t0 = time.time(); k = 0
    with RUNLOG.open("a") as fh:
        for c in conds[lo:hi]:
            st = structure(c)
            for rep in range(R_REPS):
                if (c["idx"], rep) in done:
                    continue
                t1 = time.time(); r = replicate(c, st, rep); r["sec"] = time.time() - t1
                fh.write(json.dumps(r) + "\n"); fh.flush(); k += 1
                if k % 5 == 0:
                    print(f"[{time.strftime('%H:%M:%S')}] cond {c['idx']} rep {rep} {r['sec']:.0f}s  done {k}  "
                          f"elapsed {(time.time()-t0)/3600:.2f}h", flush=True)
    print("run segment complete"); return 0


def pooled(recs, arm_key):
    """arm_key(rec) -> loss sum of the arm; returns (Delta, A, rho) from pooled sums over the given replicates."""
    L0 = sum(r["sums"]["BASE"] for r in recs); LG = sum(r["sums"]["G1"] for r in recs); LA = sum(arm_key(r) for r in recs)
    A = 1.0 - LG / L0; D = 1.0 - LA / LG
    return D, A, (D / A if A > 0 else float("nan"))


def cmd_analyze(a):
    f = check_frozen(); conds = conditions()
    recs = [json.loads(l) for l in RUNLOG.read_text().splitlines() if l.strip()]
    by = {}
    for r in recs:
        by.setdefault(r["idx"], []).append(r)
    missing = [c["idx"] for c in conds if len(by.get(c["idx"], [])) != R_REPS]
    if missing and not a.partial:
        raise SystemExit(f"REFUSED: {len(missing)} conditions incomplete; analysis runs only on the full grid")
    rng = np.random.default_rng(seed_for("analyze"))
    out = {"freeze": f, "conditions": []}
    for c in conds:
        rs = by.get(c["idx"], [])
        if not rs:
            continue
        row = dict(c)
        arms = {"NOVA-S": lambda r: r["sums"]["NOVA-S"], "NOVA-V": lambda r: r["sums"]["NOVA-V"], "NOVA-U": lambda r: r["sums"]["NOVA-U"]}
        for gm in GAMMAS:
            arms[f"policy@{gm}"] = (lambda g_: (lambda r: r["policy"][str(g_)]))(gm)
        row["arms"] = {}
        idx = np.arange(len(rs)); boot = [rng.integers(0, len(rs), len(rs)) for _ in range(B_BOOT)]
        for name, key in arms.items():
            D, A, rho = pooled(rs, key)
            bs = np.array([pooled([rs[i] for i in b], key) for b in boot])
            e = dict(Delta=D, A=A, rho=rho, Delta_lcb95=float(np.percentile(bs[:, 0], 5)),
                     rho_lcb95=float(np.percentile(bs[:, 2], 5)), rho_ci95=[float(np.percentile(bs[:, 2], 2.5)), float(np.percentile(bs[:, 2], 97.5))])
            if c["rho"] > 0:
                e["kappa"] = rho / c["rho"]; e["kappa_ci95"] = [x / c["rho"] for x in e["rho_ci95"]]
            row["arms"][name] = e
        # screened vs unscreened, paired
        row["paired_screen_gain"] = {}
        for gm in GAMMAS[:-1]:
            diffs = []
            for b in boot:
                sub = [rs[i] for i in b]
                diffs.append(pooled(sub, arms[f"policy@{gm}"])[2] - pooled(sub, arms["NOVA-S"])[2])
            d0 = row["arms"][f"policy@{gm}"]["rho"] - row["arms"]["NOVA-S"]["rho"]
            row["paired_screen_gain"][str(gm)] = dict(diff=d0, ci95=[float(np.percentile(diffs, 2.5)), float(np.percentile(diffs, 97.5))])
        folds = sum(r["folds"] for r in rs)
        row["p_dep"] = sum(r["departures"] for r in rs) / folds
        row["p_ben_given_dep"] = (sum(r["beneficial"] for r in rs) / max(1, sum(r["departures"] for r in rs)))
        row["opportunity_rate"] = sum(r["opportunities"] for r in rs) / folds
        row["screen_exclusion_of_opportunities"] = {str(gm): (sum(r["screen_excluded"][str(gm)] for r in rs) / max(1, sum(r["opportunities"] for r in rs))) for gm in GAMMAS}
        row["oracle_regret_mean"] = {arm: float(np.mean([o[arm] for r in rs for o in r["oracle_regret"]])) for arm in ("G1", "NOVA-S", "NOVA-V")}
        row["pop_approx_gap"] = float((structure(c)["V"][structure(c)["jstar"]] - structure(c)["V"][structure(c)["jG"]]) / structure(c)["R0"])
        out["conditions"].append(row)
    # prespecified decisions
    get = lambda key: next((r for r in out["conditions"] if (r["mech"], r["rho"], r["n"], r["A"]) == key), None)
    pc, pr = get(POSITIVE_CONTROL), get(PRIMARY)
    out["decisions"] = {
        "positive_control_regime_viability": (pc is not None and pc["arms"]["NOVA-S"]["rho_lcb95"] > 0.10),
        "primary_recoverability": (pr is not None and pr["arms"]["NOVA-S"]["rho_lcb95"] > 0.0),
        "validity_M1_M3_breach": [r["idx"] for r in out["conditions"] if r["mech"] in ("M1", "M3") and r["arms"]["NOVA-S"]["rho"] >= 0.10],
        "near_tie_effective_conditions": [(r["idx"], gm) for r in out["conditions"] for gm, v in r["paired_screen_gain"].items() if v["ci95"][0] > 0],
        "grid_max_rho_NOVA_S_upward_biased": max((r["arms"]["NOVA-S"]["rho"], r["idx"]) for r in out["conditions"]),
    }
    (OUT / "nova_d_result.json").write_text(json.dumps(out, indent=1)); print(json.dumps(out["decisions"], indent=1)); return 0


def main():
    ap = argparse.ArgumentParser(); sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    sub.add_parser("freeze").set_defaults(fn=cmd_freeze)
    r = sub.add_parser("run"); r.add_argument("--cond", default=None, help="a:b slice of condition indices"); r.set_defaults(fn=cmd_run)
    z = sub.add_parser("analyze"); z.add_argument("--partial", action="store_true"); z.set_defaults(fn=cmd_analyze)
    a = ap.parse_args(); sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
