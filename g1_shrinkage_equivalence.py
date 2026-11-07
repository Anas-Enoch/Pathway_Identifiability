#!/usr/bin/env python3
"""Exact equivalence: frozen joint shrinkage == matched ridge-conditioned scoring (pre-freeze derivation check).

Frozen rule on standardised Z = [O, T, P] (mu = tr S / d = 1):
    S_a = (1 - a) S + a I ;  C_a = S_a,TP - S_a,TO (S_a,OO + eps I)^-1 S_a,OP ;  v_a,j = [S_a,PP - S_a,PO (S_a,OO + eps I)^-1 S_a,OP]_jj
    score_a,j = sum_k C_a,kj^2 / v_a,j
Claim (derived):  with lam = (a + eps) / (1 - a),
    C_a = (1 - a) * C_ridge(lam),            C_ridge(lam) = S_TP - S_TO (S_OO + lam I)^-1 S_OP
    v_a,j = (1 - a) * v_ridge,j(lam) + a,    v_ridge,j(lam) = [S_PP - S_PO (S_OO + lam I)^-1 S_OP]_jj
    score_a,j = (1 - a) * sum_k C_ridge,kj^2 / (v_ridge,j + a / (1 - a))
The factor (1 - a) does not affect the argmax. Two separable effects: ridge lam on the O-regression (under-conditioning,
leaks unconditional covariance), and an additive variance floor a / (1 - a) (reorders only when conditional variances
differ across candidates). This script tests the identity numerically on (1) random covariances and (2) the exact
NOVA-D population covariances, and isolates the two effects. Read-only; nothing frozen is touched or changed.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nova_d_run as R                    # noqa: E402
import nova_d_exact_ordering as X         # noqa: E402

EPS = 1e-6
o, t, p = R.D_O, R.D_T, R.P
io, it, ip = slice(0, o), slice(o, o + t), slice(o + t, o + t + p)


def standardise(Sigma):
    sd = np.sqrt(np.diag(Sigma)); return Sigma / np.outer(sd, sd)


def frozen_joint(S, a):
    """Exactly the frozen algebra (after standardisation), returning (C, v, score)."""
    d = S.shape[0]; Sa = (1 - a) * S + a * (np.trace(S) / d) * np.eye(d)
    sol = np.linalg.solve(Sa[io, io] + EPS * np.eye(o), Sa[io, ip])
    C = Sa[it, ip] - Sa[it, io] @ sol; v = np.diag(Sa[ip, ip] - Sa[ip, io] @ sol)
    return C, v, (C * C).sum(0) / np.maximum(v, 1e-12)


def ridge_conditioned(S, lam, floor):
    """Decoupled form: ridge lam on the O-regression, additive floor on the residual variance."""
    sol = np.linalg.solve(S[io, io] + lam * np.eye(o), S[io, ip])
    C = S[it, ip] - S[it, io] @ sol; v = np.diag(S[ip, ip] - S[ip, io] @ sol) + floor
    return C, v, (C * C).sum(0) / np.maximum(v, 1e-12)


def matched(S, a):
    """The ridge-conditioned estimator matched to joint shrinkage a; returns scores scaled by (1 - a)."""
    lam, floor = (a + EPS) / (1 - a), a / (1 - a)
    C, v, sc = ridge_conditioned(S, lam, floor)
    return (1 - a) * C, (1 - a) * v, (1 - a) * sc


def main():
    rng = np.random.default_rng(7); worst = 0.0; argmax_mism = 0; ntest = 0
    # (1) random SPD covariances of the right block shape, standardised
    for i in range(300):
        d = o + t + p; A = rng.normal(size=(d, d + 5)); S = standardise(A @ A.T / (d + 5) + 1e-3 * np.eye(d))
        for a in (0.05, 0.0683, 0.1279, 0.3704, 0.6):
            C1, v1, s1 = frozen_joint(S, a); C2, v2, s2 = matched(S, a)
            rel = max(np.max(np.abs(C1 - C2)) / np.max(np.abs(C1)), np.max(np.abs(v1 - v2)) / np.max(np.abs(v1)),
                      np.max(np.abs(s1 - s2)) / np.max(np.abs(s1)))
            worst = max(worst, rel); argmax_mism += int(np.argmax(s1) != np.argmax(s2)); ntest += 1
    print(f"(1) random covariances: {ntest} (S, a) pairs; max relative discrepancy {worst:.3e}; argmax mismatches {argmax_mism}")
    # (2) exact NOVA-D population covariances, operative alphas
    worst2 = 0.0; mism2 = 0; n2 = 0
    for c in R.conditions():
        S = standardise(X.population_cov(R.structure(c)))
        for a in (0.05, 0.0683, 0.1279, 0.3704):
            _, _, s1 = frozen_joint(S, a); _, _, s2 = matched(S, a)
            worst2 = max(worst2, np.max(np.abs(s1 - s2)) / np.max(np.abs(s1))); mism2 += int(np.argmax(s1) != np.argmax(s2)); n2 += 1
    print(f"(2) NOVA-D exact covariances: {n2} pairs; max relative discrepancy {worst2:.3e}; argmax mismatches {mism2}")
    # (3) isolate the two effects on the NOVA-D family at the operative alphas: ridge only vs floor only
    print("(3) which effect flips the argmax (vs Gaussian argmax at lam = 0, floor = 0)?")
    print(f"    {'alpha':>7} {'flips: both':>12} {'ridge only':>11} {'floor only':>11}  of 72")
    for a in (0.05, 0.0683, 0.1279, 0.3704):
        lam, floor = (a + EPS) / (1 - a), a / (1 - a); fb = fr = ff = 0
        for c in R.conditions():
            S = standardise(X.population_cov(R.structure(c))); g = int(np.argmax(ridge_conditioned(S, 0.0, 0.0)[2]))
            fb += int(np.argmax(ridge_conditioned(S, lam, floor)[2]) != g)
            fr += int(np.argmax(ridge_conditioned(S, lam, 0.0)[2]) != g)
            ff += int(np.argmax(ridge_conditioned(S, 0.0, floor)[2]) != g)
        print(f"    {a:7.4f} {fb:>12} {fr:>11} {ff:>11}")
    ok = worst < 1e-10 and worst2 < 1e-10 and argmax_mism == 0 and mism2 == 0
    print("EQUIVALENCE", "HOLDS" if ok else "FAILS"); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
