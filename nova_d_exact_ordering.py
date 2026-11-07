#!/usr/bin/env python3
"""Exact population ordering of the frozen G1 rule on the NOVA-D family (read-only; no sampling; no simulation).

For each of the 72 frozen NOVA-D conditions, builds the EXACT population covariance of Z = [O, T, P] from the
condition's structure (nova_d_run.structure, which is deterministic in the frozen salt) and applies the frozen
p5_gate3a.g1_pick algebra to it: standardise to unit variances, shrink S <- (1-a) S + a (tr S / d) I, condition on O
with the frozen 1e-6 ridge on S_OO, score_j = sum_k C_kj^2 / v_j. Reports the resulting population argmax for the
frozen rule at its alpha floor (0.05) and for the same rule at alpha = 0, against the Gaussian ranking G_j and the
nonlinear optimum j*.

Population covariance (Hermite term h(u/s) has E[u h] = 0, Var h = 1, independent across candidates):
    Sigma_OO = Sig                       Sigma_OP = Sig a'            Sigma_OT = Sig b'
    Sigma_PP = a Sig a' + s^2 I          Sigma_TP = b Sig a' + s^2 lam
    Sigma_TT = b Sig b' + s^2 lam lam' + nu nu' + sigma_eps^2 I
This script changes nothing: it imports the frozen implementation for its structure definitions only and writes a
single JSON record plus a Markdown table.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nova_d_run as R   # noqa: E402  (structure(), conditions(), constants; nothing executed)

RIDGE = 1e-6             # the frozen g1_pick's ridge on the observed block, retained exactly


def population_cov(st):
    Sig, a, b, lam, nu = st["Sig"], st["a"], st["b"], st["lam"], st["nu"]
    s2, se2 = R.S2, st["s_eps2"]
    Soo = Sig; Sop = Sig @ a.T; Sot = Sig @ b.T
    Spp = a @ Sig @ a.T + s2 * np.eye(R.P)
    Stp = b @ Sig @ a.T + s2 * lam
    Stt = b @ Sig @ b.T + s2 * (lam @ lam.T) + nu @ nu.T + se2 * np.eye(R.D_T)
    return np.block([[Soo, Sot, Sop], [Sot.T, Stt, Stp], [Sop.T, Stp.T, Spp]])


def frozen_rule_scores(Sigma, alpha):
    """Frozen g1_pick algebra applied to an exact covariance."""
    d = Sigma.shape[0]
    sd = np.sqrt(np.diag(Sigma)); S = Sigma / np.outer(sd, sd)                # standardise
    S = (1 - alpha) * S + alpha * (np.trace(S) / d) * np.eye(d)                # shrink
    o, t = R.D_O, R.D_T; io, it, ip = slice(0, o), slice(o, o + t), slice(o + t, d)
    Soo = S[io, io] + RIDGE * np.eye(o)                                        # frozen ridge retained
    sol = np.linalg.solve(Soo, S[io, ip])
    C = S[it, ip] - S[it, io] @ sol
    v = np.diag(S[ip, ip] - S[ip, io] @ sol)
    return (C * C).sum(axis=0) / np.maximum(v, 1e-12)


def main():
    rows = []
    for c in R.conditions():
        st = R.structure(c); Sigma = population_cov(st)
        sG = st["G"]; sF = frozen_rule_scores(Sigma, 0.05); s0 = frozen_rule_scores(Sigma, 0.0)
        srt = lambda s: np.sort(s)[::-1]
        rows.append(dict(
            idx=c["idx"], mech=c["mech"], rho_pop=c["rho"], n=c["n"], A=c["A"], gamma_rel=c["gam"],
            jstar=int(st["jstar"]), gaussian_argmax=int(np.argmax(sG)),
            frozen_alpha05_argmax=int(np.argmax(sF)), alpha0_argmax=int(np.argmax(s0)),
            gaussian_top2_relgap=float((srt(sG)[0] - srt(sG)[1]) / srt(sG)[0]),
            frozen_alpha05_top2_relgap=float((srt(sF)[0] - srt(sF)[1]) / srt(sF)[0]),
            frozen_alpha05_top4=np.argsort(-sF)[:4].tolist(), alpha0_top4=np.argsort(-s0)[:4].tolist(),
            flip=bool(int(np.argmax(sF)) != int(np.argmax(sG))),
            frozen_limit_pick_is_jstar=bool(int(np.argmax(sF)) == int(st["jstar"])),
            # exact variance check: trace of Sigma_TT must equal D_T * c_marg (equal marginal variances by construction)
            marginal_var_check=float(np.max(np.abs(np.diag(Sigma)[R.D_O:R.D_O + R.D_T] / st["c_marg"] - 1.0))),
        ))
    flips = [r["idx"] for r in rows if r["flip"]]
    summary = dict(
        n_conditions=len(rows), n_flips=len(flips), flipped_conditions=flips,
        flips_by_gamma_rel={g: sum(1 for r in rows if r["flip"] and r["gamma_rel"] == g) for g in (0.02, 0.20)},
        m2_flips_land_on_jstar=all(r["frozen_limit_pick_is_jstar"] for r in rows if r["flip"] and r["mech"] == "M2"),
        alpha0_equals_gaussian_everywhere=all(r["alpha0_argmax"] == r["gaussian_argmax"] for r in rows),
        max_marginal_var_check=max(r["marginal_var_check"] for r in rows),
        positive_control_43_flip=rows[43]["flip"], primary_31_flip=rows[31]["flip"],
        primary_31_frozen_top2_relgap=rows[31]["frozen_alpha05_top2_relgap"],
        condition_3_frozen_top4=rows[3]["frozen_alpha05_top4"], ridge=RIDGE,
    )
    out = dict(summary=summary, rows=rows)
    (HERE / "nova_d_out").mkdir(exist_ok=True)
    (HERE / "nova_d_out" / "exact_ordering.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"\n{'idx':>3} {'mech':4} {'rho':>4} {'n':>5} {'A':>4} {'Grel':>5} | {'G':>2} {'frz':>3} {'a0':>3} {'j*':>2} | {'gapG':>7} {'gapF':>7} | flip")
    for r in rows:
        print(f"{r['idx']:>3} {r['mech']:4} {r['rho_pop']:>4} {r['n']:>5} {r['A']:>4} {r['gamma_rel']:>5} | {r['gaussian_argmax']:>2} {r['frozen_alpha05_argmax']:>3} {r['alpha0_argmax']:>3} {r['jstar']:>2} | {r['gaussian_top2_relgap']:7.4f} {r['frozen_alpha05_top2_relgap']:7.4f} | {'FLIP' if r['flip'] else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
