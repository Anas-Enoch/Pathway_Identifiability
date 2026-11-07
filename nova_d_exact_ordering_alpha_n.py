#!/usr/bin/env python3
"""Exact population ordering of the frozen G1 rule at the shrinkage level IN FORCE at each condition's n.

Companion to nova_d_exact_ordering.py (which evaluates the rule at its alpha floor 0.05, the n -> infinity limit).
The frozen rule sets alpha_n = min(1, max(0.05, d/(d + n_train))) with d = 60 and n_train the outer-fold training
size int(0.8 n). In the NOVA-D grid that is 0.3704 (n = 128), 0.1279 (n = 512) and 0.0683 (n = 1024); the floor binds
only for n > 1140. The population object that finite-sample picks at a given n scatter around is therefore the rule
applied to the exact covariance at alpha_n, not at the floor. Read-only; no sampling; nothing frozen is touched.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nova_d_run as R                      # noqa: E402
import nova_d_exact_ordering as X           # noqa: E402  (population_cov, frozen_rule_scores; ridge retained)

D = R.D_O + R.D_T + R.P


def alpha_at(n: int) -> float:
    n_train = int(0.8 * n)
    return min(1.0, max(0.05, D / (D + n_train)))


def main():
    rows = []
    for c in R.conditions():
        st = R.structure(c); Sig = X.population_cov(st); a_n = alpha_at(c["n"])
        sG, sN, sF = st["G"], X.frozen_rule_scores(Sig, a_n), X.frozen_rule_scores(Sig, 0.05)
        srt = lambda s: np.sort(s)[::-1]
        rows.append(dict(
            idx=c["idx"], mech=c["mech"], rho_pop=c["rho"], n=c["n"], A=c["A"], gamma_rel=c["gam"], alpha_n=a_n,
            jstar=int(st["jstar"]), gaussian_argmax=int(np.argmax(sG)),
            argmax_at_alpha_n=int(np.argmax(sN)), argmax_at_floor=int(np.argmax(sF)),
            top5_at_alpha_n=np.argsort(-sN)[:5].tolist(), scores_top5_at_alpha_n=np.round(srt(sN)[:5], 4).tolist(),
            top2_relgap_at_alpha_n=float((srt(sN)[0] - srt(sN)[1]) / srt(sN)[0]),
            flip_at_alpha_n=bool(int(np.argmax(sN)) != int(np.argmax(sG))),
            flip_at_floor=bool(int(np.argmax(sF)) != int(np.argmax(sG))),
            alpha_n_pick_is_jstar=bool(int(np.argmax(sN)) == int(st["jstar"])),
        ))
    fl_n = [r["idx"] for r in rows if r["flip_at_alpha_n"]]; fl_f = [r["idx"] for r in rows if r["flip_at_floor"]]
    summary = dict(
        alpha_by_n={n: alpha_at(n) for n in R.N_GRID},
        n_flips_at_alpha_n=len(fl_n), flipped_at_alpha_n=fl_n,
        n_flips_at_floor=len(fl_f), flipped_at_floor=fl_f,
        flips_at_alpha_n_by_gamma_rel={g: sum(1 for r in rows if r["flip_at_alpha_n"] and r["gamma_rel"] == g) for g in (0.02, 0.20)},
        flips_at_alpha_n_by_n={n: sum(1 for r in rows if r["flip_at_alpha_n"] and r["n"] == n) for n in R.N_GRID},
        m2_flips_at_alpha_n_land_on_jstar=all(r["alpha_n_pick_is_jstar"] for r in rows if r["flip_at_alpha_n"] and r["mech"] == "M2"),
        primary_31=dict(alpha_n=rows[31]["alpha_n"], argmax=rows[31]["argmax_at_alpha_n"], is_jstar=rows[31]["alpha_n_pick_is_jstar"],
                        top2_relgap=rows[31]["top2_relgap_at_alpha_n"]),
        positive_control_43=dict(alpha_n=rows[43]["alpha_n"], argmax=rows[43]["argmax_at_alpha_n"], is_jstar=rows[43]["alpha_n_pick_is_jstar"],
                                 top2_relgap=rows[43]["top2_relgap_at_alpha_n"]),
        condition_3=dict(alpha_n=rows[3]["alpha_n"], top5=rows[3]["top5_at_alpha_n"], scores=rows[3]["scores_top5_at_alpha_n"]),
    )
    out = dict(summary=summary, rows=rows)
    (HERE / "nova_d_out").mkdir(exist_ok=True)
    (HERE / "nova_d_out" / "exact_ordering_alpha_n.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(summary, indent=1))
    print(f"\n{'idx':>3} {'mech':4} {'rho':>4} {'n':>5} {'A':>4} {'a_n':>6} | {'G':>2} {'a_n':>3} {'flr':>3} {'j*':>2} | {'gap@a_n':>8} | flip@a_n  top5@a_n")
    for r in rows:
        print(f"{r['idx']:>3} {r['mech']:4} {r['rho_pop']:>4} {r['n']:>5} {r['A']:>4} {r['alpha_n']:6.3f} | {r['gaussian_argmax']:>2} {r['argmax_at_alpha_n']:>3} {r['argmax_at_floor']:>3} {r['jstar']:>2} | {r['top2_relgap_at_alpha_n']:8.4f} | {'FLIP' if r['flip_at_alpha_n'] else '    '}  {r['top5_at_alpha_n']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
