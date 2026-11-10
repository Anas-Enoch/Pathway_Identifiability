#!/usr/bin/env python3
"""
nhanes_fa_stability_independent_check.py  v1.0   Independent recomputation of the NHANES-FA stability v0.4 screen.

Reads ONLY the saved aggregate outputs (prediction table, screen report, development fit aggregate). It imports
nothing from nhanes_fa_stability.py, opens no participant-level file and no NHANES file, and reruns no simulation.
Every prediction set is recomputed from the saved integer histograms with integer arithmetic:

    L = min{c : 40 * cum(c) >= B},      U = min{c : 40 * cum(c) >= 39 * B}      (the 2.5 % / 97.5 % rule)

which removes the floating-point comparison used by quantile_set() in runner v0.1.0 (cumsum of proportions >= 0.975
can evaluate to 0.9749999999999999 when the cumulative count is exactly 0.975 * B).

Sections 1-4 verify the saved record. Section 5 is DESCRIPTIVE and POST HOC: it was not declared in protocol v0.4
and decides nothing.

  python3 -I nhanes_fa_stability_independent_check.py PREDICTION.json SCREEN.json DEV_FIT.json [--out result.json]
"""
import argparse, hashlib, itertools, json, sys

U_CAP = 5
N_LO, N_HI, N_DEV, B_FROZEN = 1200, 2400, 1725, 2000
ATTAINABLE_BLOCK = (0, 4, 6, 7, 8, 9, 10)


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def cum(h):
    out, s = [], 0
    for k in h:
        s += int(k); out.append(s)
    return out


def q_int(h):
    c = cum(h); b = c[-1]
    lo = next(i for i, x in enumerate(c) if 40 * x >= b)
    hi = next(i for i, x in enumerate(c) if 40 * x >= 39 * b)
    return lo, hi, b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prediction"); ap.add_argument("screen"); ap.add_argument("fit"); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    pred = json.load(open(a.prediction)); sc = json.load(open(a.screen)); fit = json.load(open(a.fit))
    att_c = {sum(x) for x in itertools.product(ATTAINABLE_BLOCK, repeat=5)}
    res = {"check": "nhanes_fa_stability_independent_check.py v1.0",
           "inputs_sha256": {"prediction": sha(a.prediction), "screen": sha(a.screen), "fit": sha(a.fit)}}
    fails = []

    # 1. structure
    st = {"complete_flag": pred["complete"], "B": pred["B"], "verdict_range": pred["verdict_range"],
          "screen_cites_prediction": sc["prediction_sha256"] == res["inputs_sha256"]["prediction"],
          "prediction_cites_fit": pred["fit_sha256"] == res["inputs_sha256"]["fit"]}
    for m in ("MG", "ME"):
        t = pred["table"][m]; sizes = sorted(int(k) for k in t)
        st[m] = {"n_sizes": len(sizes), "contiguous_1200_2400": sizes == list(range(N_LO, N_HI + 1)),
                 "hist_sum_ne_B": sum(1 for n in sizes if sum(t[str(n)]["hist_C"]) != B_FROZEN),
                 "mass_on_unattainable_C": sum(k for n in sizes for c, k in enumerate(t[str(n)]["hist_C"]) if c not in att_c),
                 "cycles": sum(sum(t[str(n)]["hist_C"]) for n in sizes)}
        if not (st[m]["contiguous_1200_2400"] and st[m]["hist_sum_ne_B"] == 0 and st[m]["mass_on_unattainable_C"] == 0):
            fails.append(f"structure {m}")
    if not (st["screen_cites_prediction"] and st["prediction_cites_fit"] and st["complete_flag"]):
        fails.append("hash chain / completeness")
    res["1_structure"] = st

    # 2. integer-rule prediction sets vs saved
    q = {}
    for m in ("MG", "ME"):
        t = pred["table"][m]; disc, udist, c4, c0 = [], {}, {}, {}
        for n in range(N_LO, N_HI + 1):
            r = t[str(n)]; lo, hi, b = q_int(r["hist_C"]); c = cum(r["hist_C"])
            udist[hi] = udist.get(hi, 0) + 1; c4[n] = c[4]; c0[n] = c[0]
            if (lo, hi) != (r["L"], r["U"]):
                disc.append({"n": n, "saved": [r["L"], r["U"]], "integer_rule": [lo, hi],
                             "cum_below_correct_U": c[hi - 1] if hi else 0, "cum_at_correct_U": c[hi]})
        nmax = max(c4, key=c4.get); nmin = min(c4, key=c4.get)
        q[m] = {"U_distribution_integer_rule": {str(k): v for k, v in sorted(udist.items())}, "discrepancies": disc,
                "count_C_le_4": {"needed_for_U_le_4": 39 * B_FROZEN // 40, "max": c4[nmax], "at_n_max": nmax, "min": c4[nmin], "at_n_min": nmin,
                                 "mean": round(sum(c4.values()) / len(c4), 1), "at_n_dev": c4[N_DEV]},
                "count_C_eq_0": {"max": max(c0.values()), "min": min(c0.values()), "at_n_dev": c0[N_DEV]},
                "by_size_band": {f"{lo}-{hi}": {"mean_P_C0": round(sum(c0[n] for n in range(lo, hi + 1)) / ((hi - lo + 1) * B_FROZEN), 4),
                                                "mean_P_C_le_4": round(sum(c4[n] for n in range(lo, hi + 1)) / ((hi - lo + 1) * B_FROZEN), 4)}
                                 for lo, hi in ((1200, 1499), (1500, 1799), (1800, 2099), (2100, 2400))}}
    res["2_prediction_sets"] = q

    # 3. screen decisions recomputed
    mg = pred["table"]["MG"]; d = mg[str(N_DEV)]; lo, hi, _ = q_int(d["hist_C"])
    assessable = [n for n in range(N_LO, N_HI + 1) if q_int(mg[str(n)]["hist_C"])[1] <= U_CAP]
    modal = [max(range(9), key=lambda j, b=b: (d["pick_counts"][b][j], -j)) for b in range(5)]
    s3 = {"a_dev_C_in_set": {"L": lo, "U": hi, "dev_C": fit["dev_C"], "pass": lo <= fit["dev_C"] <= hi},
          "b_assessable_sizes": len(assessable), "c_modal_picks": {"predicted": modal, "development": fit["full_fit_picks"], "pass": modal == fit["full_fit_picks"]},
          "d_gaps_positive": {"gaps": fit["full_fit_signed_gap"], "pass": all(g > 0 for g in fit["full_fit_signed_gap"])}}
    s3["freeze_permitted"] = bool(s3["a_dev_C_in_set"]["pass"] and assessable and s3["c_modal_picks"]["pass"] and s3["d_gaps_positive"]["pass"])
    s3["agrees_with_saved_screen"] = (s3["freeze_permitted"] == sc["freeze_permitted"] and len(assessable) == sc["b_sharpness"]["n_assessable"]
                                      and s3["a_dev_C_in_set"]["pass"] == sc["a_dev_C_in_set"]["pass"] and (lo, hi) == (sc["a_dev_C_in_set"]["L"], sc["a_dev_C_in_set"]["U"])
                                      and s3["c_modal_picks"]["pass"] == sc["c_modal_picks"]["pass"] and s3["d_gaps_positive"]["pass"] == sc["d_full_fit_gap_positive"]["pass"])
    if not s3["agrees_with_saved_screen"]: fails.append("screen decisions")
    res["3_screen"] = s3

    # 4. declared per-block quantities at the development size (protocol s2 items 2-3)
    blk = {}
    for m in ("MG", "ME"):
        r = pred["table"][m][str(N_DEV)]; h = r["hist_C"]
        blk[m] = {"hist_C_nonzero": {str(c): k for c, k in enumerate(h) if k}, "P_C0": h[0] / sum(h),
                  "modal_prob_by_block": r["modal_prob"], "per_fold_reversal_prob_by_block": r["ds_reversal_prob"]}
    blk["fixed_pairs"] = fit["pick_names"]
    blk["full_fit_top_two_scores"] = [sorted(s, reverse=True)[:2] for s in fit["full_fit_scores"]]
    blk["full_fit_relative_gap"] = [round((a_ - b_) / a_, 4) for a_, b_ in blk["full_fit_top_two_scores"]]
    res["4_blocks_at_n_dev"] = blk

    # 5. DESCRIPTIVE, POST HOC (not in protocol v0.4; decides nothing): how often would the secondary model M_E have
    #    produced a count outside the primary model's set? Shown for the realised sets (U = 6) and for a hypothetical U = 4.
    me = pred["table"]["ME"]; out6, out4 = {}, {}
    for n in range(N_LO, N_HI + 1):
        c = cum(me[str(n)]["hist_C"]); u = q_int(mg[str(n)]["hist_C"])[1]
        out6[n] = (c[-1] - c[u]) / c[-1]; out4[n] = (c[-1] - c[4]) / c[-1]
    res["5_descriptive_post_hoc"] = {
        "label": "post hoc, descriptive; plug-in Monte Carlo proportions; not a declared analysis",
        "P_ME_C_above_MG_U": {"at_n_dev": round(out6[N_DEV], 4), "min": round(min(out6.values()), 4), "max": round(max(out6.values()), 4),
                              "mean": round(sum(out6.values()) / len(out6), 4)},
        "P_ME_C_above_4_hypothetical": {"at_n_dev": round(out4[N_DEV], 4), "min": round(min(out4.values()), 4), "max": round(max(out4.values()), 4),
                                        "mean": round(sum(out4.values()) / len(out4), 4)}}

    res["verdict"] = "PASS: saved screen reproduced" if not fails else f"FAIL: {fails}"
    txt = json.dumps(res, indent=1)
    if a.out: open(a.out, "w").write(txt + "\n")
    print(txt)
    return 0 if not fails else 1


if __name__ == "__main__":
    sys.exit(main())
