#!/usr/bin/env python3
"""
nhanes_fa_stability_margin_check.py  v1.1   Two POST HOC, DESCRIPTIVE checks on the NHANES-FA stability v0.4 screen.
Neither was declared in protocol v0.4; neither changes the screen decision, which uses the saved counts.

(a) Selection-adjusted sharpness margin (needs only the saved prediction table).
    The best single-size count of Gaussian cycles with C <= 4 is the maximum of 1,201 independently seeded Monte Carlo
    counts, so it overstates the fitted model's probability at that size. A binomial GLM (logit link, quadratic in n)
    is fitted to the 1,201 counts; the Pearson dispersion checks that the scatter about the trend is binomial Monte
    Carlo noise. Reported: trend value and shortfall from the 1,950 needed, at chosen sizes; the normal-approximation
    chance that a fresh 2,000-cycle run at that size reaches 1,950; the expected number of sizes that would reach 1,950
    by Monte Carlo chance alone under the trend.

(b) Reproduction and per-block attribution at n = 1,725 (needs the runner and the frozen modules; --runner-dir).
    Re-runs the 2,000 cycles of a model from its frozen seed, checks the histogram and the fixed-pair reversal counts
    against the saved table, and tabulates the per-cycle vector (c_1..c_5). The Gaussian model needs only the aggregate
    fit; the empirical model also needs the participant-level rows (--rows dev_rows.npy), which stay local.
    v1.1 adds, per block, the number of cycles in which exactly k of the five folds picked the development leader
    (k = 0..5), the cycles that are unanimous on a different candidate (they contribute 0 to C), and the share of
    cycles in which all 25 picks equal the development picks. All v1.0 fields are unchanged.

  python3 -I nhanes_fa_stability_margin_check.py PREDICTION.json DEV_FIT.json [--runner-dir DIR] [--rows dev_rows.npy] [--out result.json]
"""
import argparse, collections, hashlib, importlib.util, json, math, os, sys
import numpy as np

B, N_LO, N_HI, N_DEV, NEED = 2000, 1200, 2400, 1725, 1950


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def upper_tail(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def part_a(pred):
    ns = np.arange(N_LO, N_HI + 1)
    c4 = np.array([sum(pred["table"]["MG"][str(n)]["hist_C"][:5]) for n in ns], float)
    x = (ns - 1800) / 600.0
    X = np.column_stack([np.ones_like(x), x, x * x])
    beta = np.array([math.log(c4.mean() / (B - c4.mean())), 0.0, 0.0])
    for _ in range(100):
        eta = X @ beta; p = 1 / (1 + np.exp(-eta)); w = B * p * (1 - p)
        new = np.linalg.solve(X.T @ (w[:, None] * X), X.T @ (w * (eta + (c4 - B * p) / w)))
        done = np.max(np.abs(new - beta)) < 1e-12; beta = new
        if done: break
    p = 1 / (1 + np.exp(-(X @ beta))); se = np.sqrt(B * p * (1 - p)); z = (c4 - B * p) / se
    i = int(np.argmax(c4))
    pass_p = np.array([upper_tail((NEED - 0.5 - B * a) / s) for a, s in zip(p, se)])
    at = {}
    for n in (N_DEV, 2000, int(ns[i]), N_HI):
        k = n - N_LO
        at[str(n)] = {"observed_count": int(c4[k]), "trend_count": round(float(B * p[k]), 1), "trend_P_C_le_4": round(float(p[k]), 4),
                      "binomial_se_cycles": round(float(se[k]), 2), "observed_minus_trend_in_se": round(float(z[k]), 2),
                      "shortfall_observed": int(NEED - c4[k]), "shortfall_trend": round(float(NEED - B * p[k]), 1),
                      "shortfall_trend_in_se": round(float((NEED - B * p[k]) / se[k]), 2),
                      "chance_fresh_run_reaches_1950": float(f"{pass_p[k]:.1e}")}
    return {"label": "post hoc, descriptive; fitted-model Monte Carlo quantities; decides nothing",
            "model": "binomial GLM, logit link, quadratic in n, 1,201 Gaussian counts of C <= 4, independent seeds per size",
            "pearson_dispersion": round(float(np.sum(z ** 2) / (len(ns) - 3)), 3),
            "best_single_size": int(ns[i]), "largest_standardised_residual": round(float(z.max()), 2),
            "at_sizes": at,
            "expected_sizes_reaching_1950_by_chance": round(float(pass_p.sum()), 2),
            "chance_at_least_one_size_reaches_1950": round(float(1 - np.prod(1 - pass_p)), 2),
            "window_mean_count": {"2223-2323": round(float(c4[(ns >= 2223) & (ns <= 2323)].mean()), 1),
                                  "2350-2400": round(float(c4[(ns >= 2350) & (ns <= 2400)].mean()), 1)}}


def part_b(pred, fit, runner_dir, rows_path):
    spec = importlib.util.spec_from_file_location("nhanes_fa_stability_asrun", os.path.join(runner_dir, "nhanes_fa_stability.py"))
    st = importlib.util.module_from_spec(spec); spec.loader.exec_module(st)
    cols = st.block_cols(fit["allvars"]); pairs = fit["fixed_pairs"]; leaders = fit["full_fit_picks"]; out = {"runner_sha256": sha(os.path.join(runner_dir, "nhanes_fa_stability.py"))}
    models = [("MG", None)] + ([("ME", np.load(rows_path))] if rows_path else [])
    for model, rows in models:
        rng = np.random.default_rng(st.seed_for(model, N_DEV)); hist = [0] * 51; pat = collections.Counter(); rev = [0] * 5
        k_hist = [[0] * 6 for _ in range(5)]; unan_other = [0] * 5; all_dev = 0
        for _ in range(B):
            s = st.cycle_stats(st.draw(model, rng, N_DEV, fit, rows), cols, pairs)
            hist[s["C"]] += 1; pat[tuple(s["c"])] += 1; same = True
            for b in range(5):
                rev[b] += sum(1 for v in s["ds"][b] if v < 0)
                k = sum(1 for j in s["picks"][b] if j == leaders[b]); k_hist[b][k] += 1; same = same and k == 5
                if k == 0 and s["c"][b] == 0: unan_other[b] += 1
            all_dev += same
        saved = pred["table"][model][str(N_DEV)]; nz = B - hist[0]
        one_block = [sum(k for c, k in pat.items() if sum(1 for v in c if v) == 1 and c[b]) for b in range(5)]
        out[model] = {"histogram_identical_to_saved": hist == saved["hist_C"],
                      "reversal_counts_identical_to_saved": rev == [round(v * 5 * B) for v in saved["ds_reversal_prob"]],
                      "cycles_with_disagreement": nz, "cycles_with_disagreement_in_exactly_one_block_by_block": one_block,
                      "cycles_with_disagreement_in_two_or_more_blocks": nz - sum(one_block),
                      "patterns_c1_to_c5": [[list(c), k] for c, k in sorted(pat.items(), key=lambda kv: (-kv[1], kv[0]))],
                      "cycles_by_folds_picking_development_leader_k0_to_k5_by_block": k_hist,
                      "cycles_unanimous_on_a_different_candidate_by_block": unan_other,
                      "cycles_with_C_eq_0": hist[0], "cycles_with_all_25_picks_equal_to_development": int(all_dev)}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prediction"); ap.add_argument("fit"); ap.add_argument("--runner-dir", default=None)
    ap.add_argument("--rows", default=None); ap.add_argument("--out", default=None)
    a = ap.parse_args()
    pred = json.load(open(a.prediction)); fit = json.load(open(a.fit))
    res = {"check": "nhanes_fa_stability_margin_check.py v1.1",
           "inputs_sha256": {"prediction": sha(a.prediction), "fit": sha(a.fit)},
           "a_selection_adjusted_margin": part_a(pred)}
    if a.runner_dir:
        res["b_reproduction_and_attribution_at_n_dev"] = part_b(pred, fit, a.runner_dir, a.rows)
    txt = json.dumps(res, indent=1)
    if a.out: open(a.out, "w").write(txt + "\n")
    print(txt)
    print(f"[environment, not part of the result] numpy {np.__version__}, python {sys.version.split()[0]}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
