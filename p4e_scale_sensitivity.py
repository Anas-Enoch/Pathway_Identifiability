#!/usr/bin/env python3
"""Paper IV, 186-study primary: a post hoc sensitivity analysis for degenerate target scales.

POST HOC. This is not confirmatory and cannot replace the frozen primary.

Input: record_results.csv from paper4_predictor_benchmark. The owner verified all 393 files
against its SHA256SUMS.txt (393 OK); this copy came through a terminal paste. To confirm the
copy is intact, the script first reproduces the archived primary exactly (mean G1 - RIG of
0.019297341311863, with 51 positive, 38 negative and 97 zero), and stops if it cannot.

The exclusion rules are declared here, before any contrast is computed. Each uses only the RMSE
level of the archived G1 and RIG policies, never the contrast:
  A (numerical resolution): max(RMSE_G1, RMSE_RIG) > 1e9. Above this level a double-precision
    pooled RMSE cannot resolve selector differences of order 1e-7 or smaller.
  B (scale plausibility):   max(RMSE_G1, RMSE_RIG) > 10. A design-standardized error of more
    than ten design SDs points to a near-degenerate design-half target scale.
Contrast: d = RMSE_G1 - RMSE_RIG, so positive values favour RIG (manuscript sign).
Interval: a percentile bootstrap over records, B = 10000, seed 20261017. This approximates the
frozen equal-study bootstrap; it does not reproduce it exactly.

Usage: python3 -I p4e_scale_sensitivity.py PATH/TO/record_results.csv
"""
import csv
import hashlib
import json
import random
import statistics as st
import sys

VERSION = "1.0.0"
TARGET_MEAN = 0.019297341311863003
RULES = {"A_resolution_gt_1e9": 1e9, "B_plausibility_gt_10": 10.0}


def summarize(d, B=10000, seed=20261017):
    rng = random.Random(seed)
    n = len(d)
    boots = sorted(st.mean(rng.choices(d, k=n)) for _ in range(B))
    return {"n_records": n, "mean": st.mean(d), "median": st.median(d),
            "positive": sum(x > 0 for x in d), "negative": sum(x < 0 for x in d),
            "zero": sum(x == 0 for x in d),
            "bootstrap_ci95": [boots[int(0.025 * B)], boots[int(0.975 * B) - 1]]}


def main(path):
    raw = open(path, "rb").read()
    rows = list(csv.DictReader(raw.decode().splitlines()))
    rec = {}
    for r in rows:
        if r["policy"] in ("archive_G1", "archive_RIG"):
            rec.setdefault(r["record"], {})[r["policy"]] = float(r["rmse"])
    assert len(rec) == 186 and all(len(v) == 2 for v in rec.values())
    d_all = {k: v["archive_G1"] - v["archive_RIG"] for k, v in rec.items()}
    full = summarize(list(d_all.values()))
    if not (abs(full["mean"] - TARGET_MEAN) < 1e-12 and
            (full["positive"], full["negative"], full["zero"]) == (51, 38, 97)):
        sys.exit(f"reproduction failed: {full}")
    level = {k: max(v.values()) for k, v in rec.items()}
    lv = sorted(level.values())
    out = {"version": VERSION, "label": "POST_HOC_SENSITIVITY_NOT_CONFIRMATORY",
           "input_sha256": hashlib.sha256(raw).hexdigest(),
           "reproduced_primary": full,
           "rmse_level_distribution": {
               "min": lv[0], "q25": lv[46], "median": st.median(lv), "q75": lv[139],
               "q90": lv[167], "max": lv[-1],
               "n_gt_2": sum(x > 2 for x in lv), "n_gt_10": sum(x > 10 for x in lv),
               "n_gt_1e3": sum(x > 1e3 for x in lv), "n_gt_1e9": sum(x > 1e9 for x in lv)},
           "rules": {}}
    for name, thr in RULES.items():
        excl = sorted(k for k in rec if level[k] > thr)
        kept = [d_all[k] for k in rec if level[k] <= thr]
        out["rules"][name] = {"threshold": thr, "n_excluded": len(excl),
                              "excluded_nonzero_contrasts": {k: d_all[k] for k in excl if d_all[k] != 0},
                              "excluded_records": excl if len(excl) <= 40 else excl[:40] + ["..."],
                              "summary": summarize(kept)}
    return out


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    print(json.dumps(main(sys.argv[1]), indent=2))
