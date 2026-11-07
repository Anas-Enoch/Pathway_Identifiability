#!/usr/bin/env python3
"""G1-R replicate-level output integrity check (read-only; required before `g1r_run.py analyze`).

Rejects: malformed records; duplicate or missing (condition, replicate) keys; picks outside [0, 20); non-finite values;
regrets inconsistent with the recorded picks (recomputed exactly from the frozen structure); failed replicate gates
(A1 discrepancies above GATE_TOL on C, v or s); freeze hashes not matching the current files or the self-test report.
Writes g1r_out/integrity.json; exit 0 only on PASS.
"""
from __future__ import annotations
import json, math, sys
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import g1r_run as R   # noqa: E402  (import only)

REQUIRED = ("idx", "rep", "picks", "regret", "A1_C", "A1_v", "A1_s", "lambda_c", "endpoint_hit", "delta_lw", "complete_shrinkage", "A4_exact_top_tie", "kappa_hat")


def main():
    rep = {"ok": True, "problems": []}; bad = rep["problems"].append
    try:
        f = R.check_frozen()
    except SystemExit as e:
        bad(f"freeze check: {e}"); f = None
    conds = {c["idx"]: c for c in R.conditions()}; want = {(i, r) for i in conds for r in range(R.R_REPS)}
    lines = R.RUNLOG.read_text().splitlines() if R.RUNLOG.exists() else []
    recs, unparsable, malformed = [], [], []
    for ln, line in enumerate(lines, 1):
        if not line.strip(): continue
        try:
            d = json.loads(line)
        except json.JSONDecodeError as e:
            unparsable.append({"line": ln, "error": str(e)[:80]}); continue
        if not isinstance(d, dict) or any(k not in d for k in REQUIRED) or set(d["picks"]) != set(R.ARMS) or set(d["regret"]) != set(R.ARMS):
            malformed.append(ln); continue
        recs.append((ln, d))
    rep.update(lines=len(lines), parsed=len(recs), unparsable=unparsable[:50], unparsable_count=len(unparsable), malformed_lines=malformed[:50], malformed_count=len(malformed))
    if unparsable: bad(f"{len(unparsable)} unparsable line(s)")
    if malformed: bad(f"{len(malformed)} malformed record(s)")
    keys = Counter((d["idx"], d["rep"]) for _, d in recs); dups = {k: n for k, n in keys.items() if n > 1}
    if dups:
        text_by_key = {}
        for (ln, d), line in zip(recs, [l for l in lines if l.strip()]): text_by_key.setdefault((d["idx"], d["rep"]), set()).add(line)
        rep["duplicates_byte_identical"] = {str(k): len(text_by_key[k]) == 1 for k in dups}; bad(f"{len(dups)} duplicated key(s)")
    have = set(keys); missing = sorted(want - have); extra = sorted(have - want)
    rep.update(distinct_keys=len(have), missing=len(missing), extra=len(extra), conditions_incomplete=sorted(i for i in conds if sum(1 for (j, _) in have if j == i) != R.R_REPS))
    if missing or extra or rep["conditions_incomplete"]: bad(f"completeness: {len(missing)} missing, {len(extra)} extra, {len(rep['conditions_incomplete'])} incomplete condition(s)")
    st_cache, invalid, nonfinite, regret_mismatch, gate_fail = {}, [], [], [], []
    for ln, d in recs:
        i = d["idx"]
        if i not in conds: invalid.append((ln, "idx")); continue
        if i not in st_cache: st_cache[i] = R.structure(conds[i])
        st = st_cache[i]
        for arm, j in d["picks"].items():
            if not (isinstance(j, int) and 0 <= j < R.P): invalid.append((ln, arm)); continue
            if abs(d["regret"][arm] - R.regret(st, j)) > 1e-12: regret_mismatch.append((ln, arm))
        vals = [d["regret"][a] for a in R.ARMS] + [d["A1_C"], d["A1_v"], d["A1_s"], d["lambda_c"], d["delta_lw"], d["kappa_hat"]]
        if not all(isinstance(v, (int, float)) and math.isfinite(v) for v in vals): nonfinite.append(ln)
        if max(d["A1_C"], d["A1_v"], d["A1_s"]) > R.GATE_TOL: gate_fail.append(ln)
    rep.update(invalid_picks=invalid[:50], invalid_count=len(invalid), nonfinite_lines=nonfinite[:50], nonfinite_count=len(nonfinite),
               regret_mismatches=regret_mismatch[:50], regret_mismatch_count=len(regret_mismatch), gate_failures=gate_fail[:50], gate_failure_count=len(gate_fail))
    for cnt, label in ((invalid, "invalid pick(s)"), (nonfinite, "record(s) with non-finite values"), (regret_mismatch, "regret(s) inconsistent with pick"), (gate_fail, "record(s) with failed A1 gate")):
        if cnt: bad(f"{len(cnt)} {label}")
    rep["ok"] = not rep["problems"]
    R.OUT.mkdir(exist_ok=True); (R.OUT / "integrity.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: v for k, v in rep.items() if not k.endswith("s") or isinstance(v, int)}, indent=1)); print("INTEGRITY", "PASS" if rep["ok"] else "FAIL")
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
