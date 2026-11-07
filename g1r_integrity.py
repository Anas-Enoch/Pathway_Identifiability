#!/usr/bin/env python3
"""G1-R replicate-level output integrity check v0.2.0 (read-only; required before `g1r_run.py analyze`).

Validates TYPES AND DOMAINS before any arithmetic: integer ids and picks (booleans excluded), picks in [0, 20),
regrets finite and non-negative and equal to the exact regret recomputed from the frozen structure, A1 discrepancies
finite and in [0, GATE_TOL], lambda_c a member of the frozen search set, delta_lw finite in [0, 1], kappa_hat finite
and >= 1, flags actual booleans. Rejects unparsable or malformed records (including non-dict `picks`), duplicate or
missing (condition, replicate) keys, and any freeze mismatch.

Report handling: any existing report is invalidated (removed) before checking; every exception yields a FAIL report;
the completed report is written atomically (temp file + os.replace). The report is bound to the replicate log, the
freeze record and this checker by SHA-256, which `analyze` re-verifies. Exit 0 only on PASS.
"""
from __future__ import annotations
import hashlib, json, math, os, sys, traceback
from collections import Counter
from pathlib import Path

CHECKER_VERSION = "0.2.0"
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import g1r_run as R   # noqa: E402  (import only)

REQUIRED = ("idx", "rep", "picks", "regret", "A1_C", "A1_v", "A1_s", "lambda_c", "endpoint_hit", "delta_lw", "complete_shrinkage", "A4_exact_top_tie", "kappa_hat")
BOOL_FIELDS = ("endpoint_hit", "complete_shrinkage", "A4_exact_top_tie")


def is_int(x):      return isinstance(x, int) and not isinstance(x, bool)
def is_num(x):      return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)
def is_bool(x):     return isinstance(x, bool)


def validate(d, conds):
    """Return a list of field-level problems for one parsed record (empty list = valid). No arithmetic on invalid fields."""
    p = []
    if not isinstance(d, dict): return ["record is not an object"]
    for k in REQUIRED:
        if k not in d: p.append(f"missing {k}")
    if p: return p
    if not (is_int(d["idx"]) and d["idx"] in conds): p.append("idx not an integer condition index")
    if not (is_int(d["rep"]) and 0 <= d["rep"] < R.R_REPS): p.append("rep not an integer in range")
    for name in ("picks", "regret"):
        if not (isinstance(d[name], dict) and set(d[name]) == set(R.ARMS)): p.append(f"{name} not an object keyed by the arms")
    if p: return p
    for arm, j in d["picks"].items():
        if not (is_int(j) and 0 <= j < R.P): p.append(f"pick {arm} invalid")
    for arm, r in d["regret"].items():
        if not (is_num(r) and r >= 0.0): p.append(f"regret {arm} not a finite non-negative number")
    for k in ("A1_C", "A1_v", "A1_s"):
        if not (is_num(d[k]) and 0.0 <= d[k] <= R.GATE_TOL): p.append(f"{k} not finite in [0, GATE_TOL]")
    if not (is_num(d["lambda_c"]) and float(d["lambda_c"]) in R.LAMBDA_SET): p.append("lambda_c not in the frozen search set")
    if not (is_num(d["delta_lw"]) and 0.0 <= d["delta_lw"] <= 1.0): p.append("delta_lw not finite in [0, 1]")
    if not (is_num(d["kappa_hat"]) and d["kappa_hat"] >= 1.0): p.append("kappa_hat not finite and >= 1")
    for k in BOOL_FIELDS:
        if not is_bool(d[k]): p.append(f"{k} not a boolean")
    return p


def write_report(rep):
    R.OUT.mkdir(exist_ok=True); tmp = R.INTEGRITY.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(rep, indent=1, default=str)); os.replace(tmp, R.INTEGRITY)


def main():
    # invalidate any previous report before doing anything else
    try:
        if R.INTEGRITY.exists(): R.INTEGRITY.unlink()
    except OSError as e:
        print(f"cannot remove previous report: {e}"); return 2
    rep = {"ok": False, "checker_version": CHECKER_VERSION, "problems": []}
    try:
        rep["checker_sha256"] = R.sha_file(Path(__file__))
        try:
            f = R.check_frozen(); rep["freeze_sha256"] = R.sha_file(R.FREEZE)
            if f.get("g1r_integrity.py") != rep["checker_sha256"]: rep["problems"].append("checker hash differs from the frozen one")
        except SystemExit as e:
            rep["problems"].append(f"freeze check: {e}"); f = None
        if not R.RUNLOG.exists():
            rep["problems"].append("no replicate log"); write_report(rep); print("INTEGRITY FAIL"); return 1
        rep["replicates_sha256"] = R.sha_file(R.RUNLOG)
        conds = {c["idx"]: c for c in R.conditions()}; want = {(i, r) for i in conds for r in range(R.R_REPS)}
        lines = R.RUNLOG.read_text().splitlines()
        recs, unparsable, invalid = [], [], []
        for ln, line in enumerate(lines, 1):
            if not line.strip(): continue
            try:
                d = json.loads(line)
            except json.JSONDecodeError as e:
                unparsable.append({"line": ln, "error": str(e)[:80]}); continue
            probs = validate(d, conds)
            if probs: invalid.append({"line": ln, "problems": probs[:6]}); continue
            recs.append((ln, d))
        rep.update(lines=len(lines), valid_records=len(recs), unparsable_count=len(unparsable), unparsable=unparsable[:50],
                   invalid_count=len(invalid), invalid=invalid[:50])
        if unparsable: rep["problems"].append(f"{len(unparsable)} unparsable line(s)")
        if invalid: rep["problems"].append(f"{len(invalid)} invalid record(s) (type/domain)")
        keys = Counter((d["idx"], d["rep"]) for _, d in recs); dups = {k: n for k, n in keys.items() if n > 1}
        if dups: rep["duplicate_count"] = len(dups); rep["problems"].append(f"{len(dups)} duplicated key(s)")
        have = set(keys); missing = sorted(want - have); extra = sorted(have - want)
        incomplete = sorted(i for i in conds if sum(1 for (j, _) in have if j == i) != R.R_REPS)
        rep.update(distinct_keys=len(have), missing=len(missing), extra=len(extra), conditions_incomplete=incomplete)
        if missing or extra or incomplete: rep["problems"].append(f"completeness: {len(missing)} missing, {len(extra)} extra, {len(incomplete)} incomplete condition(s)")
        st_cache, mism = {}, []
        for ln, d in recs:
            st = st_cache.setdefault(d["idx"], R.structure(conds[d["idx"]]))
            for arm, j in d["picks"].items():
                if abs(d["regret"][arm] - R.regret(st, j)) > 1e-12: mism.append({"line": ln, "arm": arm})
        rep.update(regret_mismatch_count=len(mism), regret_mismatches=mism[:50])
        if mism: rep["problems"].append(f"{len(mism)} regret(s) inconsistent with pick")
        rep["ok"] = not rep["problems"]
    except Exception:  # noqa: BLE001  — any failure is a FAIL report, never a missing one
        rep["ok"] = False; rep["problems"].append("checker exception: " + traceback.format_exc().strip().splitlines()[-1])
    write_report(rep)
    print(json.dumps({k: v for k, v in rep.items() if k not in ("unparsable", "invalid", "regret_mismatches")}, indent=1, default=str))
    print("INTEGRITY", "PASS" if rep["ok"] else "FAIL"); return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
