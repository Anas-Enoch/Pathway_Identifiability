#!/usr/bin/env python3
"""NOVA-D output integrity check (read-only; run BEFORE `nova_d_run.py analyze`).

Does not modify the frozen implementation or any output. Verifies, per the owner's requirement of 5 Nov 2026:
  1. every JSONL record parses;
  2. every (condition, replicate) key is unique — and if a duplicate exists, whether the copies are byte-identical;
  3. exactly 72 x 100 = 7,200 distinct keys; no condition missing any replicate;
  4. every record is consistent with the frozen configuration: condition index in range, 5 folds, all arm and
     policy sums present and finite, jstar / jG equal to the frozen structure's, gamma list and picks of length 5;
  5. the frozen files still match FREEZE.json.
Exit 0 only if every check passes. Prints a JSON report either way; writes nova_d_out/integrity.json.
"""
from __future__ import annotations
import json, math, sys
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nova_d_run as R   # noqa: E402  (import only; nothing is executed or written)


def main():
    rep = {"ok": True, "problems": []}
    bad = rep["problems"].append
    # 5. freeze still valid
    f = json.loads(R.FREEZE.read_text())
    for k, p in (("protocol_sha256", HERE / R.PROTOCOL), ("nova_d_run.py", HERE / "nova_d_run.py"), ("nova_d_calib.py", HERE / "nova_d_calib.py")):
        if R.sha_file(p) != f[k]:
            bad(f"freeze mismatch: {p.name}")
    conds = {c["idx"]: c for c in R.conditions()}
    want = {(i, r) for i in conds for r in range(R.R_REPS)}
    # 1. parse
    lines = R.RUNLOG.read_text().splitlines()
    recs, unparsable = [], []
    for ln, line in enumerate(lines, 1):
        if not line.strip():
            continue
        try:
            recs.append((ln, json.loads(line)))
        except json.JSONDecodeError as e:
            unparsable.append({"line": ln, "error": str(e)[:80]})
    rep["lines"] = len(lines); rep["parsed"] = len(recs); rep["unparsable"] = unparsable
    if unparsable:
        bad(f"{len(unparsable)} unparsable line(s)")
    # 2. uniqueness
    keys = Counter((d["idx"], d["rep"]) for _, d in recs)
    dups = {k: n for k, n in keys.items() if n > 1}
    rep["duplicate_keys"] = len(dups)
    if dups:
        ident = {}
        for k in dups:
            copies = [line for (ln, d), line in zip(recs, [l for l in lines if l.strip()]) if (d["idx"], d["rep"]) == k]
            ident[str(k)] = len(set(copies)) == 1
        rep["duplicates_byte_identical"] = ident
        bad(f"{len(dups)} duplicated key(s); identical copies: {sum(ident.values())}/{len(ident)}")
    # 3. completeness
    have = set(keys)
    missing = sorted(want - have); extra = sorted(have - want)
    rep["distinct_keys"] = len(have); rep["missing"] = len(missing); rep["extra"] = len(extra)
    per_cond = Counter(i for i, _ in have)
    rep["conditions_incomplete"] = sorted(i for i in conds if per_cond.get(i, 0) != R.R_REPS)
    if missing or extra or rep["conditions_incomplete"]:
        bad(f"completeness: {len(missing)} missing, {len(extra)} extra, {len(rep['conditions_incomplete'])} condition(s) incomplete")
    # 4. configuration consistency
    st_cache, incons = {}, []
    for ln, d in recs:
        i = d.get("idx")
        if i not in conds:
            incons.append((ln, "idx out of range")); continue
        if i not in st_cache:
            st_cache[i] = R.structure(conds[i])
        st = st_cache[i]
        probs = []
        if d.get("folds") != 5: probs.append("folds != 5")
        if d.get("jstar") != st["jstar"] or d.get("jG") != st["jG"]: probs.append("jstar/jG mismatch")
        for a in ("BASE", "G1", "NOVA-V", "NOVA-U", "NOVA-S"):
            v = d.get("sums", {}).get(a)
            if v is None or not math.isfinite(v) or v <= 0: probs.append(f"sum {a}")
        for g in R.GAMMAS:
            v = d.get("policy", {}).get(str(g))
            if v is None or not math.isfinite(v) or v <= 0: probs.append(f"policy {g}")
        for key in ("picks", "t", "gamma_obs", "oracle_regret"):
            if len(d.get(key, [])) != 5: probs.append(f"{key} length")
        if probs:
            incons.append((ln, i, d.get("rep"), probs))
    rep["inconsistent_records"] = incons[:50]; rep["inconsistent_count"] = len(incons)
    if incons:
        bad(f"{len(incons)} record(s) inconsistent with the frozen configuration")
    rep["ok"] = not rep["problems"]
    (R.OUT / "integrity.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps({k: v for k, v in rep.items() if k != "inconsistent_records"}, indent=1))
    print("INTEGRITY", "PASS" if rep["ok"] else "FAIL")
    return 0 if rep["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
