"""Extraction validator for docrec T1 (schema + verbatim-quote + count-form checks only).
Usage: python3 -I validate.py <records.jsonl> [batch_k.json]
Prints VALID or the problems for each record. It never prints a grade."""
import sys, json, os, runpy
DOCREC = "/home/claude/nova/docrec/docrec.py"
HERE = os.path.dirname(os.path.abspath(__file__))
STORE = "/mnt/user-data/uploads/github.pathway/stageA_closeout/docrec_out/store"
d = runpy.run_path(DOCREC)
u = json.load(open(os.path.join(HERE, "universe_blind.json")))
want = None
if len(sys.argv) > 2:
    want = [r["accession"] for r in json.load(open(sys.argv[2]))]
seen, bad = [], 0
for i, line in enumerate(open(sys.argv[1], encoding="utf-8"), 1):
    if not line.strip():
        continue
    try:
        rec = json.loads(line)
    except ValueError as e:
        print(f"line {i}: NOT JSON: {e}"); bad += 1; continue
    acc = rec.get("accession") if isinstance(rec, dict) else None
    if acc not in u["studies"]:
        print(f"line {i}: accession {acc!r} not in the sealed tranche"); bad += 1; continue
    if acc in seen:
        print(f"line {i}: duplicate {acc}"); bad += 1
    seen.append(acc)
    us = u["studies"][acc]
    views = {n: d["norm_ws"](open(os.path.join(STORE, acc, n), encoding="utf-8", errors="replace").read())
             for n, m in us["files"].items() if m["usable"]}
    q3 = rec.get("Q3") if isinstance(rec.get("Q3"), dict) else {}
    c = q3.get("count")
    fake = c if (q3.get("scope") != "none" and isinstance(c, int) and not isinstance(c, bool) and c >= 1) else 1
    g, reason, probs, flags = d["grade_record"](rec, fake, us, views, "T1")
    if g == "INVALID":
        bad += 1
        print(f"{acc}: INVALID")
        for p in probs:
            print(f"    - {p}")
    else:
        print(f"{acc}: VALID" + (f"  flags={flags}" if flags else ""))
if want is not None:
    miss = sorted(set(want) - set(seen)); extra = sorted(set(seen) - set(want))
    if miss: print("MISSING records:", miss); bad += 1
    if extra: print("EXTRA records (not in this batch):", extra); bad += 1
print("ALL VALID" if bad == 0 else f"{bad} problem(s)")
