#!/usr/bin/env python3
"""
stageA_universe_audit.py  v0.1.1   Stage A — A17: MetaboLights public-universe verification (listing completeness).

Question (one): does the REST LIST enumerate the complete public accession universe, or a partial / paginated set?
Method (fixed before execution; no study-level file and no value is ever requested):
  1. LIST   {REST}/studies                       R times (R = 3, >= PAUSE s apart): bytes, SHA-256, response headers,
            JSON shape, accession tokens; pagination-/total-like keys detected, never acted on.
  2. FTPROOT {FTP}/                              2 times: names-only HTML index of the public mirror root; folder names
            matching ^MTBLS\\d+$ are the mirror's own public set (the source of record since A14).
  3. Compare: REST set vs mirror set; repeat stability of each; duplicates; non-conforming tokens.
Verdict rules (fixed; §7.5a of the protocol):
  ESTABLISHED_BY_CONCORDANCE  every LIST repeat byte-identical, every FTPROOT repeat identical, REST set == mirror set,
                              no pagination-/total-like key, no duplicates, no non-conforming token.
  BOUNDED                     repeats stable, sets differ: census universe U = REST ∩ mirror; band [|∩|, |∪|];
                              every accession in the symmetric difference is listed by name, carried UNRESOLVED,
                              never screened unless a later dated amendment resolves it.
  UNRESOLVED                  any repeat unstable, any pagination-/total-like key, FTPROOT unavailable, duplicates,
                              or a non-conforming token in either source. Census remains blocked.
Outputs: stageA_universe_manifest.csv (accession, in_rest, in_mirror, status; sorted numerically), stageA_universe_report.json,
         stageA_universe_fetch_log.txt, SHA256SUMS.txt. Exit 0 only on ESTABLISHED_BY_CONCORDANCE or BOUNDED; else 1.
Lock: three hashes (this file, the protocol, stageA_universe_authorisation.txt) + owner_signoff, as for the schema probe.
No redirect is ever followed; a non-200, redirect, oversize or transport failure is a recorded halt (UNRESOLVED).
"""
from __future__ import annotations

import argparse, hashlib, http.client, json, re, socket, sys, time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import HTTPRedirectHandler, Request, build_opener

VERSION = "0.1.1"
PROTOCOL = "EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0.6.9.md"
AUTH_FILE = "stageA_universe_authorisation.txt"
AUTH_LINE = re.compile(r"LIVE UNIVERSE AUDIT AUTHORISED (\d{4}-\d{2}-\d{2}) protocol=(EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0\.\d+(?:\.\d+)?\.md) owner=(\S+)")
REST = "https://www.ebi.ac.uk/metabolights/ws"
FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"
R_LIST, R_ROOT, PAUSE = 3, 2, 2.0
MAX_BYTES, TIMEOUT = 20_000_000, 120
# v0.1.1 (A-5, after A17 run 1 2026-10-25): the mirror root index streamed ~722 KB in ~3 min and the second repeat's read
# timed out at 120 s. FTPROOT gets its own per-read timeout; nothing is retried; a timeout is still a recorded halt.
TIMEOUT_ROOT = 600
ACC = re.compile(r"^MTBLS\d+$")
PAGINATION_KEY = re.compile(r"^(total|totalElements|totalPages|page|pageNumber|pageSize|size|limit|offset|next|previous|last|first|count|numberOfElements|hasNext|hasMore)$", re.I)
HREF = re.compile(r'href="([^"?#]+)"', re.I)
HEADER_KEEP = re.compile(r"^(content-type|content-length|content-range|link|x-total-count|x-pagination.*|etag|last-modified|date|cache-control|server|transfer-encoding|x-ratelimit.*|retry-after)$", re.I)


class Halt(Exception):
    pass


class Refused(Exception):
    pass


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


class Redirect:
    def __init__(self, code: int, location: str):
        self.code, self.location = code, location


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def url_for(kind: str) -> str:
    if kind == "LIST":
        return f"{REST}/studies"
    if kind == "FTPROOT":
        return f"{FTP}/"
    raise Refused(f"REFUSED: unknown request kind {kind!r} (only LIST and FTPROOT exist)")


class Transport:
    def __init__(self, offline: dict | None = None, log=None):
        self.offline, self.log = offline, (log if log is not None else [])
        self.opener = build_opener(_NoRedirect())

    def get(self, kind: str) -> tuple[bytes, dict]:
        url = url_for(kind)
        if self.offline is not None:
            v = self.offline.get(url)
            if callable(v):
                v = v()
            if v is None:
                v = HTTPError(url, 404, "mock 404", None, None)
            if isinstance(v, Redirect):
                self.log.append((url, v.code, 0, f"REDIRECT-NOT-FOLLOWED location={v.location}"))
                raise Halt(f"{url}: HTTP {v.code} redirect not followed")
            if isinstance(v, HTTPError):
                self.log.append((url, v.code, 0, "")); raise Halt(f"{url}: {v}")
            if isinstance(v, Exception):
                self.log.append((url, "ERR", 0, f"TRANSPORT-FAILURE {type(v).__name__}: {str(v)[:120]}"))
                raise Halt(f"{url}: transport failure before any response ({type(v).__name__}: {v})")
            b, hdr = (v if isinstance(v, tuple) else (v, {}))
        else:
            try:
                t0 = time.monotonic()
                with self.opener.open(Request(url, headers={"User-Agent": f"stageA_universe_audit/{VERSION}"}), timeout=(TIMEOUT_ROOT if kind == "FTPROOT" else TIMEOUT)) as r:
                    if r.status != 200:
                        raise Halt(f"{url}: HTTP {r.status}")
                    hdr = {k: str(v)[:200] for k, v in r.headers.items() if HEADER_KEEP.match(k)}
                    cl = r.headers.get("Content-Length")
                    if cl and cl.isdigit() and int(cl) > MAX_BYTES:
                        self.log.append((url, 200, 0, f"OVERSIZE-REFUSED-UNREAD content-length={cl}"))
                        raise Halt(f"{url}: oversize response refused before reading the body")
                    b = r.read(MAX_BYTES + 1)
            except HTTPError as e:
                if e.code in (301, 302, 303, 307, 308):
                    loc = str(e.headers.get("Location", ""))[:200] if e.headers else ""
                    self.log.append((url, e.code, 0, f"REDIRECT-NOT-FOLLOWED location={loc}"))
                    raise Halt(f"{url}: HTTP {e.code} redirect not followed")
                self.log.append((url, e.code, 0, "")); raise Halt(f"{url}: {e}")
            except (URLError, http.client.HTTPException, socket.timeout, OSError) as e:
                self.log.append((url, "ERR", 0, f"TRANSPORT-FAILURE {type(e).__name__}: {str(e)[:120]}"))
                raise Halt(f"{url}: transport failure before any response ({type(e).__name__}: {e})")
        if len(b) > MAX_BYTES:
            raise Halt(f"{url}: exceeds {MAX_BYTES} bytes")
        el = f"{time.monotonic() - t0:.1f}s" if self.offline is None else "offline"
        self.log.append((url, 200, len(b), sha_bytes(b) + f"\telapsed={el}"))
        return b, hdr


# ---------------------------------------------------------------- parsing (shape-strict; nothing inferred)
def walk_keys(o, out: set, depth=0):
    if depth > 6:
        return
    if isinstance(o, dict):
        for k, v in o.items():
            out.add(str(k)); walk_keys(v, out, depth + 1)
    elif isinstance(o, list):
        for v in o[:50]:
            walk_keys(v, out, depth + 1)


def parse_list(b: bytes) -> dict:
    try:
        o = json.loads(b.decode("utf-8"))
    except Exception as e:
        raise Halt(f"LIST: malformed JSON ({type(e).__name__})")
    keys: set = set(); walk_keys(o, keys)
    if isinstance(o, dict) and isinstance(o.get("content"), list):
        items = o["content"]
    elif isinstance(o, list):
        items = o
    else:
        raise Halt(f"LIST: unexpected JSON shape; top-level keys {sorted(keys)[:20]}")
    tokens = [str(x) if not isinstance(x, dict) else str(x.get("accession", x.get("studyIdentifier", x.get("id", "")))) for x in items]
    return {"tokens": tokens, "keys": sorted(keys), "pagination_like_keys": sorted(k for k in keys if PAGINATION_KEY.match(k))}


def parse_root(b: bytes) -> tuple[list[str], list[str]]:
    """(relative folder names, navigation hrefs). v0.1.1 (A-5): an absolute href (starts with '/') or a parent link ('../',
    '?C=' sort links are already excluded by HREF) is the index's own navigation, not a study folder; it is listed, never
    counted and never a 'non-conforming token' (A17 run 1 flagged the parent link 'studies')."""
    html = b.decode("utf-8", "replace")
    hrefs = [h for h in HREF.findall(html) if h.endswith("/")]
    nav = sorted({h for h in hrefs if h.startswith("/") or h.startswith("..") or h.startswith("http")})
    names = [h.rstrip("/").rsplit("/", 1)[-1] for h in hrefs if h not in nav]
    if not names:
        raise Halt("FTPROOT: index returned 200 but no folder href recognised")
    return names, nav


def classify(tokens: list[str]) -> dict:
    conforming = [t for t in tokens if ACC.match(t)]
    return {"n_tokens": len(tokens), "n_conforming": len(conforming), "n_unique": len(set(conforming)),
            "duplicates": sorted({t for t in conforming if conforming.count(t) > 1})[:50],
            "non_conforming": sorted({t for t in tokens if not ACC.match(t)})[:50],
            "min": min((int(t[5:]) for t in conforming), default=None), "max": max((int(t[5:]) for t in conforming), default=None)}


def numkey(a: str) -> int:
    return int(a[5:])


# ---------------------------------------------------------------- audit
def audit(t: Transport, pause: float = PAUSE) -> dict:
    rep = {"version": VERSION, "list_url": url_for("LIST"), "root_url": url_for("FTPROOT"),
           "list_repeats": [], "root_repeats": [], "halted": None, "verdict": "UNRESOLVED", "reasons": []}
    rest_sets, root_sets = [], []
    try:
        for i in range(R_LIST):
            if i and pause:
                time.sleep(pause)
            b, hdr = t.get("LIST"); p = parse_list(b); c = classify(p["tokens"])
            rep["list_repeats"].append({"sha256": sha_bytes(b), "bytes": len(b), "headers": hdr, "json_keys": p["keys"][:40],
                                        "pagination_like_keys": p["pagination_like_keys"], **c})
            rest_sets.append(set(x for x in p["tokens"] if ACC.match(x)))
        for i in range(R_ROOT):
            if pause:
                time.sleep(pause)
            b, hdr = t.get("FTPROOT"); names, nav = parse_root(b); c = classify(names)
            rep["root_repeats"].append({"sha256": sha_bytes(b), "bytes": len(b), "headers": hdr, "navigation_hrefs_ignored": nav, **c})
            root_sets.append(set(x for x in names if ACC.match(x)))
        rest, root = rest_sets[0], root_sets[0]
        both, either = rest & root, rest | root
        rep["comparison"] = {"n_rest": len(rest), "n_mirror": len(root), "n_both": len(both), "n_either": len(either),
                             "rest_only": sorted(rest - root, key=numkey), "mirror_only": sorted(root - rest, key=numkey),
                             "list_repeats_identical_bytes": len({r["sha256"] for r in rep["list_repeats"]}) == 1,
                             "list_repeats_identical_sets": all(s == rest for s in rest_sets),
                             "root_repeats_identical_sets": all(s == root for s in root_sets)}
        cmp = rep["comparison"]; reasons = []
        if not cmp["list_repeats_identical_sets"]:
            reasons.append("LIST accession set differs between repeats")
        if not cmp["root_repeats_identical_sets"]:
            reasons.append("mirror root set differs between repeats")
        if any(r["pagination_like_keys"] for r in rep["list_repeats"]):
            reasons.append("pagination-/total-like key present in LIST: " + ", ".join(sorted({k for r in rep["list_repeats"] for k in r["pagination_like_keys"]})))
        if any(r["duplicates"] for r in rep["list_repeats"] + rep["root_repeats"]):
            reasons.append("duplicate accession in a listing")
        if any(r["non_conforming"] for r in rep["list_repeats"] + rep["root_repeats"]):
            reasons.append("non-conforming token in a listing: " + ", ".join(sorted({x for r in rep["list_repeats"] + rep["root_repeats"] for x in r["non_conforming"]})[:10]))
        if any(h.lower() in ("link", "content-range", "x-total-count") or h.lower().startswith("x-pagination") for r in rep["list_repeats"] for h in r["headers"]):
            reasons.append("pagination-like response header present on LIST")
        rep["reasons"] = reasons
        if reasons:
            rep["verdict"] = "UNRESOLVED"
        elif rest == root:
            rep["verdict"] = "ESTABLISHED_BY_CONCORDANCE"
        else:
            rep["verdict"] = "BOUNDED"
        rep["census_universe"] = {"definition": "REST ∩ mirror (both interfaces list the accession as public)",
                                  "n": len(both), "band": [len(both), len(either)],
                                  "unresolved_accessions": sorted(either - both, key=numkey)} if rep["verdict"] != "UNRESOLVED" else None
        rep["manifest"] = [(a, a in rest, a in root) for a in sorted(either, key=numkey)]
    except (Halt, Refused) as e:
        rep["halted"] = str(e); rep["verdict"] = "UNRESOLVED"; rep["reasons"].append(f"halted: {e}")
        # v0.1.1 (A-5): a halt after at least one LIST and one FTPROOT succeeded still writes the comparison and the
        # manifest, marked PROVISIONAL; the verdict stays UNRESOLVED and no census may read a provisional manifest.
        if rest_sets and root_sets:
            rest, root = rest_sets[0], root_sets[0]; both, either = rest & root, rest | root
            rep["provisional_comparison"] = {"basis": "first LIST and first FTPROOT only; NOT a verdict", "n_rest": len(rest), "n_mirror": len(root),
                                             "n_both": len(both), "n_either": len(either), "rest_only": sorted(rest - root, key=numkey),
                                             "mirror_only": sorted(root - rest, key=numkey)}
            rep["manifest"] = [(a, a in rest, a in root) for a in sorted(either, key=numkey)]; rep["manifest_provisional"] = True
    return rep


def write_outputs(rep: dict, log: list, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    man = rep.pop("manifest", [])
    prov = rep.get("manifest_provisional", False)
    lines = ["accession,in_rest,in_mirror,status"] + [f"{a},{int(r)},{int(m)},{'PROVISIONAL_' if prov else ''}{'UNIVERSE' if r and m else 'UNRESOLVED'}" for a, r, m in man]
    (out / "stageA_universe_manifest.csv").write_text("\n".join(lines) + "\n")
    rep["manifest_rows"] = len(man); rep["manifest_sha256"] = sha_bytes((out / "stageA_universe_manifest.csv").read_bytes())
    rep["run_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    (out / "stageA_universe_report.json").write_text(json.dumps(rep, indent=2) + "\n")
    (out / "stageA_universe_fetch_log.txt").write_text("".join(f"{u}\t{s}\t{n}\t{h}\n" for u, s, n, h in log))
    (out / "SHA256SUMS.txt").write_text("".join(f"{sha_bytes((out / f).read_bytes())}  {f}\n" for f in ("stageA_universe_manifest.csv", "stageA_universe_report.json", "stageA_universe_fetch_log.txt")))


def exit_code_for(rep: dict) -> int:
    return 0 if rep["verdict"] in ("ESTABLISHED_BY_CONCORDANCE", "BOUNDED") else 1


# ---------------------------------------------------------------- lock / authorisation (same contract as the schema probe)
def auth_status(auth_path: Path | None) -> dict:
    if not auth_path or not auth_path.exists():
        return {"present": False, "valid": False, "line": None}
    m = AUTH_LINE.fullmatch(auth_path.read_text().strip())
    return {"present": True, "valid": bool(m), "line": auth_path.read_text().strip()[:200],
            "date": m.group(1) if m else None, "protocol": m.group(2) if m else None, "owner": m.group(3) if m else None}


def sha_file(p: Path | None):
    return sha_bytes(p.read_bytes()) if p and p.exists() else None


def lock_manifest(protocol_path: Path | None, auth_path: Path | None) -> dict:
    here = Path(__file__).resolve()
    return {"version": VERSION, "protocol": PROTOCOL,
            "requests": {"LIST": url_for("LIST"), "LIST_repeats": R_LIST, "FTPROOT": url_for("FTPROOT"), "FTPROOT_repeats": R_ROOT, "pause_s": PAUSE,
                         "nothing_else": "no study folder, no study file, no value; redirects never followed"},
            "verdict_rules": ["ESTABLISHED_BY_CONCORDANCE: stable repeats, REST set == mirror set, no pagination-like key/header, no duplicates/non-conforming tokens",
                              "BOUNDED: stable repeats, sets differ; universe = REST ∩ mirror; band [|∩|,|∪|]; difference listed, UNRESOLVED, never screened",
                              "UNRESOLVED: anything else; census stays blocked"],
            "outputs": ["stageA_universe_manifest.csv", "stageA_universe_report.json", "stageA_universe_fetch_log.txt", "SHA256SUMS.txt"],
            "hashes": {here.name: sha_file(here), PROTOCOL: sha_file(protocol_path), AUTH_FILE: sha_file(auth_path)},
            "authorisation_file": auth_status(auth_path), "owner_signoff": None,
            "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def cmd_plan(a):
    m = lock_manifest(Path(a.protocol) if a.protocol else None, Path(a.auth) if a.auth else None)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    (out / "stageA_universe_lock.json").write_text(json.dumps(m, indent=2) + "\n")
    print(json.dumps(m, indent=2)); print("[plan] written; owner adds owner_signoff before any run")


def cmd_run(a):
    lock = Path(a.out) / "stageA_universe_lock.json"
    if not lock.exists():
        sys.exit("[run] REFUSED: no stageA_universe_lock.json")
    m = json.loads(lock.read_text())
    if not m.get("owner_signoff"):
        sys.exit("[run] REFUSED: lock not signed by owner")
    now = lock_manifest(Path(a.protocol), Path(a.auth))["hashes"]
    if any(v is None for v in now.values()) or now != m.get("hashes"):
        sys.exit(f"[run] REFUSED: hashes differ from the lock or a locked file is missing\n locked={m.get('hashes')}\n now={now}")
    st = auth_status(Path(a.auth))
    if not st["valid"] or st["protocol"] != PROTOCOL or Path(a.protocol).name != PROTOCOL:
        sys.exit(f"[run] REFUSED: authorisation invalid or not naming {PROTOCOL}")
    t = Transport(); rep = audit(t); write_outputs(rep, t.log, Path(a.out))
    print(json.dumps({k: rep[k] for k in ("verdict", "reasons", "comparison", "census_universe", "halted") if k in rep}, indent=2))
    sys.exit(exit_code_for(rep))


# ---------------------------------------------------------------- selftest (offline mock only)
def _mock(n=3584, mirror_extra=(), rest_extra=()):
    accs = [f"MTBLS{i}" for i in range(1, n + 1)]
    rest = {url_for("LIST"): (json.dumps({"content": accs + list(rest_extra)}).encode(), {"Content-Type": "application/json"})}
    root_html = "<html>" + "".join(f'<a href="{a}/">{a}</a>' for a in accs + list(mirror_extra)) + '<a href="README.txt">r</a></html>'
    rest[url_for("FTPROOT")] = (root_html.encode(), {"Content-Type": "text/html"})
    return rest


def cmd_selftest(a):
    ok = True

    def check(c, msg):
        nonlocal ok
        print(("PASS " if c else "FAIL ") + msg); ok &= bool(c)

    for kind in ("FILES", "FETCH", "DOWNLOAD", "MTBLS1"):
        try:
            url_for(kind); check(False, f"guard refuses {kind}")
        except Refused:
            check(True, f"guard refuses {kind} (only LIST and FTPROOT exist)")
    t = Transport(offline=_mock()); r = audit(t, pause=0)
    check(r["verdict"] == "ESTABLISHED_BY_CONCORDANCE" and r["comparison"]["n_rest"] == 3584 and r["census_universe"]["n"] == 3584 and exit_code_for(r) == 0,
          "concordant REST and mirror -> ESTABLISHED_BY_CONCORDANCE, exit 0")
    check(len(t.log) == R_LIST + R_ROOT and all(l[1] == 200 for l in t.log), f"exactly {R_LIST} LIST + {R_ROOT} FTPROOT requests, nothing else")
    r = audit(Transport(offline=_mock(mirror_extra=("MTBLS9001",), rest_extra=("MTBLS9002",))), pause=0)
    check(r["verdict"] == "BOUNDED" and r["census_universe"]["band"] == [3584, 3586] and r["census_universe"]["unresolved_accessions"] == ["MTBLS9001", "MTBLS9002"]
          and r["comparison"]["rest_only"] == ["MTBLS9002"] and r["comparison"]["mirror_only"] == ["MTBLS9001"] and exit_code_for(r) == 0,
          "discordant sets -> BOUNDED; universe = intersection; band and difference listed")
    off = _mock(); off[url_for("LIST")] = (json.dumps({"content": [f"MTBLS{i}" for i in range(1, 3585)], "totalElements": 4000, "page": 0}).encode(), {})
    r = audit(Transport(offline=off), pause=0)
    check(r["verdict"] == "UNRESOLVED" and any("pagination" in x for x in r["reasons"]) and exit_code_for(r) == 1, "pagination-like key -> UNRESOLVED, exit 1")
    off = _mock(); off[url_for("LIST")] = (json.dumps({"content": [f"MTBLS{i}" for i in range(1, 3585)]}).encode(), {"Link": '<...?page=2>; rel="next"'})
    check(audit(Transport(offline=off), pause=0)["verdict"] == "UNRESOLVED", "pagination-like Link header -> UNRESOLVED")
    calls = {"n": 0}
    def unstable():
        calls["n"] += 1; return (json.dumps({"content": [f"MTBLS{i}" for i in range(1, 3585 - (calls["n"] % 2))]}).encode(), {})
    off = _mock(); off[url_for("LIST")] = unstable
    check(audit(Transport(offline=off), pause=0)["verdict"] == "UNRESOLVED", "LIST set differing between repeats -> UNRESOLVED")
    off = _mock(); off[url_for("LIST")] = (json.dumps({"content": ["MTBLS1", "MTBLS1", "MTBLS2"]}).encode(), {})
    check(audit(Transport(offline=off), pause=0)["verdict"] == "UNRESOLVED", "duplicate accession -> UNRESOLVED")
    off = _mock(); off[url_for("LIST")] = (json.dumps({"content": ["MTBLS1", "MTBLS2", "MTBLS3_old"]}).encode(), {})
    r = audit(Transport(offline=off), pause=0); check(r["verdict"] == "UNRESOLVED" and "MTBLS3_old" in " ".join(r["reasons"]), "non-conforming token -> UNRESOLVED, named")
    off = _mock(); off.pop(url_for("FTPROOT"))
    r = audit(Transport(offline=off), pause=0); check(r["verdict"] == "UNRESOLVED" and "404" in r["halted"], "FTPROOT 404 -> halt, UNRESOLVED")
    off = _mock(); off[url_for("FTPROOT")] = Redirect(301, FTP)
    tr = Transport(offline=off); r = audit(tr, pause=0)
    check(r["verdict"] == "UNRESOLVED" and "redirect" in r["halted"] and not any(l[0] == FTP for l in tr.log), "FTPROOT redirect -> halt, destination never requested")
    off = _mock(); off[url_for("LIST")] = http.client.RemoteDisconnected("x")
    r = audit(Transport(offline=off), pause=0); check("transport failure" in r["halted"] and exit_code_for(r) == 1, "transport failure -> recorded halt, exit 1")
    off = _mock(); off[url_for("LIST")] = (b"{not json", {})
    check("malformed" in audit(Transport(offline=off), pause=0)["halted"], "malformed LIST -> halt")
    off = _mock(); off[url_for("FTPROOT")] = (off[url_for("FTPROOT")][0].replace(b"<html>", b'<html><a href="/pub/databases/metabolights/studies/">Parent</a><a href="../">up</a>'), {})
    r = audit(Transport(offline=off), pause=0)
    check(r["verdict"] == "ESTABLISHED_BY_CONCORDANCE" and r["root_repeats"][0]["navigation_hrefs_ignored"] == ["../", "/pub/databases/metabolights/studies/"]
          and r["root_repeats"][0]["non_conforming"] == [], "A-5: parent/absolute hrefs are navigation, listed, never counted or flagged")
    off = _mock(mirror_extra=("MTBLS9001",)); calls2 = {"n": 0}; good = off[url_for("FTPROOT")]
    def root_then_timeout():
        calls2["n"] += 1
        if calls2["n"] == 1:
            return good
        raise_ = socket.timeout("The read operation timed out"); return raise_
    off[url_for("FTPROOT")] = root_then_timeout
    tt = Transport(offline=off); r = audit(tt, pause=0)
    check(r["verdict"] == "UNRESOLVED" and "transport failure" in r["halted"] and r["provisional_comparison"]["mirror_only"] == ["MTBLS9001"]
          and r.get("manifest_provisional") and exit_code_for(r) == 1, "A-5: second root repeat times out -> UNRESOLVED, provisional comparison kept, exit 1")
    import tempfile as _tf
    dp = Path(_tf.mkdtemp()); write_outputs(r, tt.log, dp); mp = (dp / "stageA_universe_manifest.csv").read_text().splitlines()
    check(mp[1] == "MTBLS1,1,1,PROVISIONAL_UNIVERSE" and mp[-1] == "MTBLS9001,0,1,PROVISIONAL_UNRESOLVED", "A-5: provisional manifest rows are marked PROVISIONAL")
    off = _mock(); off[url_for("FTPROOT")] = (b"<html>Index of /</html>", {})
    check("no folder href" in audit(Transport(offline=off), pause=0)["halted"], "root index without folder hrefs -> halt, nothing inferred")
    import tempfile
    d = Path(tempfile.mkdtemp()); t = Transport(offline=_mock(mirror_extra=("MTBLS9001",))); r = audit(t, pause=0); write_outputs(r, t.log, d)
    man = (d / "stageA_universe_manifest.csv").read_text().splitlines()
    check(man[0] == "accession,in_rest,in_mirror,status" and man[1] == "MTBLS1,1,1,UNIVERSE" and man[-1] == "MTBLS9001,0,1,UNRESOLVED" and len(man) == 3586
          and (d / "SHA256SUMS.txt").exists() and json.loads((d / "stageA_universe_report.json").read_text())["manifest_rows"] == 3585,
          "manifest sorted numerically with status; SHA256SUMS written; report carries manifest hash and run timestamp")
    check(auth_status(None)["valid"] is False, "authorisation absent -> invalid")
    ap = d / AUTH_FILE; ap.write_text(f"LIVE UNIVERSE AUDIT AUTHORISED 2026-10-25 protocol={PROTOCOL} owner=TEST\n")
    check(auth_status(ap)["valid"] and auth_status(ap)["protocol"] == PROTOCOL, "authorisation line recognised")
    ap.write_text(f"LIVE UNIVERSE AUDIT AUTHORISED 2026-10-25 protocol={PROTOCOL} owner=TEST\nextra\n")
    check(not auth_status(ap)["valid"], "second line rejects the authorisation")
    pf = d / PROTOCOL; pf.write_text("x"); ap.write_text(f"LIVE UNIVERSE AUDIT AUTHORISED 2026-10-25 protocol={PROTOCOL} owner=TEST\n")
    m1 = lock_manifest(pf, ap); check(m1["owner_signoff"] is None and all(m1["hashes"].values()), "lock manifest: three hashes, unsigned")
    pf.write_text("y"); check(lock_manifest(pf, ap)["hashes"] != m1["hashes"], "any byte change after lock is detected")
    print("SELFTEST PASSED" if ok else "SELFTEST FAILED"); sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    for name in ("plan", "run"):
        p = sub.add_parser(name)
        p.add_argument("--protocol", default=PROTOCOL); p.add_argument("--auth", default=AUTH_FILE); p.add_argument("--out", default="stageA_universe_out")
    a = ap.parse_args()
    {"selftest": cmd_selftest, "plan": cmd_plan, "run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
