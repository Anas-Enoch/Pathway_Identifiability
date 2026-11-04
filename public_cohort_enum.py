#!/usr/bin/env python3
"""
public_cohort_enum.py  v0.1.0   Discovery-slot enumeration for PUBLIC_COHORT_SCREEN_v0.1.md s2 (frozen f53d8b6).

Six families, three slots each, metadata only. For every family it records the UTC timestamp, endpoint, exact
search expression, total hits, the archived raw index response, the deduplicated canonical ids, their SHA-256 sort
keys 'PUBLIC-COHORT-v0.1|<family>|<canonical-id>', and the first three. A family whose index cannot be enumerated
through the endpoint below is marked UNSEARCHABLE with the response archived; its slots are not refilled.
It never requests a data file: only index/search endpoints and, for MetaboLights, the i_Investigation.txt of a
new accession (the Stage A reading view) to read 'Study Public Release Date'.

  python3 -I public_cohort_enum.py selftest
  python3 -I public_cohort_enum.py run --out screen_enum [--family dryad ...] \
        --known-mtbls <file with the Stage A census accessions> --known-st <file with the B0 universe study ids>
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path

VERSION = "0.1.0"
SALT = "PUBLIC-COHORT-v0.1|"
SLOTS = 3
MTBLS_AFTER = "2026-10-28"      # Stage A census date (protocol s2)
ST_AFTER = "2026-10-17"         # B0 census date (protocol s2)
UA = "public_cohort_enum/" + VERSION
TIMEOUT = 120
EXPR = {
    "dryad":    "(metabolomics OR lipidomics) AND (human OR participants OR cohort)",
    "zenodo":   "(metabolomics OR lipidomics) AND (human OR participants OR cohort)",
    "figshare": ":title: metabolomics",
    "ngdc":     "metabolomics human",
    "metabolights": f"public accessions absent from the Stage A census universe, release date > {MTBLS_AFTER}",
    "workbench":    f"study ids absent from the B0 universe, summary release_date > {ST_AFTER}",
}
HUMAN = re.compile(r"\b(human|participants?|patients?|subjects?|cohort|volunteers?|homo sapiens)\b", re.I)
METAB = re.compile(r"\b(metabolom\w*|lipidom\w*)\b", re.I)


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def key(family, cid):
    return hashlib.sha256((SALT + family + "|" + cid).encode()).hexdigest()


def http(url, data=None, headers=None):
    req = urllib.request.Request(url, data=data, headers={"User-Agent": UA, **(headers or {})})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return r.status, r.read()


def now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


# ---------------------------------------------------------------- family enumerators: return (records, raw_pages)
# each record: {"id": canonical id, "title": ..., "url": ..., "text": title+abstract for the topical filter, "date": ...}
def enum_dryad(fetch):
    recs, pages = [], []
    q = urllib.parse.quote(EXPR["dryad"])
    page = 1
    while True:
        url = f"https://datadryad.org/api/v2/search?q={q}&per_page=100&page={page}"
        st, raw = fetch(url); pages.append({"url": url, "status": st, "sha256": sha_bytes(raw), "bytes": len(raw)})
        js = json.loads(raw)
        items = js.get("_embedded", {}).get("stash:datasets", [])
        for d in items:
            doi = d.get("identifier", "")
            recs.append({"id": doi, "title": d.get("title", ""), "url": "https://datadryad.org/dataset/" + urllib.parse.quote(doi, safe=""),
                         "text": d.get("title", "") + " " + d.get("abstract", ""), "date": d.get("publicationDate", "")})
        total = js.get("total", 0)
        if page * 100 >= total or not items or page >= 10:
            break
        page += 1
    return recs, pages, total


def enum_zenodo(fetch):
    recs, pages = [], []
    q = urllib.parse.quote(EXPR["zenodo"])
    page, total = 1, 0
    while True:
        url = f"https://zenodo.org/api/records?q={q}&type=dataset&size=100&page={page}"
        st, raw = fetch(url); pages.append({"url": url, "status": st, "sha256": sha_bytes(raw), "bytes": len(raw)})
        js = json.loads(raw)
        hits = js.get("hits", {}).get("hits", [])
        total = js.get("hits", {}).get("total", 0)
        for h in hits:
            md = h.get("metadata", {})
            doi = h.get("doi") or md.get("doi") or f"zenodo:{h.get('id')}"
            recs.append({"id": doi, "title": md.get("title", ""), "url": h.get("links", {}).get("self_html") or f"https://zenodo.org/records/{h.get('id')}",
                         "text": md.get("title", "") + " " + re.sub("<[^>]+>", " ", md.get("description", "")), "date": md.get("publication_date", "")})
        if page * 100 >= total or not hits or page >= 10:
            break
        page += 1
    return recs, pages, total


def enum_figshare(fetch):
    recs, pages = [], []
    page, total = 1, 0
    while True:
        url = "https://api.figshare.com/v2/articles/search"
        body = json.dumps({"search_for": EXPR["figshare"], "item_type": 3, "page_size": 100, "page": page}).encode()
        st, raw = fetch(url, data=body, headers={"Content-Type": "application/json"})
        pages.append({"url": url + " [POST]", "body": body.decode(), "status": st, "sha256": sha_bytes(raw), "bytes": len(raw)})
        items = json.loads(raw)
        for a in items:
            doi = a.get("doi") or f"figshare:{a.get('id')}"
            recs.append({"id": doi, "title": a.get("title", ""), "url": a.get("url_public_html") or a.get("url", ""),
                         "text": a.get("title", ""), "date": a.get("published_date", "")})
        total += len(items)
        if len(items) < 100 or page >= 10:
            break
        page += 1
    return recs, pages, total


def enum_ngdc(fetch):
    """NGDC OMIX has no documented public search API known to this script. One attempt at the site's list endpoint;
    anything that is not JSON with records marks the family UNSEARCHABLE (protocol s2)."""
    url = "https://ngdc.cncb.ac.cn/omix/api/omix/list?keyword=" + urllib.parse.quote(EXPR["ngdc"]) + "&pageSize=100&pageNum=1"
    st, raw = fetch(url)
    page = {"url": url, "status": st, "sha256": sha_bytes(raw), "bytes": len(raw)}
    try:
        js = json.loads(raw)
    except Exception:
        return None, [page], 0
    items = js.get("data", {}).get("list") if isinstance(js.get("data"), dict) else js.get("data")
    if not isinstance(items, list):
        return None, [page], 0
    recs = [{"id": it.get("accession") or it.get("omixAccession") or "", "title": it.get("title", ""),
             "url": "https://ngdc.cncb.ac.cn/omix/release/" + (it.get("accession") or it.get("omixAccession") or ""),
             "text": (it.get("title", "") + " " + it.get("description", "")), "date": it.get("releaseDate", "")} for it in items]
    return recs, [page], len(recs)


def enum_metabolights(fetch, known):
    url = "https://www.ebi.ac.uk/metabolights/ws/studies"
    st, raw = fetch(url); pages = [{"url": url, "status": st, "sha256": sha_bytes(raw), "bytes": len(raw)}]
    js = json.loads(raw)
    acc = js.get("content") if isinstance(js, dict) else js
    acc = [a if isinstance(a, str) else a.get("accession", "") for a in acc]
    new = sorted(a for a in acc if a.startswith("MTBLS") and a not in known)
    recs = []
    for a in new:
        iurl = f"https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/{a}/i_Investigation.txt"
        try:
            st2, raw2 = fetch(iurl)
        except Exception as e:  # noqa: BLE001
            pages.append({"url": iurl, "status": "error", "error": str(e)[:120]}); continue
        pages.append({"url": iurl, "status": st2, "sha256": sha_bytes(raw2), "bytes": len(raw2)})
        txt = raw2.decode("utf-8", "replace")
        m = re.search(r'^Study Public Release Date\t"?([0-9-]+)"?', txt, re.M)
        t = re.search(r'^Study Title\t"?(.*?)"?$', txt, re.M)
        d = re.search(r'^Study Description\t"?(.*?)"?$', txt, re.M)
        rel = m.group(1) if m else ""
        if rel > MTBLS_AFTER:
            recs.append({"id": a, "title": t.group(1) if t else "", "url": f"https://www.ebi.ac.uk/metabolights/{a}",
                         "text": (t.group(1) if t else "") + " " + (d.group(1) if d else ""), "date": rel})
    return recs, pages, len(new)


def enum_workbench(fetch, known):
    url = "https://www.metabolomicsworkbench.org/rest/study/study_id/ST/species"
    st, raw = fetch(url); pages = [{"url": url, "status": st, "sha256": sha_bytes(raw), "bytes": len(raw)}]
    js = json.loads(raw)
    rows = list(js.values()) if isinstance(js, dict) else js
    ids = sorted({str(r.get("Study ID") or r.get("study_id") or "") for r in rows if isinstance(r, dict)})
    new = [s for s in ids if s.startswith("ST") and s not in known]
    recs = []
    for s in new:
        surl = f"https://www.metabolomicsworkbench.org/rest/study/study_id/{s}/summary"
        try:
            st2, raw2 = fetch(surl)
        except Exception as e:  # noqa: BLE001
            pages.append({"url": surl, "status": "error", "error": str(e)[:120]}); continue
        pages.append({"url": surl, "status": st2, "sha256": sha_bytes(raw2), "bytes": len(raw2)})
        try:
            sm = json.loads(raw2)
        except Exception:
            continue
        sm = list(sm.values())[0] if isinstance(sm, dict) and sm and all(isinstance(v, dict) for v in sm.values()) else sm
        rel = str(sm.get("release_date", "") or sm.get("Release Date", ""))
        if rel[:10] > ST_AFTER:
            recs.append({"id": s, "title": sm.get("study_title", ""), "url": f"https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID={s}",
                         "text": sm.get("study_title", "") + " " + sm.get("study_summary", "") + " " + sm.get("species", "") if isinstance(sm.get("species", ""), str) else sm.get("study_title", ""),
                         "date": rel})
    return recs, pages, len(new)


FAMILIES = ("ngdc", "dryad", "zenodo", "figshare", "metabolights", "workbench")


def select(family, recs):
    """Topical filter (metabolomics/lipidomics AND human/cohort terms), dedupe by canonical id, hash sort, first SLOTS."""
    seen, elig = set(), []
    for r in recs:
        cid = (r["id"] or "").strip().lower()
        if not cid or cid in seen:
            continue
        seen.add(cid)
        if METAB.search(r["text"]) and HUMAN.search(r["text"]):
            elig.append(dict(r, id=cid, sort_key=key(family, cid)))
    elig.sort(key=lambda r: (r["sort_key"], r["id"]))
    return elig, elig[:SLOTS]


def run_family(family, fetch, known_mtbls, known_st):
    rec = {"family": family, "expression": EXPR[family], "timestamp_utc": now(), "version": VERSION}
    try:
        if family == "dryad": recs, pages, total = enum_dryad(fetch)
        elif family == "zenodo": recs, pages, total = enum_zenodo(fetch)
        elif family == "figshare": recs, pages, total = enum_figshare(fetch)
        elif family == "ngdc": recs, pages, total = enum_ngdc(fetch)
        elif family == "metabolights": recs, pages, total = enum_metabolights(fetch, known_mtbls)
        else: recs, pages, total = enum_workbench(fetch, known_st)
    except Exception as e:  # noqa: BLE001
        rec.update(status="UNSEARCHABLE", error=f"{type(e).__name__}: {str(e)[:200]}", pages=[])
        return rec
    rec["pages"] = pages
    if recs is None:
        rec.update(status="UNSEARCHABLE", reason="index endpoint returned no parseable record list")
        return rec
    elig, chosen = select(family, recs)
    rec.update(status="ENUMERATED", total_hits=total, indexed_records=len(recs), topical_eligible=len(elig),
               ordered_ids=[(r["sort_key"][:16], r["id"]) for r in elig],
               selected=[{k: r[k] for k in ("id", "title", "url", "date", "sort_key")} for r in chosen])
    return rec


def cmd_run(a):
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    known_mtbls = set(re.findall(r"MTBLS\d+", Path(a.known_mtbls).read_text())) if a.known_mtbls else set()
    known_st = set(re.findall(r"ST\d{6}", Path(a.known_st).read_text())) if a.known_st else set()
    fams = a.family or list(FAMILIES)
    archive = {}

    def fetch(url, data=None, headers=None):
        st, raw = http(url, data, headers)
        archive[sha_bytes(raw)] = raw
        return st, raw
    results = []
    for f in fams:
        if f in ("metabolights", "workbench") and not (known_mtbls if f == "metabolights" else known_st):
            results.append({"family": f, "status": "UNSEARCHABLE", "reason": "known-universe file not supplied; new accessions cannot be defined",
                            "timestamp_utc": now(), "expression": EXPR[f]}); continue
        r = run_family(f, fetch, known_mtbls, known_st)
        results.append(r)
        (out / f"screen_enum_{f}.json").write_text(json.dumps(r, indent=2) + "\n")
        print(f"{f}: {r['status']}", json.dumps({k: r.get(k) for k in ('total_hits', 'indexed_records', 'topical_eligible')}),
              [s["id"] for s in r.get("selected", [])])
    (out / "raw").mkdir(exist_ok=True)
    for h, raw in archive.items():
        (out / "raw" / f"{h}.bin").write_bytes(raw)
    man = {"script": "public_cohort_enum.py", "version": VERSION, "script_sha256": sha_bytes(Path(__file__).read_bytes()),
           "families": fams, "created_utc": now(), "results": [{k: v for k, v in r.items() if k != "pages"} for r in results]}
    (out / "screen_enum_manifest.json").write_text(json.dumps(man, indent=2) + "\n")
    print("manifest", sha_bytes((out / "screen_enum_manifest.json").read_bytes()))


# ---------------------------------------------------------------- selftest with mocked responses
def cmd_selftest(a):
    ok = True

    def check(c, m):
        nonlocal ok
        print(("PASS " if c else "FAIL ") + m); ok &= bool(c)
    # hash order deterministic and independent of input order
    recs = [{"id": f"10.5061/dryad.{i}", "title": "Human plasma metabolomics cohort", "url": "", "text": "Human plasma metabolomics cohort", "date": ""} for i in range(10)]
    e1, c1 = select("dryad", recs); e2, c2 = select("dryad", list(reversed(recs)))
    check([r["id"] for r in c1] == [r["id"] for r in c2] and len(c1) == 3, "hash-ordered selection deterministic, three slots")
    # topical filter: non-human and non-metabolomics records excluded; duplicates collapsed
    recs2 = recs[:2] + [{"id": "x1", "title": "Mouse liver metabolomics", "url": "", "text": "mouse liver metabolomics", "date": ""},
                        {"id": "x2", "title": "Human cohort proteomics", "url": "", "text": "human cohort proteomics", "date": ""},
                        dict(recs[0])]
    e3, _ = select("dryad", recs2)
    check(len(e3) == 2, "topical filter keeps only human metabolomics/lipidomics; duplicate id collapsed")
    # mocked endpoints
    def mock(url, data=None, headers=None):
        if "datadryad" in url:
            return 200, json.dumps({"total": 2, "_embedded": {"stash:datasets": [
                {"identifier": "doi:10.5061/dryad.a", "title": "Human serum metabolomics", "abstract": "cohort", "publicationDate": "2026-09-01"},
                {"identifier": "doi:10.5061/dryad.b", "title": "Mouse metabolomics", "abstract": "", "publicationDate": "2026-09-01"}]}}).encode()
        if "ngdc" in url:
            return 200, b"<html>not json</html>"
        if "ws/studies" in url:
            return 200, json.dumps({"content": ["MTBLS1", "MTBLS99999"]}).encode()
        if "i_Investigation" in url:
            return 200, b'Study Title\t"Human urine metabolomics"\nStudy Description\t"cohort"\nStudy Public Release Date\t2026-11-01\n'
        if "species" in url:
            return 200, json.dumps({"1": {"Study ID": "ST000001"}, "2": {"Study ID": "ST009999"}}).encode()
        if "summary" in url:
            return 200, json.dumps({"study_id": "ST009999", "study_title": "Human plasma metabolomics", "study_summary": "cohort", "release_date": "2026-10-20"}).encode()
        raise RuntimeError("unexpected url " + url)
    r = run_family("dryad", mock, set(), set()); check(r["status"] == "ENUMERATED" and [s["id"] for s in r["selected"]] == ["doi:10.5061/dryad.a"], "dryad mock: one topical record selected")
    r = run_family("ngdc", mock, set(), set()); check(r["status"] == "UNSEARCHABLE", "ngdc mock: non-JSON index -> UNSEARCHABLE")
    r = run_family("metabolights", mock, {"MTBLS1"}, set()); check(r["status"] == "ENUMERATED" and [s["id"] for s in r["selected"]] == ["mtbls99999"], "metabolights mock: new accession with release date after census selected")
    r = run_family("workbench", mock, set(), {"ST000001"}); check([s["id"] for s in r["selected"]] == ["st009999"], "workbench mock: new study after B0 date selected")
    def mock_fail(url, data=None, headers=None):
        raise RuntimeError("connection refused")
    r = run_family("zenodo", mock_fail, set(), set()); check(r["status"] == "UNSEARCHABLE", "network failure -> UNSEARCHABLE, not substituted")
    print("SELFTEST", "PASSED" if ok else "FAILED"); sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    r = sub.add_parser("run"); r.add_argument("--out", default="screen_enum"); r.add_argument("--family", action="append", choices=FAMILIES)
    r.add_argument("--known-mtbls"); r.add_argument("--known-st")
    a = ap.parse_args()
    cmd_selftest(a) if a.cmd == "selftest" else cmd_run(a)


if __name__ == "__main__":
    main()
