#!/usr/bin/env python3
"""docrec.py — documentary recovery (§3.7) executor for DOC_RECOVERY_PROTOCOL_v0.2.md.

Standard library only. Run with `python3 -I docrec.py <command>` from the folder that holds this file, the protocol
and stageA_documentary_candidates_v1.0.csv.

Commands
  selftest                       offline checks (no network)
  plan  --salt-option a|b [--beacon-not-before ISO] --abstract yes|no --manual yes|no
                                 write docrec_lock.json (unsigned)
  order --salt SALT              dry-run of the hash order (no network)
  fetch --tranche T1|A35         owner-run, network; needs the signed lock and an authorisation line
  register --tranche T           register owner-supplied texts (DR-9) before sealing
  seal --tranche T               write the evidence universe (the owner signs it before extraction)
  grade --tranche T              verify quotes, grade, apply adjudication, write results and report

Access boundary: MetaboLights public mirror (directory index, i_ files; s_/a_ files for A35 only) and Europe PMC REST
(search, PMC fullTextXML). Redirects are never followed; bodies shorter than Content-Length are refused; every request is
journaled before the next one starts.
"""
import argparse
import csv
import hashlib
import html.entities
import io
import json
import math
import os
import re
import sys
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

VERSION = "0.1.1"
PROTOCOL = "DOC_RECOVERY_PROTOCOL_v0.2.md"
CANDIDATES = "stageA_documentary_candidates_v1.0.csv"
CANDIDATES_SHA = "c1b5a4893356c517b64ed51c8777c372c52dfbc4210bc86998e7b63be18ab0d7"
LOCK = "docrec_lock.json"
AUTH_FILE = "docrec_authorisation.txt"
OUTDIR = "docrec_out"
T1_SIZE = 40
S3_MAX = 5
S3_CAP = 1000
MAX_BYTES = 20 * 1024 * 1024
TIMEOUT = 120
PAUSE = 0.5
OUTAGE_K = 5          # consecutive indeterminate failures on one host -> stop (exit 2)
MAX_ATTEMPTS = 3      # a study still indeterminate after 3 fetch invocations is finalised as-is
MIN_QUOTE_TOKENS = 3
POOL_M = 371
FUTILITY_POOL = 94
ALPHA = 0.05

MTBLS = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/"
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest/"
BEACON = "https://beacon.nist.gov/beacon/2.0/pulse/time/"  # "pulse at a specific (Unix ms) time, or next closest"
ACC = r"MTBLS[1-9]\d*"
SALT_RX = re.compile(r"^(stageA-docrec-v1|NISTv2:([1-9]\d*)-([1-9]\d*):([0-9A-F]{128}))$")
FETCH_AUTH_RX = re.compile(r"^DOCREC FETCH AUTHORISED (\d{4}-\d{2}-\d{2}) tranche=(T1|A35) candidates=([0-9a-f]{64}) "
                           r"protocol=(\S+) lock=([0-9a-f]{64}) salt=(\S+) owner=(\S+)$")
SEAL_AUTH_RX = re.compile(r"^DOCREC SEAL SIGNED (\d{4}-\d{2}-\d{2}) tranche=(T1|A35) universe=([0-9a-f]{64}) owner=(\S+)$")
UTC_RX = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")
VIEW_MARKERS = ("[TITLE]", "[ABSTRACT]", "[SECTION:", "[TABLE:", "[FRONT]", "[METHODS_PDF]", "[VIEW_HEURISTIC_PDF]")

ALLOWED = [
    ("index", re.compile(r"^" + re.escape(MTBLS) + r"(" + ACC + r")/$")),
    ("isa_i", re.compile(r"^" + re.escape(MTBLS) + r"(" + ACC + r")/i_[^/?#]+\.txt$")),
    ("isa_sa", re.compile(r"^" + re.escape(MTBLS) + r"(" + ACC + r")/[sa]_[^/?#]+\.txt$")),
    ("search", re.compile(r"^" + re.escape(EPMC) + r"search\?[^#/]+$")),
    ("fulltext", re.compile(r"^" + re.escape(EPMC) + r"PMC[1-9]\d*/fullTextXML$")),
]


class Refused(Exception):
    pass


class Outage(Exception):
    pass


def now():
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha_file(p):
    with open(p, "rb") as f:
        return sha_bytes(f.read())


def norm_ws(s):
    return re.sub(r"\s+", " ", unicodedata.normalize("NFC", s or "")).strip()


def one_line(s, n=160):
    return re.sub(r"\s+", " ", str(s))[:n]


def guard(url, allow_sa=False):
    """Return the URL class or raise Refused. s_/a_ files only when allow_sa (stratum A)."""
    path = urllib.parse.urlsplit(url).path  # traversal is a property of the path; queries may carry %2F (DOIs, cursors)
    low = path.lower()
    if ".." in path or "%2f" in low or "%2e%2e" in low or "\\" in path:
        raise Refused(f"path traversal: {url}")
    for kind, rx in ALLOWED:
        if rx.match(url):
            if kind == "isa_sa" and not allow_sa:
                raise Refused(f"s_/a_ file outside stratum A: {url}")
            return kind
    raise Refused(f"not on the allow-list: {url}")


def mirror_url(acc, name=None):
    if not re.fullmatch(ACC, acc):
        raise Refused(f"bad accession {acc}")
    if name is None:
        return f"{MTBLS}{acc}/"
    if "/" in name or name.startswith("."):
        raise Refused(f"bad file name {name!r}")
    return f"{MTBLS}{acc}/{urllib.parse.quote(name)}"


def search_url(query, page_size, cursor="*"):
    params = [("query", query), ("resultType", "core"), ("synonym", "FALSE"), ("format", "json"),
              ("pageSize", str(page_size)), ("cursorMark", cursor)]
    return EPMC + "search?" + urllib.parse.urlencode(params, quote_via=urllib.parse.quote)


def fulltext_url(pmcid):
    if not re.fullmatch(r"PMC[1-9]\d*", pmcid or ""):
        raise Refused(f"bad PMCID {pmcid!r}")
    return f"{EPMC}{pmcid}/fullTextXML"


# ----------------------------------------------------------------------------------------------- transport
class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None  # urllib then raises HTTPError with the 3xx code: recorded, never followed


_TEST_OPENER = None  # selftest only; always None in a real run
_TEST_SLEEP = None


class Transport:
    """get() returns (ok, status, body, note, definite). definite = the outcome is a property of the resource
    (200, redirect, 400, 404, 410, oversize); indeterminate = everything else (transport failure, 5xx, 429, 401/403/407/408…,
    short body)."""

    def __init__(self, log_path, opener=None, sleep=None):
        self.opener = opener or _TEST_OPENER or urllib.request.build_opener(NoRedirect)
        self.sleep = sleep or _TEST_SLEEP or time.sleep
        self.fail_run = {}
        new = not os.path.exists(log_path)
        self.log = open(log_path, "a", encoding="utf-8")
        if new:
            self.log.write("utc\turl\tstatus\tbytes\tcontent_length\telapsed_s\tnote\n")
            self.log.flush()

    def header(self, text):
        self.log.write(f"# {text}\n")
        self.log.flush()

    def _rec(self, url, status, nbytes, cl, elapsed, note):
        self.log.write(f"{now()}\t{url}\t{status}\t{nbytes}\t{cl}\t{elapsed:.2f}\t{note}\n")
        self.log.flush()

    def get(self, url, allow_sa=False):
        guard(url, allow_sa)
        host = urllib.parse.urlsplit(url).netloc
        status, body, note = None, b"", ""
        for attempt in (1, 2):
            t0 = time.time()
            status, body, note, cl = None, b"", "", ""
            try:
                req = urllib.request.Request(url, headers={"User-Agent": f"docrec/{VERSION} (value-blind documentary recovery)"})
                with self.opener.open(req, timeout=TIMEOUT) as r:
                    status = getattr(r, "status", None) or r.getcode()
                    cl = r.headers.get("Content-Length") or ""
                    body = r.read(MAX_BYTES + 1)
                if len(body) > MAX_BYTES:
                    note, body = "OVERSIZE", b""
                elif cl.isdigit() and int(cl) != len(body):
                    note = f"SHORT_BODY({len(body)}/{cl})"
            except urllib.error.HTTPError as e:
                status, note = e.code, f"HTTP_{e.code}"
            except Exception as e:  # noqa: BLE001 - every transport failure is recorded, never raised
                note = f"TRANSPORT_FAILED:{type(e).__name__}:{one_line(e)}"
            self._rec(url, status, len(body), cl, time.time() - t0, note)
            retry = note.startswith("SHORT_BODY") or note.startswith("TRANSPORT_FAILED") or status in (429, 500, 502, 503, 504)
            if retry and attempt == 1:
                self.sleep(PAUSE * 4)
                continue
            break
        ok = status == 200 and note == ""
        # definite: a property of the resource. 401/403/407/408 etc. may be a firewall, proxy or timeout: indeterminate.
        definite = ok or note == "OVERSIZE" or status in (301, 302, 303, 307, 308, 400, 404, 410)
        self.fail_run[host] = 0 if definite else self.fail_run.get(host, 0) + 1
        self.sleep(PAUSE)
        if self.fail_run[host] >= OUTAGE_K:
            raise Outage(f"{OUTAGE_K} consecutive indeterminate failures on {host}")
        return ok, status, (body if ok else b""), note, definite

    def close(self):
        self.log.close()


# ----------------------------------------------------------------------------------------------- parsing
def index_names(html_bytes):
    names = []
    for href in re.findall(r"""href=["']([^"']+)["']""", html_bytes.decode("utf-8", "replace")):
        n = urllib.parse.unquote(href)
        if "/" in n.rstrip("/") or n.startswith("?") or n.startswith("."):
            continue
        names.append(n)
    return names


def norm_doi(v):
    v = (v or "").strip().strip('"').strip()
    v = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:\s*)", "", v, flags=re.I).strip().lower()
    return v if re.fullmatch(r"10\.\d{4,9}/\S+", v) else None


def parse_investigation(text):
    """Publication identifiers from Study- and Investigation-level ISA fields (§3 S2), plus file names (A35)."""
    out = {"pmids": [], "dois": [], "study_files": [], "assay_files": []}
    for row in csv.reader(io.StringIO(text), delimiter="\t"):
        if not row:
            continue
        key = row[0].strip().strip('"')
        vals = [v.strip().strip('"').strip() for v in row[1:]]
        vals = [v for v in vals if v]
        if key in ("Study PubMed ID", "Investigation PubMed ID"):
            out["pmids"] += [v for v in vals if re.fullmatch(r"[1-9]\d{0,8}", v)]
        elif key in ("Study Publication DOI", "Investigation Publication DOI"):
            out["dois"] += [d for d in (norm_doi(v) for v in vals) if d]
        elif key == "Study File Name":
            out["study_files"] += vals
        elif key == "Study Assay File Name":
            out["assay_files"] += vals
    return {k: list(dict.fromkeys(v)) for k, v in out.items()}


def art_ids(r):
    ids = set()
    if r.get("doi"):
        ids.add("doi:" + str(r["doi"]).lower())
    if r.get("pmid"):
        ids.add("pmid:" + str(r["pmid"]))
    if r.get("pmcid"):
        ids.add("pmcid:" + str(r["pmcid"]).upper())
    ids.add(f"{r.get('source', '')}:{r.get('id', '')}")
    return ids


def art_slim(r):
    return {k: r.get(k) for k in ("id", "source", "pmid", "pmcid", "doi", "title", "firstPublicationDate",
                                  "isOpenAccess", "inEPMC")}


def _sort_key(r):
    return (r.get("firstPublicationDate") or "9999-99-99", str(r.get("source") or ""), str(r.get("id") or ""))


def s3_select(results, s2_ids):
    """Sort, then drop any record sharing an identifier (DOI, PMID, PMCID, source:id) with S2 or an earlier record;
    keep the first S3_MAX."""
    seen, uniq = set(s2_ids), []
    for r in sorted((r for r in results if isinstance(r, dict)), key=_sort_key):
        ids = art_ids(r)
        dup = bool(ids & seen)
        seen |= ids
        if not dup:
            uniq.append(r)
    return uniq[:S3_MAX]


def accession_in_text(acc, text):
    return re.search(r"(?<![A-Za-z0-9])" + re.escape(acc) + r"(?![A-Za-z0-9])", text) is not None


# ----------------------------------------------------------------------------------------------- reading views
# Top-level sections enter by a methods-type title at the START of the title (or a methods sec-type); their
# descendants inherit inclusion. Any section whose title or type matches EXCL_RX is dropped with all descendants.
METH_ANY_RX = re.compile(r"\bmethods?\b|\bmethodology\b", re.I)  # "Research design and methods", "STAR★Methods"
METH_START_RX = re.compile(r"^\s*(\d+(\.\d+)*\.?\s*|[IVX]+\.\s*)?(online\s+)?(materials?|methods?|methodology|experimental|research\s+design|"
                           r"patients?|participants?|subjects?|volunteers?|study\s+(design|population|participants|subjects|"
                           r"cohorts?|setting)|cohorts?|populations?|recruitment|enrol?lment|sample\s+(collection|preparation|"
                           r"acquisition)|specimens?|ethic|data\s+availability|availability\s+of\s+data)", re.I)
SECTYPE_RX = re.compile(r"method|subject|data-availability|data_availability", re.I)
EXCL_RX = re.compile(r"result|discussion|conclusion|supplementar|finding|outcome", re.I)      # top level
EXCL_NESTED_RX = re.compile(r"result|discussion|conclusion|supplementar", re.I)              # inside an included section
ABSTRACT_LABEL = r"(results?|findings|conclusions?|interpretation|discussion)"
ABSTRACT_STOP_RX = re.compile(r"(?:^|\s)" + ABSTRACT_LABEL + r"\s*:", re.I)
TAB_INC = re.compile(r"characteristic|demographic|baseline|participant|patient|subject|clinical", re.I)
TAB_EXC = re.compile(r"metabol|lipid|feature|compound|fold|differential|signific|biomarker|pathway|concentration|"
                     r"abundance|level|biochem|laborator", re.I)
SKIP_IN_TEXT = {"fig", "fig-group", "table-wrap", "table-wrap-group", "table", "supplementary-material",
                "inline-supplementary-material", "boxed-text", "media", "graphic", "inline-graphic", "disp-formula",
                "alternatives", "chem-struct-wrap"}
_ENT_RX = re.compile(r"&(?!(amp|lt|gt|quot|apos|#\d+|#x[0-9a-fA-F]+);)([A-Za-z][A-Za-z0-9]*);")


def _local(tag):
    return tag.rsplit("}", 1)[-1] if isinstance(tag, str) else ""


def _text(el):
    """Text of el without floats (figures, tables, supplements, boxes) nested anywhere inside it."""
    parts = []

    def rec(e):
        if e.text:
            parts.append(e.text)
        for c in e:
            if _local(c.tag) not in SKIP_IN_TEXT:
                rec(c)
            if c.tail:
                parts.append(c.tail)
    rec(el)
    return norm_ws("".join(parts))


def parse_xml(xml_bytes):
    try:
        return ET.fromstring(xml_bytes)
    except ET.ParseError:
        s = xml_bytes.decode("utf-8", "replace")
        s = re.sub(r"<!DOCTYPE[^>]*>", "", s)

        def ent(m):
            name = m.group(2)
            return f"&#{html.entities.name2codepoint[name]};" if name in html.entities.name2codepoint else "&amp;" + name + ";"
        s = _ENT_RX.sub(ent, s)
        s = re.sub(r"&(?![A-Za-z#][A-Za-z0-9]*;)", "&amp;", s)  # bare ampersands ("AT&T")
        return ET.fromstring(s.encode("utf-8"))


def _title(sec):
    t = next((c for c in sec if _local(c.tag) == "title"), None)
    return _text(t) if t is not None else ""


def view_from_xml(xml_bytes, include_abstract):
    root = parse_xml(xml_bytes)
    out = []
    for el in root.iter():
        if _local(el.tag) == "article-title":
            out.append("[TITLE] " + _text(el))
            break
    if include_abstract:
        for el in root.iter():
            if _local(el.tag) != "abstract":
                continue
            secs = [c for c in el if _local(c.tag) == "sec"]
            if secs:  # structured abstract: its Results/Conclusions parts are dropped
                keep = [_text(c) for c in el if _local(c.tag) not in ("sec", "title")]
                keep += [_text(s) for s in secs if not EXCL_RX.search(_title(s))]
            else:     # labelled paragraphs ("<bold>Results:</bold> …") or one paragraph with inline labels
                keep = []
                for c in el:
                    if _local(c.tag) == "title":
                        continue
                    t_ = _text(c)
                    lead = next(iter(c), None)
                    if lead is not None and _local(lead.tag) in ("bold", "italic", "sc", "b", "strong") and not (c.text or "").strip() \
                            and re.fullmatch(ABSTRACT_LABEL + r"\s*[:.]?", _text(lead), re.I):
                        break                      # "<bold>Results</bold> …": this and all later paragraphs dropped
                    m = ABSTRACT_STOP_RX.search(t_)
                    if m and m.start() <= 1:
                        break                      # this and all later paragraphs are results/conclusions
                    keep.append(t_[:m.start()] if m else t_)
                    if m:
                        break
                if not [c for c in el if _local(c.tag) != "title"]:
                    t_ = _text(el)
                    m = ABSTRACT_STOP_RX.search(t_)
                    keep = [t_[:m.start()] if m else t_]
            txt = norm_ws(" ".join(keep))
            if txt:
                out.append("[ABSTRACT] " + txt)

    def walk(el, inherited, path, top):
        for ch in el:
            tag = _local(ch.tag)
            if tag not in ("sec", "notes"):
                continue
            title = _title(ch)
            stype = ch.get("sec-type") or ch.get("notes-type") or ""
            if not title and not stype:
                walk(ch, inherited, path, top)  # untitled, untyped wrapper (e.g. back > notes > sec): transparent
                continue
            excl = EXCL_NESTED_RX if inherited else EXCL_RX
            if excl.search(title) or excl.search(stype):
                continue  # excluded together with all its descendants
            inc = inherited or (top and (bool(METH_START_RX.match(title)) or bool(METH_ANY_RX.search(title))
                                         or bool(SECTYPE_RX.search(stype))))
            p = path + [title or stype or "untitled"]
            if inc:
                paras = [_text(c) for c in ch if _local(c.tag) in ("p", "list", "def-list", "disp-quote")]
                paras = [x for x in paras if x]
                if paras:
                    out.append(f"[SECTION: {' > '.join(p)}] " + " ".join(paras))
            walk(ch, inc, p, False)

    for part in list(root):
        if _local(part.tag) in ("body", "back"):
            walk(part, False, [], True)
    for el in root.iter():
        if _local(el.tag) == "table-wrap":
            lab = next((_text(c) for c in el if _local(c.tag) == "label"), "")
            cap = next((_text(c) for c in el if _local(c.tag) == "caption"), "")
            head = f"{lab} {cap}".strip()
            if TAB_INC.search(head) and not TAB_EXC.search(head):
                rows = []
                for tr in el.iter():
                    if _local(tr.tag) == "tr":
                        rows.append(" | ".join(_text(c) for c in tr if _local(c.tag) in ("td", "th")))
                out.append(f"[TABLE: {head}] " + " // ".join(rows))
    return "\n".join(out) + "\n"


H_INTRO = re.compile(r"(introduction|background)", re.I)
H_METH = re.compile(r"(materials?\s+and\s+methods|methods?\s+and\s+materials?|patients\s+and\s+methods|subjects\s+and\s+methods|"
                    r"research\s+design\s+and\s+methods|(online\s+)?methods?|experimental(\s+section|\s+procedures)?|study\s+design|"
                    r"study\s+population|study\s+participants)", re.I)
H_END = re.compile(r"(results?|discussion|results\s+and\s+discussion|conclusions?|references|acknowledge?ments?|funding|"
                   r"author\s+contributions|supplementary(\s+(materials?|information|data))?|figure\s+legends?|figures?|tables?)", re.I)


def _is_head(line, rx):
    """A heading is a whole short line, optionally numbered, in title or upper case, that IS one of the phrases."""
    s = line.strip()
    if not s or len(s) > 60:
        return False
    s2 = re.sub(r"^(\d+(\.\d+)*\.?|[IVX]+\.)\s*", "", s).rstrip(":").strip()
    if not s2 or not (s2 == s2.upper() or s2[0].isupper()):
        return False
    return rx.fullmatch(s2) is not None


def view_from_pdf_text(text, include_abstract):
    lines = text.splitlines()
    out = ["[VIEW_HEURISTIC_PDF]"]
    i_intro = next((i for i, l in enumerate(lines) if _is_head(l, H_INTRO)), None)
    if include_abstract and i_intro is not None:
        out.append("[FRONT] " + norm_ws(" ".join(lines[:i_intro])))
    start = i_intro + 1 if i_intro is not None else 0
    i_meth = next((i for i in range(start, len(lines)) if _is_head(lines[i], H_METH)), None)
    if i_meth is not None:
        j = next((k for k in range(i_meth + 1, len(lines)) if _is_head(lines[k], H_END)), len(lines))
        out.append("[METHODS_PDF] " + norm_ws(" ".join(lines[i_meth:j])))
    return "\n".join(out) + "\n"


# ----------------------------------------------------------------------------------------------- candidates, order
def load_candidates(path):
    with open(path, "rb") as f:
        b = f.read()
    if sha_bytes(b) != CANDIDATES_SHA:
        raise Refused(f"{path}: SHA-256 {sha_bytes(b)} != frozen {CANDIDATES_SHA}")
    return list(csv.DictReader(io.StringIO(b.decode("utf-8"))))


def order_rows(rows, salt):
    if not SALT_RX.match(salt):
        raise Refused(f"salt does not match the frozen pattern: {salt!r}")
    out = []
    for stratum, origins in (("M", ("M5_SOURCE_UNIQUE", "M5_SOURCE_DEGENERATE")), ("A", ("A35_ALLQC_DIAGNOSTIC",))):
        sel = [r for r in rows if r["origin"] in origins]
        keyed = sorted(((sha_bytes(f"{salt}|{r['accession']}".encode()), r) for r in sel), key=lambda t: t[0])
        for rank, (k, r) in enumerate(keyed, 1):
            out.append({"stratum": stratum, "rank": rank, "accession": r["accession"], "origin": r["origin"],
                        "n_rows": int(r["n_rows"]), "key": k, "in_T1": int(stratum == "M" and rank <= T1_SIZE)})
    return out


ORDER_FIELDS = ["stratum", "rank", "accession", "origin", "n_rows", "key", "in_T1"]


def write_order(order, path):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=ORDER_FIELDS, lineterminator="\n")
    w.writeheader()
    w.writerows(order)
    data = buf.getvalue().encode()
    with open(path, "wb") as f:
        f.write(data)
    return sha_bytes(data)


def read_order(path):
    b = open(path, "rb").read()
    rows = []
    for r in csv.DictReader(io.StringIO(b.decode("utf-8"))):
        rows.append({"stratum": r["stratum"], "rank": int(r["rank"]), "accession": r["accession"], "origin": r["origin"],
                     "n_rows": int(r["n_rows"]), "key": r["key"], "in_T1": int(r["in_T1"])})
    return rows, sha_bytes(b)


def tranche_members(order, tranche):
    if tranche == "T1":
        return [o for o in order if o["in_T1"] == 1]
    return [o for o in order if o["stratum"] == "A"]


# ----------------------------------------------------------------------------------------------- lock / authorisation
def lock_body(salt_option, beacon_not_before, abstract, manual, base="."):
    def h(p):
        q = os.path.join(base, p)
        return sha_file(q) if os.path.exists(q) else None
    return {"docrec_version": VERSION, "protocol": PROTOCOL, "protocol_sha256": h(PROTOCOL),
            "docrec_sha256": sha_file(os.path.abspath(__file__)), "candidates": CANDIDATES,
            "candidates_sha256": h(CANDIDATES), "salt_option": salt_option, "beacon_not_before": beacon_not_before,
            "abstract_in_view": abstract, "manual_texts": manual, "t1_size": T1_SIZE, "s3_max": S3_MAX,
            "owner_signoff": None}


def check_lock(base="."):
    p = os.path.join(base, LOCK)
    if not os.path.exists(p):
        raise Refused("no docrec_lock.json (run plan, then sign)")
    lk = json.load(open(p))
    if not lk.get("owner_signoff"):
        raise Refused("docrec_lock.json is not signed (owner_signoff is null)")
    cur = {"protocol_sha256": sha_file(os.path.join(base, PROTOCOL)), "docrec_sha256": sha_file(os.path.abspath(__file__)),
           "candidates_sha256": sha_file(os.path.join(base, CANDIDATES))}
    for k, v in cur.items():
        if lk.get(k) != v:
            raise Refused(f"lock mismatch for {k}: lock {lk.get(k)} current {v}")
    if lk["candidates_sha256"] != CANDIDATES_SHA:
        raise Refused("lock candidates hash is not the frozen one")
    return lk


def _auth_lines(base="."):
    p = os.path.join(base, AUTH_FILE)
    return [l.strip() for l in open(p, encoding="utf-8")] if os.path.exists(p) else []


def check_fetch_auth(lk, tranche, base="."):
    """The LAST fetch line for the tranche governs (a corrected line may be appended)."""
    hits = [m for m in (FETCH_AUTH_RX.match(l) for l in _auth_lines(base)) if m and m.group(2) == tranche]
    if not hits:
        raise Refused(f"no 'DOCREC FETCH AUTHORISED … tranche={tranche} …' line in {AUTH_FILE}")
    m = hits[-1]
    if m.group(3) != CANDIDATES_SHA or m.group(4) != PROTOCOL:
        raise Refused("authorisation line names another candidates file or protocol")
    if m.group(5) != sha_file(os.path.join(base, LOCK)):
        raise Refused("authorisation line names another lock (lock= must be the SHA-256 of the signed docrec_lock.json)")
    salt = m.group(6)
    sm = SALT_RX.match(salt)
    if not sm:
        raise Refused("authorisation salt does not match the frozen pattern")
    if lk["salt_option"] == "a" and salt != "stageA-docrec-v1":
        raise Refused("lock says salt option (a) but the salt is not the fixed literal")
    if lk["salt_option"] == "b" and not sm.group(2):
        raise Refused(f"lock says salt option (b) but the salt is not NISTv2:<chain>-<pulse>:<outputValue> from {lk.get('beacon_url')}")
    return salt, m.group(7)


def seal_signed(tranche, base="."):
    return any(m and m.group(2) == tranche for m in (SEAL_AUTH_RX.match(l) for l in _auth_lines(base)))


def check_seal_auth(tranche, universe_sha, base="."):
    hits = [m for m in (SEAL_AUTH_RX.match(l) for l in _auth_lines(base)) if m and m.group(2) == tranche]
    if not hits or hits[-1].group(3) != universe_sha:
        raise Refused(f"the last 'DOCREC SEAL SIGNED … tranche={tranche}' line must name universe={universe_sha[:12]}…")
    return hits[-1].group(4)


# ----------------------------------------------------------------------------------------------- fetch
def fetch_study(t, o, store, include_abstract, tranche, attempt=1):
    acc = o["accession"]
    d = os.path.join(store, acc)
    os.makedirs(d, exist_ok=True)
    man = {"accession": acc, "tranche": tranche, "n_rows": o["n_rows"], "docrec_version": VERSION, "files": {},
           "s2": [], "s3": [], "notes": [], "indeterminate": [], "attempt": attempt, "complete": False}
    allow_sa = tranche == "A35"

    def save(name, data, role, **extra):
        with open(os.path.join(d, name), "wb") as f:
            f.write(data)
        man["files"][name] = {"sha256": sha_bytes(data), "role": role, **extra}

    def get(url, label):
        ok, st, body, note, definite = t.get(url, allow_sa)
        if not ok:
            man["notes"].append(f"{label} {st} {note}")
            if not definite:
                man["indeterminate"].append(label)
        return ok, body

    def get_json(url, label):
        ok, body = get(url, label)
        if not ok:
            return None, None
        try:
            return json.loads(body), body
        except ValueError as e:
            man["notes"].append(f"{label} BAD_JSON {one_line(e, 80)}")
            man["indeterminate"].append(label)  # possibly truncated without Content-Length
            return None, None

    def finish():
        if man["indeterminate"] and attempt < MAX_ATTEMPTS:
            man["complete"] = False
        else:
            if man["indeterminate"]:
                man["notes"].append(f"FINALISED_AFTER_{attempt}_ATTEMPTS")
            man["complete"] = True
        with open(os.path.join(d, "manifest.json"), "w") as f:
            json.dump(man, f, indent=1, sort_keys=True)
        return man

    ok, body = get(mirror_url(acc), "S1_INDEX")
    if not ok:
        return finish()
    save("S1_index.html", body, "S1_index")
    names = index_names(body)
    inames = sorted(n for n in names if re.fullmatch(r"i_[^/]+\.txt", n))
    if not inames:
        man["notes"].append("S1_NO_INVESTIGATION_FILE")
    inv = {"pmids": [], "dois": [], "study_files": [], "assay_files": []}
    for n in inames:
        ok, body = get(mirror_url(acc, n), f"S1_FETCH {n}")
        if not ok:
            continue
        save("S1_" + n, body, "S1", usable=True)
        p = parse_investigation(body.decode("utf-8", "replace"))
        for k in inv:
            inv[k] = list(dict.fromkeys(inv[k] + p[k]))
    man["investigation"] = inv
    if allow_sa:
        for n in inv["study_files"] + inv["assay_files"]:
            if not re.fullmatch(r"[sa]_[^/]+\.txt", n) or n not in names:
                man["notes"].append(f"S1_SA_NOT_LISTED {n}")
                continue
            ok, body = get(mirror_url(acc, n), f"S1_FETCH {n}")
            if ok:
                save("S1_" + n, body, "S1", usable=True)
    # S2: publications named in the investigation file
    s2, s2_ids = [], set()
    for kind, val in [("pmid", p_) for p_ in inv["pmids"]] + [("doi", d_) for d_ in inv["dois"]]:
        q = f"EXT_ID:{val} AND SRC:MED" if kind == "pmid" else f'DOI:"{val}"'
        js, body = get_json(search_url(q, 25), f"S2_SEARCH {kind}:{val}")
        if js is None:
            continue
        save(f"S2_search_{kind}_{sha_bytes(val.encode())[:10]}.json", body, "S2_search", query=q)
        res = (js.get("resultList") or {}).get("result") or [] if isinstance(js, dict) else []
        for r in sorted((r for r in res if isinstance(r, dict)), key=_sort_key):
            match = (str(r.get("pmid")) == val) if kind == "pmid" else (str(r.get("doi") or "").lower() == val)
            if match and not (art_ids(r) & s2_ids):
                s2_ids |= art_ids(r)
                s2.append(art_slim(r))
    # S3: accession search (DR-2 frozen parameters); not run for stratum A
    s3_all, cursor, page = [], "*", 0
    if tranche == "T1":
        while len(s3_all) < S3_CAP:
            js, body = get_json(search_url(f'"{acc}"', 1000, cursor), f"S3_SEARCH p{page + 1}")
            if js is None:
                break
            page += 1
            save(f"S3_search_p{page}.json", body, "S3_search", query=f'"{acc}"')
            res = (js.get("resultList") or {}).get("result") or [] if isinstance(js, dict) else []
            s3_all += [r for r in res if isinstance(r, dict)]
            nxt = js.get("nextCursorMark") if isinstance(js, dict) else None
            if not res or not nxt or nxt == cursor:
                break
            cursor = nxt
        man["s3_hits"] = len(s3_all)
    s3 = [art_slim(r) for r in s3_select(s3_all, s2_ids)]
    # full texts and views
    for role, arts in (("S2", s2), ("S3", s3)):
        for a in arts:
            pmcid = str(a.get("pmcid") or "")
            if not re.fullmatch(r"PMC[1-9]\d*", pmcid):
                a["fulltext"] = "NO_PMCID"
                continue
            ok, body = get(fulltext_url(pmcid), f"{role}_FULLTEXT {pmcid}")
            if not ok:
                a["fulltext"] = "UNAVAILABLE"
                continue
            save(f"{role}_{pmcid}.xml", body, f"{role}_fulltext")
            try:
                usable = True
                if role == "S3":
                    usable = accession_in_text(acc, norm_ws("".join(parse_xml(body).itertext())))
                    a["accession_in_text"] = usable
                v = view_from_xml(body, include_abstract).encode("utf-8")
                save(f"{role}_{pmcid}.view.txt", v, f"{role}_view", usable=usable)
                a["fulltext"] = "OK"
            except (ET.ParseError, ValueError) as e:
                a["fulltext"] = "VIEW_PARSE_ERROR"
                man["notes"].append(f"{role}_VIEW_PARSE_ERROR {pmcid} {one_line(e, 80)}")
                man["indeterminate"].append(f"{role}_VIEW {pmcid}")  # possibly truncated without Content-Length
    man["s2"], man["s3"] = s2, s3
    return finish()


def cmd_fetch(args):
    lk = check_lock()
    if seal_signed(args.tranche):
        raise Refused(f"tranche {args.tranche} is already sealed and signed; no further retrieval")
    salt, owner = check_fetch_auth(lk, args.tranche)
    order = order_rows(load_candidates(CANDIDATES), salt)
    os.makedirs(OUTDIR, exist_ok=True)
    order_sha = write_order(order, order_path(args.tranche))
    members = tranche_members(order, args.tranche)
    store = os.path.join(OUTDIR, "store")
    os.makedirs(store, exist_ok=True)
    t = Transport(os.path.join(OUTDIR, "docrec_fetch_log.tsv"))
    t.header(f"run {now()} docrec {VERSION} tranche {args.tranche} n={len(members)} order_sha256 {order_sha} salt {salt} owner {owner}")
    code, touched = 0, []
    try:
        for i, o in enumerate(members, 1):
            mp = os.path.join(store, o["accession"], "manifest.json")
            attempt = 1
            if os.path.exists(mp):
                prev = json.load(open(mp))
                if prev.get("complete"):
                    print(f"[{i}/{len(members)}] {o['accession']} complete — skipped (resume)")
                    continue
                attempt = int(prev.get("attempt", 0)) + 1
            touched.append(o["accession"])
            man = fetch_study(t, o, store, lk["abstract_in_view"], args.tranche, attempt)
            print(f"[{i}/{len(members)}] {o['accession']} files={len(man['files'])} s2={len(man['s2'])} s3={len(man['s3'])} "
                  f"complete={man['complete']} notes={len(man['notes'])}", flush=True)
    except Outage as e:
        # a run that ends in an outage stop does not count as an attempt for the studies it touched
        for acc in touched:
            mp = os.path.join(store, acc, "manifest.json")
            if os.path.exists(mp):
                m = json.load(open(mp))
                if m.get("indeterminate"):
                    m["complete"] = False
                    m["attempt"] = max(0, int(m.get("attempt", 1)) - 1)
                    m["notes"] = [x for x in m["notes"] if not x.startswith("FINALISED_AFTER_")] + ["RUN_ENDED_IN_OUTAGE_STOP"]
                    with open(mp, "w") as f:
                        json.dump(m, f, indent=1, sort_keys=True)
        print(f"OUTAGE STOP: {e}. Rerun later to resume.", file=sys.stderr)
        code = 2
    except Refused as e:
        print(f"REFUSED (protocol violation, run halted): {e}", file=sys.stderr)
        code = 1
    finally:
        t.close()
    if code == 0:
        inc = [o["accession"] for o in members
               if not (os.path.exists(os.path.join(store, o["accession"], "manifest.json"))
                       and json.load(open(os.path.join(store, o["accession"], "manifest.json"))).get("complete"))]
        if inc:
            print(f"{len(inc)} studies incomplete (transient failures): rerun fetch. {inc[:10]}")
            code = 3
        else:
            print(f"all {len(members)} studies complete. order_sha256 {order_sha}")
    return code


# ----------------------------------------------------------------------------------------------- register / seal
def cmd_register(args):
    lk = check_lock()
    if not lk["manual_texts"]:
        raise Refused("DR-9 = no in the signed lock: owner-supplied texts are not admissible")
    if seal_signed(args.tranche):
        raise Refused(f"tranche {args.tranche} is already sealed and signed")
    store = os.path.join(OUTDIR, "store")
    members, _ = _tranche_from_order(args.tranche)
    inc = [o["accession"] for o in members if not (os.path.exists(os.path.join(store, o["accession"], "manifest.json"))
                                                   and json.load(open(os.path.join(store, o["accession"], "manifest.json"))).get("complete"))]
    if inc:
        raise Refused(f"register only after every tranche member is complete (a later fetch would rewrite the manifest): {inc[:10]}")
    n = 0
    for acc in sorted(os.listdir(store)):
        md = os.path.join(store, acc, "manual")
        man_p = os.path.join(store, acc, "manifest.json")
        if not os.path.isdir(md) or not os.path.exists(man_p):
            continue
        man = json.load(open(man_p))
        if man.get("tranche") != args.tranche:
            continue
        for f in sorted(os.listdir(md)):
            if not f.lower().endswith(".pdf"):
                continue
            stem = f[:-4]
            tp = os.path.join(md, stem + ".txt")
            if not os.path.exists(tp):
                raise Refused(f"{acc}/manual/{f}: text rendering {stem}.txt missing (pdftotext)")
            view = view_from_pdf_text(open(tp, encoding="utf-8", errors="replace").read(), lk["abstract_in_view"]).encode("utf-8")
            vname = "S2M_" + re.sub(r"[^A-Za-z0-9._-]", "_", stem) + ".view.txt"
            with open(os.path.join(store, acc, vname), "wb") as fh:
                fh.write(view)
            man["files"][f"manual/{f}"] = {"sha256": sha_file(os.path.join(md, f)), "role": "S2M_pdf"}
            man["files"][f"manual/{stem}.txt"] = {"sha256": sha_file(tp), "role": "S2M_text"}
            man["files"][vname] = {"sha256": sha_bytes(view), "role": "S2M_view", "usable": True}
            n += 1
        with open(man_p, "w") as fh:
            json.dump(man, fh, indent=1, sort_keys=True)
    print(f"registered {n} owner-supplied text(s)")
    return 0


def last_order_sha(tranche):
    p = os.path.join(OUTDIR, "docrec_fetch_log.tsv")
    shas = []
    for l in (open(p, encoding="utf-8") if os.path.exists(p) else []):
        m = re.match(r"^# run \S+ docrec \S+ tranche (\S+) n=\d+ order_sha256 ([0-9a-f]{64}) ", l)
        if m and m.group(1) == tranche:
            shas.append(m.group(2))
    return shas[-1] if shas else None


def universe(tranche, store, members, order_sha):
    out = {"tranche": tranche, "docrec_version": VERSION, "order_sha256": order_sha, "studies": {}}
    problems = []
    for o in members:
        acc = o["accession"]
        mp = os.path.join(store, acc, "manifest.json")
        if not os.path.exists(mp):
            problems.append(f"{acc}: no manifest")
            continue
        man = json.load(open(mp))
        if man.get("tranche") != tranche or not man.get("complete"):
            problems.append(f"{acc}: incomplete or wrong tranche")
            continue
        files = {}
        for name, meta in man["files"].items():
            p = os.path.join(store, acc, name)
            if not os.path.exists(p) or sha_file(p) != meta["sha256"]:
                raise Refused(f"{acc}/{name}: missing or changed since retrieval")
            files[name] = {"sha256": meta["sha256"], "role": meta["role"], "usable": bool(meta.get("usable", False))}
        out["studies"][acc] = {"n_rows": o["n_rows"], "files": files, "notes": man["notes"]}
    if problems:
        raise Refused("the tranche is not complete: " + "; ".join(problems[:10]) + (" …" if len(problems) > 10 else ""))
    return out


def order_path(tranche):
    return os.path.join(OUTDIR, f"docrec_{tranche}_order.csv")


def _tranche_from_order(tranche):
    order, osha = read_order(order_path(tranche))
    return tranche_members(order, tranche), osha


def cmd_seal(args):
    check_lock()
    if seal_signed(args.tranche):
        raise Refused(f"tranche {args.tranche} is already sealed and signed; it cannot be resealed")
    members, osha = _tranche_from_order(args.tranche)
    if osha != last_order_sha(args.tranche):
        raise Refused(f"{order_path(args.tranche)} does not match the order recorded in the last fetch run for this tranche")
    expected = T1_SIZE if args.tranche == "T1" else sum(1 for _ in members)
    if len(members) != expected:
        raise Refused(f"tranche has {len(members)} members, expected {expected}")
    u = universe(args.tranche, os.path.join(OUTDIR, "store"), members, osha)
    data = json.dumps(u, indent=1, sort_keys=True).encode()
    p = os.path.join(OUTDIR, f"docrec_{args.tranche}_evidence_universe.json")
    with open(p, "wb") as f:
        f.write(data)
    print(f"{p}: {len(u['studies'])} studies, SHA-256 {sha_bytes(data)}")
    print(f"sign with: echo 'DOCREC SEAL SIGNED {now()[:10]} tranche={args.tranche} universe={sha_bytes(data)} owner=<name>' >> {AUTH_FILE}")
    return 0


# ----------------------------------------------------------------------------------------------- record validation
Q_ENUM = {"Q1": {"yes", "no", "not_stated"}, "Q2": {"established", "not_established"},
          "Q4": {"explicit_one", "explicit_repeated", "not_stated"}, "Q7": {"one_to_one", "shows_repeated", "none"}}
NEEDS_EVIDENCE = {"Q1": {"yes", "no"}, "Q2": {"established"}, "Q4": {"explicit_one", "explicit_repeated"},
                  "Q7": {"one_to_one", "shows_repeated"}}
V_ENUM = {"pass", "fail", "not_established"}


def _is_pub(src):
    return src.startswith(("S2_", "S3_", "S2M_"))


def _evidence_list_ok(ev):
    return isinstance(ev, list) and all(isinstance(e, dict) and all(isinstance(e.get(k), str) for k in ("source", "locator", "quote"))
                                        for e in ev)


def _int_ok(c):
    return isinstance(c, int) and not isinstance(c, bool) and c >= 1


def validate_common(rec, tranche):
    if not isinstance(rec, dict):
        return ["record is not a JSON object"]
    probs = []
    if not isinstance(rec.get("accession"), str) or not re.fullmatch(ACC, rec["accession"]):
        probs.append("accession")
    if rec.get("tranche") != tranche:
        probs.append(f"tranche field must be {tranche!r}")
    if not isinstance(rec.get("extractor"), str) or not rec["extractor"].strip():
        probs.append("extractor")
    if not isinstance(rec.get("extracted_utc"), str) or not UTC_RX.match(rec["extracted_utc"]):
        probs.append("extracted_utc (YYYY-MM-DDTHH:MM:SSZ)")
    sc = rec.get("sources_consulted")
    if not isinstance(sc, list) or not all(isinstance(x, str) for x in sc):
        probs.append("sources_consulted (list of file names)")
    return probs


def verify_evidence(items, ustudy, views):
    """Each quote: ≥ 3 tokens, no view markers, verbatim (whitespace-normalised, case-sensitive) in a usable view or
    S1 file of this study."""
    probs = []
    for it in items:
        src, quote, loc = it["source"], norm_ws(it["quote"]), it["locator"]
        meta = ustudy["files"].get(src)
        if not meta or not (meta["usable"] and (meta["role"].endswith("_view") or meta["role"] == "S1")):
            probs.append(f"unusable or unknown source {src!r}")
            continue
        if not loc.strip():
            probs.append(f"no locator for a quote from {src}")
        if len(quote.split()) < MIN_QUOTE_TOKENS:
            probs.append(f"quote shorter than {MIN_QUOTE_TOKENS} tokens: {quote!r}")
        elif any(m in quote for m in VIEW_MARKERS):
            probs.append(f"quote spans a view marker: {quote[:60]!r}")
        elif quote not in views[src]:
            probs.append(f"quote not found verbatim in {src}: {quote[:60]!r}")
    return probs


UNIT_RX = (r"(µl|μl|ul|ml|l\b|mg|µg|μg|ng|pg|g\b|kg|nm|µm|μm|mm|cm|m\b|mm|µm|nmol|µmol|mmol|mol\b|mm\b|°|%|rpm|min\b|"
           r"h\b|hr|hours?|s\b|sec|seconds?|days?|weeks?|months?|years?|bp\b|kb\b|da\b|kda|ppm|x\s*g\b|×)")
_NUM_WORDS = {w: i for i, w in enumerate("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen "
                                         "fifteen sixteen seventeen eighteen nineteen".split())}
_NUM_WORDS.update({w: 10 * (i + 2) for i, w in enumerate("twenty thirty forty fifty sixty seventy eighty ninety".split())})


def words_to_int(txt):
    """English number words up to 9999 ("One hundred and twenty-two"); None if any token is not a number word."""
    total, cur = 0, 0
    toks = [t for t in re.split(r"[\s-]+", txt.lower()) if t and t != "and"]
    if not toks:
        return None
    for t in toks:
        if t in _NUM_WORDS:
            cur += _NUM_WORDS[t]
        elif t == "hundred":
            cur = (cur or 1) * 100
        elif t == "thousand":
            total += (cur or 1) * 1000
            cur = 0
        else:
            return None
    return total + cur


def check_count(q, count):
    """The count must appear as an integer in a quote of q, or be the sum of a count_derivation whose items each
    appear verbatim in a quote. Returns (problems, flags)."""
    quotes = [norm_ws(e["quote"]) for e in q.get("evidence", [])]
    forms = {str(count), f"{count:,}"}
    if any(re.search(r"(?<![\d,.])" + re.escape(f) + r"(?![\d]|[.,]\d)(?!\s*" + UNIT_RX + r")", qq, re.I) for f in forms for qq in quotes):
        return [], []
    der = q.get("count_derivation")
    if not isinstance(der, list) or not der or not all(isinstance(d, dict) and isinstance(d.get("text"), str) and _int_ok(d.get("value"))
                                                     for d in der):
        return [f"count {count} does not appear in the quotes and no valid count_derivation is given"], []
    probs, flags = [], []
    for d in der:
        txt = norm_ws(d["text"])
        if not txt or not any(txt in qq for qq in quotes):
            probs.append(f"count_derivation text not in a quote: {txt!r}")
        elif re.fullmatch(r"\d{1,3}(,\d{3})*|\d+", txt):
            if int(txt.replace(",", "")) != d["value"]:
                probs.append(f"count_derivation {txt!r} != {d['value']}")
        elif re.search(r"\d", txt):
            probs.append(f"count_derivation text mixes digits and words: {txt!r}")
        elif words_to_int(txt) != d["value"]:
            probs.append(f"count_derivation {txt!r} does not read as {d['value']}")
        else:
            flags.append("NUMBER_WORDS")
    if sum(d["value"] for d in der) != count:
        probs.append(f"count_derivation sums to {sum(d['value'] for d in der)}, not {count}")
    return probs, sorted(set(flags))


def grade_record(rec, n_rows, ustudy, views, tranche="T1"):
    """Deterministic grade per protocol §9. Returns (grade, reason, problems, flags)."""
    probs = validate_common(rec, tranche)
    if isinstance(rec, dict):
        for q in ("Q1", "Q2", "Q4", "Q7"):
            v = rec.get(q)
            if not isinstance(v, dict) or v.get("answer") not in Q_ENUM[q] or not _evidence_list_ok(v.get("evidence", [])):
                probs.append(f"{q}: must be an object with a valid answer and an evidence list")
            elif v["answer"] in NEEDS_EVIDENCE[q] and not v.get("evidence"):
                probs.append(f"{q}: answer {v['answer']} needs evidence")
        q3 = rec.get("Q3")
        if not isinstance(q3, dict) or q3.get("scope") not in ("selected_assay", "study_enrollment", "none") \
                or not _evidence_list_ok(q3.get("evidence", [])):
            probs.append("Q3: must be an object with scope and an evidence list")
        elif q3["scope"] != "none":
            if not _int_ok(q3.get("count")):
                probs.append("Q3: count must be a positive integer")
            if not q3.get("evidence"):
                probs.append("Q3: count needs evidence")
        q6 = rec.get("Q6")
        if not isinstance(q6, dict) or not isinstance(q6.get("flags", []), list) or not all(isinstance(x, str) for x in q6.get("flags", [])) \
                or not isinstance(q6.get("confirmed_overlap", []), list) \
                or not all(isinstance(x, str) and re.fullmatch(ACC, x) and x != rec.get("accession") for x in q6.get("confirmed_overlap", [])) \
                or not _evidence_list_ok(q6.get("evidence", [])):
            probs.append("Q6: flags (list), confirmed_overlap (list of other accessions) and an evidence list")
        elif q6.get("confirmed_overlap") and not q6.get("evidence"):
            probs.append("Q6: confirmed overlap needs evidence")
    if probs:
        return "INVALID", "SCHEMA", probs, []
    flags = []
    for q in ("Q1", "Q2", "Q3", "Q4", "Q6", "Q7"):
        probs += [f"{q}: {p}" for p in verify_evidence(rec[q].get("evidence", []), ustudy, views)]
    if rec["Q3"]["scope"] != "none" and not probs:
        cp, cf = check_count(rec["Q3"], rec["Q3"]["count"])
        probs += [f"Q3: {p}" for p in cp]
        flags += cf
    if probs:
        return "INVALID", "QUOTE", probs, []
    if any(e["source"].startswith("S2M_") for q in ("Q1", "Q2", "Q3", "Q4", "Q6", "Q7") for e in rec[q].get("evidence", [])):
        flags.append("VIEW_HEURISTIC_PDF")
    acc = rec["accession"]
    q1, q2, q3, q4, q7 = rec["Q1"]["answer"], rec["Q2"]["answer"], rec["Q3"], rec["Q4"]["answer"], rec["Q7"]["answer"]

    def ev(q):
        return rec[q].get("evidence", [])

    def res(g, r):
        return g, r, [], sorted(set(flags))
    if q2 != "established":
        return res("DOC_UNRESOLVED", "Q2_NOT_ESTABLISHED")                       # rule 1
    if q4 == "explicit_repeated" or q7 == "shows_repeated":
        return res("DOC_CONTRADICTED", "REPEATED_SAMPLING")                      # rule 2
    if q1 == "no":
        return res("DOC_CONTRADICTED", "NOT_PARTICIPANT_SAMPLES")                # rule 3
    if q1 != "yes":
        return res("DOC_UNRESOLVED", "Q1_NOT_STATED")                            # rule 4
    if q3["scope"] != "selected_assay":
        return res("DOC_UNRESOLVED", "NO_ASSAY_SCOPED_COUNT")                    # rule 5
    if q3["count"] > n_rows:
        return res("DOC_UNRESOLVED", "COUNT_EXCEEDS_ROWS")                       # rule 6
    if q3["count"] < n_rows:
        return res("DOC_UNRESOLVED", "ROWS_EXCEED_COUNT")                        # rule 7
    q2_pub = [e for e in ev("Q2") if _is_pub(e["source"])]
    q3_pub = [e for e in ev("Q3") if _is_pub(e["source"])]
    if not q2_pub or not q3_pub:
        return res("DOC_UNRESOLVED", "S1_ONLY_EVIDENCE")                         # rule 8
    if all(e["source"].startswith("S3_") for e in q2_pub) and not any(accession_in_text(acc, e["quote"]) for e in q2_pub):
        return res("DOC_UNRESOLVED", "Q2_S3_WITHOUT_ACCESSION")                  # rule 8 (S3 restriction)
    if (q4 == "explicit_one" and any(_is_pub(e["source"]) for e in ev("Q4"))) or \
            (q7 == "one_to_one" and any(_is_pub(e["source"]) for e in ev("Q7"))):
        return res("DOC_RECOVERED_STRICT", "PROVISIONAL")                        # rule 9
    return res("DOC_RECOVERED_COUNT", "SENSITIVITY_ONLY")                        # rule 10


def grade_v_record(rec, ustudy, views, tranche="A35"):
    """Returns (status, problems, v3_count, flags)."""
    probs = validate_common(rec, tranche)
    if isinstance(rec, dict):
        for v in ("V1", "V2", "V3", "V4"):
            x = rec.get(v)
            if not isinstance(x, dict) or x.get("answer") not in V_ENUM or not _evidence_list_ok(x.get("evidence", [])):
                probs.append(f"{v}: must be an object with a valid answer and an evidence list")
            elif x["answer"] in ("pass", "fail") and not x.get("evidence"):
                probs.append(f"{v}: {x['answer']} needs evidence")
        if isinstance(rec.get("V3"), dict) and rec["V3"].get("answer") == "pass" and not _int_ok(rec["V3"].get("count")):
            probs.append("V3: pass needs a positive integer count")
    if probs:
        return "INVALID", probs, None, []
    flags = []
    for v in ("V1", "V2", "V3", "V4"):
        probs += [f"{v}: {p}" for p in verify_evidence(rec[v].get("evidence", []), ustudy, views)]
    if rec["V3"]["answer"] == "pass" and not probs:
        cp, cf = check_count(rec["V3"], rec["V3"]["count"])
        probs += [f"V3: {p}" for p in cp]
        flags += cf
    if probs:
        return "INVALID", probs, None, []
    ans = [rec[v]["answer"] for v in ("V1", "V2", "V3", "V4")]
    if "fail" in ans:
        return "FAILED", [], None, flags
    if all(a == "pass" for a in ans):
        return "VERIFIED", [], rec["V3"]["count"], flags
    return "UNVERIFIABLE", [], None, flags


def record_sha(rec):
    return sha_bytes(json.dumps(rec, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8"))


# ----------------------------------------------------------------------------------------------- bounds
def binom_cdf(x, n, p):
    return sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(0, x + 1))


def cp_upper(x, n, a=ALPHA):
    if x >= n:
        return 1.0
    lo, hi = x / n, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if binom_cdf(x, n, mid) > a else (lo, mid)
    return (lo + hi) / 2


def cp_lower(x, n, a=ALPHA):
    if x <= 0:
        return 0.0
    lo, hi = 0.0, x / n
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if 1 - binom_cdf(x - 1, n, mid) < a else (lo, mid)
    return (lo + hi) / 2


def hyper_cdf(x, N, K, n):
    tot = math.comb(N, n)
    return sum(math.comb(K, k) * math.comb(N - K, n - k) for k in range(max(0, n - (N - K)), min(x, K) + 1)) / tot


def hyper_upper(x, N, n, a=ALPHA):
    best = x
    for K in range(x, N - (n - x) + 1):
        if hyper_cdf(x, N, K, n) > a:
            best = K
    return best


def hyper_lower(x, N, n, a=ALPHA):
    if x <= 0:
        return 0
    best = N - (n - x)
    for K in range(N - (n - x), x - 1, -1):
        if 1 - hyper_cdf(x - 1, N, K, n) > a:
            best = K
    return best


def futility(x, n=T1_SIZE, N=POOL_M, salt_option="b"):
    up_cp, lo_cp = cp_upper(x, n) * N, cp_lower(x, n) * N
    up_h, lo_h = hyper_upper(x, N, n), hyper_lower(x, N, n)
    upper = up_h if salt_option == "b" else up_cp
    return {"x": x, "n": n, "N": N, "cp_upper_pool": round(up_cp, 1), "cp_lower_pool": round(lo_cp, 1),
            "hypergeometric_upper_K": up_h, "hypergeometric_lower_K": lo_h,
            "bound_used": "hypergeometric (design-based, simple random sample)" if salt_option == "b"
            else "Clopper-Pearson (model-based; systematic hash-distributed sample)",
            "futility_pool": FUTILITY_POOL, "DR7_futility_reached": upper < FUTILITY_POOL}


def components(accs, edges):
    parent = {a: a for a in accs}

    def f(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for a, b in edges:
        if a in parent and b in parent:
            parent[f(a)] = f(b)
    return len({f(a) for a in accs})


# ----------------------------------------------------------------------------------------------- grade command
def load_adjudication(path):
    """accession, provisional_grade, record_sha256, confirmed (yes/no), note — UTF-8 with or without BOM."""
    adj = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path, encoding="utf-8-sig")):
            r = {(k or "").strip(): (v or "").strip() for k, v in r.items()}
            adj[r.get("accession", "")] = r
    return adj


def cmd_grade(args):
    lk = check_lock()
    up = os.path.join(OUTDIR, f"docrec_{args.tranche}_evidence_universe.json")
    u_bytes = open(up, "rb").read()
    check_seal_auth(args.tranche, sha_bytes(u_bytes))
    u = json.loads(u_bytes)
    members, osha = _tranche_from_order(args.tranche)
    if osha != u["order_sha256"]:
        raise Refused(f"{order_path(args.tranche)} differs from the order sealed in the evidence universe")
    cur = universe(args.tranche, os.path.join(OUTDIR, "store"), members, osha)
    if json.dumps(cur, sort_keys=True) != json.dumps(u, sort_keys=True):
        raise Refused("store differs from the sealed evidence universe")
    ev_p = os.path.join(OUTDIR, f"docrec_{args.tranche}_evidence.jsonl")
    recs = []
    for i, l in enumerate(open(ev_p, encoding="utf-8"), 1):
        if l.strip():
            try:
                recs.append(json.loads(l))
            except ValueError as e:
                raise Refused(f"{ev_p} line {i}: not JSON ({one_line(e, 80)})")
    adj = load_adjudication(os.path.join(OUTDIR, f"docrec_{args.tranche}_adjudication.csv"))
    store = os.path.join(OUTDIR, "store")
    rows, edges, seen, used_adj = [], [], set(), set()
    for rec in recs:
        acc = rec.get("accession") if isinstance(rec, dict) else None
        if acc not in u["studies"] or acc in seen:
            raise Refused(f"evidence record for {acc!r} is outside the sealed tranche or duplicated")
        seen.add(acc)
        us = u["studies"][acc]
        views = {n: norm_ws(open(os.path.join(store, acc, n), encoding="utf-8", errors="replace").read())
                 for n, m in us["files"].items() if m["usable"]}
        rsha = record_sha(rec)
        if args.tranche == "A35":
            g, probs, cnt, flags = grade_v_record(rec, us, views, "A35")
            rows.append({"accession": acc, "status": g, "v3_count": cnt, "record_sha256": rsha, "flags": "|".join(flags),
                         "problems": " ; ".join(probs)})
            continue
        g, reason, probs, flags = grade_record(rec, us["n_rows"], us, views, "T1")
        final, freason = g, reason
        if g in ("DOC_RECOVERED_STRICT", "DOC_CONTRADICTED"):
            a = adj.get(acc)
            if a and a.get("provisional_grade") == g and a.get("record_sha256") == rsha and a.get("confirmed", "").lower() in ("yes", "no"):
                used_adj.add(acc)
                if a["confirmed"].lower() == "yes":
                    freason = "ADJUDICATED"
                else:
                    final, freason = "DOC_UNRESOLVED", "ADJUDICATION_REJECTED"
            else:
                freason = "PROVISIONAL"
        for o in rec["Q6"].get("confirmed_overlap", []) if g != "INVALID" else []:
            edges.append((acc, o))
        rows.append({"accession": acc, "provisional_grade": g, "provisional_reason": reason, "grade": final,
                     "reason": freason, "record_sha256": rsha, "flags": "|".join(flags),
                     "q6_flags": "|".join(rec["Q6"].get("flags", [])) if g != "INVALID" else "",
                     "problems": " ; ".join(probs)})
    missing = sorted(set(u["studies"]) - seen)
    out = os.path.join(OUTDIR, f"docrec_{args.tranche}_results.csv")
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["accession"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    report = {"tranche": args.tranche, "docrec_version": VERSION, "n_sealed": len(u["studies"]), "n_records": len(recs),
              "missing_records": missing, "stale_adjudication_rows": sorted(set(adj) - used_adj - {""})}
    if args.tranche == "A35":
        report["N_A35_verified"] = sum(r["status"] == "VERIFIED" for r in rows)
        report["status_counts"] = {s: sum(r["status"] == s for r in rows) for s in ("VERIFIED", "FAILED", "UNVERIFIABLE", "INVALID")}
    else:
        report["grade_counts"] = {}
        for r in rows:
            k = f"{r['grade']} ({r['reason']})"
            report["grade_counts"][k] = report["grade_counts"].get(k, 0) + 1
        pending = [r["accession"] for r in rows if r["reason"] == "PROVISIONAL"]
        invalid = [r["accession"] for r in rows if r["grade"] == "INVALID"]
        strict = [r["accession"] for r in rows if r["grade"] == "DOC_RECOVERED_STRICT" and r["reason"] == "ADJUDICATED"]
        count = [r["accession"] for r in rows if r["grade"] == "DOC_RECOVERED_COUNT"]
        report["pending_adjudication"], report["invalid"] = pending, invalid
        report["strict_overlap_collapsed"] = components(strict, edges)
        report["strict_plus_count_overlap_collapsed (sensitivity)"] = components(strict + count, edges)
        if pending or invalid or missing:
            report["futility"] = "NOT COMPUTED: records missing, invalid or awaiting adjudication"
        else:
            report["futility"] = futility(report["strict_overlap_collapsed"], salt_option=lk["salt_option"])
    rp = os.path.join(OUTDIR, f"docrec_{args.tranche}_report.json")
    with open(rp, "w") as f:
        json.dump(report, f, indent=1, sort_keys=True)
    print(json.dumps(report, indent=1, sort_keys=True))
    return 0


# ----------------------------------------------------------------------------------------------- plan / order
def beacon_url(iso):
    ms = int(datetime.strptime(iso, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc).timestamp() * 1000)
    return f"{BEACON}{ms}"


def cmd_plan(args):
    if args.salt_option == "b":
        if not args.beacon_not_before:
            raise Refused("salt option (b) needs --beacon-not-before YYYY-MM-DDTHH:MM:SSZ (a time after signing)")
        t = datetime.strptime(args.beacon_not_before, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
        if t <= datetime.now(timezone.utc):
            raise Refused("--beacon-not-before must be in the future: the pulse must not exist when the lock is signed")
    body = lock_body(args.salt_option, args.beacon_not_before if args.salt_option == "b" else None,
                     args.abstract == "yes", args.manual == "yes")
    body["beacon_url"] = beacon_url(args.beacon_not_before) if args.salt_option == "b" else None
    missing = [k for k in ("protocol_sha256", "candidates_sha256") if body[k] is None]
    if missing:
        raise Refused(f"files missing for the lock: {missing}")
    if body["candidates_sha256"] != CANDIDATES_SHA:
        raise Refused("candidates file is not the frozen one")
    with open(LOCK, "w") as f:
        json.dump(body, f, indent=1, sort_keys=True)
    print(json.dumps(body, indent=1, sort_keys=True))
    return 0


def cmd_order(args):
    order = order_rows(load_candidates(CANDIDATES), args.salt)
    for o in tranche_members(order, "T1"):
        print(o["rank"], o["accession"], o["origin"], o["n_rows"])
    return 0


# ----------------------------------------------------------------------------------------------- selftest
class _FakeResp:
    def __init__(self, status, body, cl=None):
        self.status, self._b = status, body
        self.headers = {"Content-Length": str(len(body) if cl is None else cl)}

    def read(self, k):
        return self._b[:k]

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class _FakeOpener:
    def __init__(self, script):
        self.script, self.calls = list(script), []

    def open(self, req, timeout=None):
        self.calls.append(req.full_url)
        s = self.script.pop(0)
        if isinstance(s, Exception):
            raise s
        return s


def _rec_base(**kw):
    r = {"accession": "MTBLS123", "tranche": "T1", "extractor": "assembler", "extracted_utc": "2026-10-30T10:00:00Z",
         "sources_consulted": ["S2_PMC1.view.txt"]}
    r.update(kw)
    return r


def selftest():
    import tempfile
    fails, n = [], [0]

    def check(cond, msg):
        n[0] += 1
        if not cond:
            fails.append(msg)

    def refused(fn):
        try:
            fn()
            return False
        except Refused:
            return True

    # 1. guard
    acc = "MTBLS123"
    check(guard(mirror_url(acc)) == "index", "index URL allowed")
    check(guard(mirror_url(acc, "i_Investigation.txt")) == "isa_i", "i_ file allowed")
    for bad, why in [(mirror_url(acc, "m_MTBLS123_maf.tsv"), "m_ file"), (MTBLS + "MTBLS123/FILES/x.txt", "FILES/"),
                     ("https://example.org/x", "other host"), (EPMC + "PMC1/supplementaryFiles", "supplementary files"),
                     (MTBLS + "MTBLS123/../MTBLS9/i_x.txt", "traversal"), (MTBLS + "MTBLS123%2F..%2FMTBLS9/i_x.txt", "encoded traversal"),
                     (MTBLS + "MTBLS123/%2e%2e/i_x.txt", "encoded dots")]:
        check(refused(lambda: guard(bad)), f"guard refuses {why}")
    check(refused(lambda: guard(mirror_url(acc, "s_MTBLS123.txt"))), "s_ refused outside stratum A")
    check(guard(mirror_url(acc, "s_MTBLS123.txt"), allow_sa=True) == "isa_sa", "s_ allowed for stratum A")
    check(refused(lambda: mirror_url(acc, "../x.txt")), "bad file name refused")
    su = search_url('"MTBLS123"', 1000)
    check(guard(su) == "search" and "synonym=FALSE" in su and "resultType=core" in su and "pageSize=1000" in su
          and "query=%22MTBLS123%22" in su, "S3 search URL frozen parameters")
    check(guard(fulltext_url("PMC4567")) == "fulltext", "fullTextXML allowed")
    check(guard(search_url('DOI:"10.1038/s41467-024-45838-3"', 25)) == "search", "DOI search (%2F in the query) allowed")
    check(guard(search_url('"MTBLS123"', 1000, "AoJ/abc+==")) == "search", "cursor with / allowed in the query")

    # 2. transport
    with tempfile.TemporaryDirectory() as td:
        lp = os.path.join(td, "log.tsv")
        op = _FakeOpener([_FakeResp(200, b"abc", cl=10), _FakeResp(200, b"abcdefghij")])
        t = Transport(lp, opener=op, sleep=lambda s: None)
        ok, st, body, note, definite = t.get(mirror_url(acc))
        check(ok and definite and body == b"abcdefghij" and len(op.calls) == 2, "short body refused, retried once, then accepted")
        op.script = [urllib.error.HTTPError(mirror_url(acc), 302, "Found", {}, None)]
        ok, st, body, note, definite = t.get(mirror_url(acc))
        check(not ok and st == 302 and note == "HTTP_302" and definite, "redirect recorded as HTTP_302 (definite), not followed")
        op.script = [_FakeResp(200, b"x", cl=5), _FakeResp(200, b"y", cl=5)]
        ok, st, body, note, definite = t.get(mirror_url(acc))
        check(not ok and note.startswith("SHORT_BODY") and not definite, "persistent short body: indeterminate")
        op.script = [OSError("dns down"), OSError("dns down")]
        ok, st, body, note, definite = t.get(mirror_url(acc))
        check(not ok and note.startswith("TRANSPORT_FAILED:OSError:dns down") and not definite, "transport failure: reason logged, indeterminate")
        op.script = [urllib.error.HTTPError(EPMC, 404, "NF", {}, None)]
        ok, st, body, note, definite = t.get(fulltext_url("PMC9"))
        check(not ok and st == 404 and definite, "404: definite")
        op.script = [urllib.error.HTTPError(EPMC, 403, "Forbidden", {}, None)]
        ok, st, body, note, definite = t.get(fulltext_url("PMC9"))
        check(not ok and st == 403 and not definite, "re-review #6: 403 (firewall/proxy) is indeterminate")
        check(refused(lambda: t.get(mirror_url(acc, "m_x.tsv"))), "transport guard")
        # outage: 5 consecutive indeterminate failures on one host raise Outage (5xx counts)
        t.fail_run["www.ebi.ac.uk"] = 0
        op.script = [urllib.error.HTTPError(EPMC, 503, "x", {}, None)] * 20
        outage = False
        try:
            for _ in range(6):
                t.get(search_url("x", 25))
        except Outage:
            outage = True
        check(outage and t.fail_run["www.ebi.ac.uk"] == OUTAGE_K, "5xx counted; Outage raised at the 5th consecutive failure")
        t.close()
        log = open(lp).read().splitlines()
        check(log[0].startswith("utc\turl") and len(log) == 1 + 2 + 1 + 2 + 2 + 1 + 1 + 2 * OUTAGE_K, f"journal: one line per attempt ({len(log)})")

    # 3. parsing
    idx = b'<a href="?C=N">x</a><a href="/pub/">up</a><a href="i_Investigation.txt">i</a><a href=\'s_MTBLS123%20b.txt\'>s</a><a href="FILES/">F</a>'
    names = index_names(idx)
    check("i_Investigation.txt" in names and "s_MTBLS123 b.txt" in names and "FILES/" in names and "?C=N" not in names,
          "directory index names decoded once (double and single quotes)")
    inv_txt = ('STUDY PUBLICATIONS\n"Study PubMed ID"\t"31234567"\t""\t"n/a"\n'
               '"Study Publication DOI"\t"https://doi.org/10.1038/S41467-020-1"\t"doi:10.3390/metabo1"\t"none"\t"10.9/x"\n'
               'Study File Name\t"s_MTBLS123.txt"\nStudy Assay File Name\t"a_x.txt"\t"a_y.txt"\n')
    inv = parse_investigation(inv_txt)
    check(inv["pmids"] == ["31234567"] and inv["dois"] == ["10.1038/s41467-020-1", "10.3390/metabo1"]
          and inv["study_files"] == ["s_MTBLS123.txt"] and inv["assay_files"] == ["a_x.txt", "a_y.txt"], f"investigation parse {inv}")
    res = [{"id": "3", "source": "MED", "pmid": "3", "doi": "10.1/c", "firstPublicationDate": "2021-01-01"},
           {"id": "1", "source": "MED", "pmid": "1", "doi": "10.1/A", "firstPublicationDate": "2019-05-01"},
           {"id": "1b", "source": "PPR", "doi": "10.1/a", "firstPublicationDate": "2019-06-01"},
           {"id": "2", "source": "MED", "pmid": "2", "firstPublicationDate": "2020-01-01"},
           {"id": "2x", "source": "PMC", "pmid": "2", "pmcid": "PMC77", "firstPublicationDate": "2020-02-01"},
           {"id": "9", "source": "MED", "pmid": "9", "firstPublicationDate": "2017-01-01"},
           {"id": "4", "source": "MED", "pmid": "4", "firstPublicationDate": "2022-01-01"},
           {"id": "5", "source": "MED", "pmid": "5", "firstPublicationDate": "2023-01-01"},
           {"id": "6", "source": "MED", "pmid": "6"}, "not-a-dict"]
    sel = s3_select(res, {"pmid:9"})
    check([r["id"] for r in sel] == ["1", "2", "3", "4", "5"], f"S3 multi-identifier dedupe/sort/first-5 {[r['id'] for r in sel]}")
    check(accession_in_text("MTBLS12", "data in MTBLS12.") and not accession_in_text("MTBLS12", "see MTBLS123")
          and not accession_in_text("MTBLS12", "xMTBLS12"), "accession boundary match")

    # 4. views
    jats = b"""<article><front><article-meta><title-group><article-title>A study</article-title></title-group>
      <abstract><sec><title>Background</title><p>We studied 100 adults.</p></sec><sec><title>Results</title><p>Lactate rose 2-fold.</p></sec></abstract>
      </article-meta></front><body>
      <sec><title>Introduction</title><p>Background text.</p></sec>
      <sec sec-type="methods"><title>Materials and Methods</title><p>Serum was collected once from each of 100 participants&nbsp;at baseline.
        <fig><caption><p>Kynurenine elevated in cases</p></caption></fig> AT&T kits were used.
        <table-wrap><caption><p>Metabolite concentrations</p></caption><table><tr><td>88.1</td></tr></table></table-wrap></p>
        <sec><title>LC-MS analysis</title><p>Samples were run in one batch.</p></sec></sec>
      <sec><title>Metabolic differences between patient groups</title><p>Taurine 3.1 vs 1.2 uM.</p></sec>
      <sec><title>Findings</title><sec><title>Study participants</title><p>Of 120 enrolled, 100 analysed.</p></sec></sec>
      <sec><title>Results</title><p>Metabolite X was elevated (p=0.01).</p>
        <table-wrap><label>Table 2</label><caption><p>Differential metabolites</p></caption><table><tr><td>X</td><td>2.1</td></tr></table></table-wrap>
      </sec>
      <sec sec-type="supplementary-material"><title>Supplementary Material</title><p>Table S1 lists 3.3 uM values.</p></sec>
      <table-wrap><label>Table 1</label><caption><p>Baseline characteristics of participants</p></caption>
        <table><tr><th>n</th><td>100</td></tr></table></table-wrap>
      <table-wrap><label>Table 3</label><caption><p>Clinical and biochemical characteristics</p></caption><table><tr><td>glucose 8.1</td></tr></table></table-wrap>
      </body><back><sec sec-type="data-availability"><title>Data availability</title><p>Data are in MTBLS123.</p></sec></back></article>"""
    v = view_from_xml(jats, True)
    nv = norm_ws(v)
    check("[TITLE] A study" in v and "We studied 100 adults." in v and "Lactate" not in v, "view: title + abstract, structured-abstract Results dropped")
    check("collected once from each of 100 participants at baseline" in nv and "AT&T kits were used." in nv, "view: methods text, entities and bare & repaired")
    check("Kynurenine" not in v and "88.1" not in v and "Metabolite concentrations" not in v, "view: figures/tables nested in a methods paragraph excluded")
    check("one batch" in v, "view: methods subsection inherited")
    check("Taurine" not in v, "view: top-level section that only mentions 'patient' is not a methods section")
    check("Of 120 enrolled" not in v, "view: children of 'Findings' excluded")
    check("Metabolite X was elevated" not in v and "Differential metabolites" not in v and "2.1" not in v, "view: results excluded")
    check("3.3 uM" not in v, "view: supplementary-material section excluded")
    check("Baseline characteristics" in v and "glucose 8.1" not in v, "view: tables filtered by caption (biochemical excluded)")
    check("Data are in MTBLS123." in v and "Background text" not in v, "view: data availability in, introduction out")
    check("We studied 100 adults." not in view_from_xml(jats, False), "view: abstract excluded when DR-8 = no")
    check(view_from_xml(b"<article><body><p>AT&T data</p></body></article>", True) is not None, "bare ampersand does not crash")
    j2 = b"""<article><front><article-meta><abstract><p><bold>Background</bold> Adults were studied.</p>
      <p><bold>Results</bold> Citrate rose 2.1-fold.</p></abstract></article-meta></front><body>
      <sec><title>RESEARCH DESIGN AND METHODS</title><p>We enrolled 75 adults.</p>
        <sec><title>Outcome measures</title><p>The primary outcome was HbA1c at entry.</p></sec></sec>
      <sec><title>STAR&#9733;Methods</title><p>Key resources table text.</p></sec></body>
      <back><notes><sec><title>Data availability</title><p>Deposited as MTBLS9.</p></sec></notes></back></article>"""
    v2b = view_from_xml(j2, True)
    check("Adults were studied." in v2b and "Citrate" not in v2b, "re-review #5: bold-labelled abstract Results dropped")
    check("We enrolled 75 adults." in v2b and "primary outcome was HbA1c" in v2b, "re-review #7: ADA heading and nested 'Outcome measures' kept")
    check("Key resources table text." in v2b, "re-review #7: STAR-Methods kept")
    check("Deposited as MTBLS9." in v2b, "re-review #7: back > notes > sec Data availability kept")
    v3b = view_from_xml(b"<article><front><article-meta><abstract><p>Background: 40 adults. Results: lactate rose. Conclusions: x.</p>"
                        b"</abstract></article-meta></front><body/></article>", True)
    check("40 adults." in v3b and "lactate" not in v3b, "inline 'Results:' label truncates the abstract")
    pdf_txt = ("Title of paper\nAbstract: 80 patients.\nIntroduction\nintro text\n2. Methods\nEach of 80 patients gave one urine sample.\n"
               "participants with diabetes had higher\nglucose (8.1 vs 5.2 mmol/L)\n3. Results\nX up.\n")
    pv = view_from_pdf_text(pdf_txt, True)
    check("[VIEW_HEURISTIC_PDF]" in pv and "Abstract: 80 patients." in pv and "one urine sample" in pv and "X up." not in pv
          and "intro text" not in pv, "PDF view heuristic")
    pdf3 = "Intro\nIntroduction\nx\nResults\nTable 2\nParticipants\ncitrate 88.1 vs 60.2\nDiscussion\ny\nMethods\nEach of 30 donors gave one sample.\n"
    pv3 = view_from_pdf_text(pdf3, True)
    check("citrate 88.1" not in pv3 and "Each of 30 donors gave one sample." in pv3, "re-review #4: bare 'Participants' table header is not a heading")
    pdf2 = "Abstract text\nINTRODUCTION\nx\nRESULTS\nparticipants with diabetes had higher\nglucose 8.1\nMETHODS\nEach of 50 subjects gave one sample.\nReferences\n1. Ref 9.9\n"
    pv2 = view_from_pdf_text(pdf2, True)
    check("Each of 50 subjects gave one sample." in pv2 and "glucose 8.1" not in pv2 and "Ref 9.9" not in pv2,
          "PDF: lower-case wrapped line is not a heading; methods ends at References")

    # 5. order
    rows = [{"accession": f"MTBLS{i}", "origin": "M5_SOURCE_UNIQUE" if i % 7 else "A35_ALLQC_DIAGNOSTIC", "n_rows": "60"} for i in range(1, 121)]
    o1, o2 = order_rows(rows, "stageA-docrec-v1"), order_rows(list(reversed(rows)), "stageA-docrec-v1")
    check(o1 == o2, "order independent of input order")
    check(sum(o["in_T1"] for o in o1) == 40 and all(o["stratum"] == "M" for o in o1 if o["in_T1"]), "T1 = first 40 of stratum M")
    check(next(o["key"] for o in o1 if o["accession"] == "MTBLS1") == sha_bytes(b"stageA-docrec-v1|MTBLS1"), "key = SHA-256(salt|accession)")
    check(all(o1[i]["key"] <= o1[i + 1]["key"] for i in range(len(o1) - 1) if o1[i]["stratum"] == o1[i + 1]["stratum"]), "ascending keys")
    good_b = "NISTv2:2-1594646:" + "A" * 128
    check(SALT_RX.match(good_b) and not SALT_RX.match("NISTv2:2-1594646:" + "a" * 128) and not SALT_RX.match("anything"), "salt pattern")
    check([o["accession"] for o in order_rows(rows, good_b) if o["in_T1"]] != [o["accession"] for o in o1 if o["in_T1"]],
          "different salt, different sample")
    check(beacon_url("2026-10-30T12:00:00Z") == BEACON + "1793361600000", "beacon URL in Unix milliseconds")

    # 6. grading
    U = {"files": {"S1_i_Investigation.txt": {"role": "S1", "usable": True},
                   "S2_PMC1.view.txt": {"role": "S2_view", "usable": True},
                   "S3_PMC2.view.txt": {"role": "S3_view", "usable": True},
                   "S3_PMC3.view.txt": {"role": "S3_view", "usable": False},
                   "S2M_paper.view.txt": {"role": "S2M_view", "usable": True},
                   "S2_PMC1.xml": {"role": "S2_fulltext", "usable": False}}}
    V = {"S1_i_Investigation.txt": norm_ws("Study Description: serum from 100 patients. Factor: visit baseline, 6 months"),
         "S2_PMC1.view.txt": norm_ws("[SECTION: Methods] Serum samples from 100 patients were analysed by LC-MS. "
                                     "Each participant provided a single fasting sample. Of 120 enrolled patients, 100 were analysed. "
                                     "Samples were collected at baseline and after 6 months. We included 60 cases and 40 controls. "
                                     "One hundred patients were enrolled."),
         "S3_PMC2.view.txt": norm_ws("[SECTION: Methods] We reanalysed MTBLS123 serum data from 100 patients. Serum data from 100 patients."),
         "S3_PMC3.view.txt": "irrelevant text here",
         "S2M_paper.view.txt": norm_ws("[METHODS_PDF] Serum samples from 100 patients were analysed by LC-MS. Each participant provided a single fasting sample.")}

    def E(src, q):
        return [{"source": src, "locator": "Methods", "quote": q}]

    def rec(**kw):
        r = _rec_base(Q1={"answer": "yes", "evidence": E("S2_PMC1.view.txt", "Serum samples from 100 patients")},
                      Q2={"answer": "established", "evidence": E("S2_PMC1.view.txt", "Serum samples from 100 patients were analysed by LC-MS.")},
                      Q3={"count": 100, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "100 patients were analysed by LC-MS")},
                      Q4={"answer": "explicit_one", "evidence": E("S2_PMC1.view.txt", "Each participant provided a single fasting sample.")},
                      Q6={"flags": [], "confirmed_overlap": [], "evidence": []}, Q7={"answer": "none", "evidence": []})
        r.update(kw)
        return r

    def g(r, n=100):
        return grade_record(r, n, U, V)[:2]
    check(g(rec()) == ("DOC_RECOVERED_STRICT", "PROVISIONAL"), f"strict case {g(rec())}")
    check(g(rec(Q4={"answer": "not_stated", "evidence": []})) == ("DOC_RECOVERED_COUNT", "SENSITIVITY_ONLY"), "count case")
    check(g(rec(Q3={"count": 120, "scope": "study_enrollment", "evidence": E("S2_PMC1.view.txt", "Of 120 enrolled patients")})) ==
          ("DOC_UNRESOLVED", "NO_ASSAY_SCOPED_COUNT"), "owner's 120-enrolled / 100-rows case is unresolved")
    check(g(rec(), n=150) == ("DOC_UNRESOLVED", "ROWS_EXCEED_COUNT"), "rows > count: unresolved, not contradicted (C2)")
    check(g(rec(), n=80) == ("DOC_UNRESOLVED", "COUNT_EXCEEDS_ROWS"), "count > rows: unresolved")
    rep = rec(Q4={"answer": "explicit_repeated", "evidence": E("S2_PMC1.view.txt", "Samples were collected at baseline and after 6 months.")})
    check(g(rep) == ("DOC_CONTRADICTED", "REPEATED_SAMPLING"), "explicit repetition: contradicted")
    check(g(dict(rep, Q2={"answer": "not_established", "evidence": []})) == ("DOC_UNRESOLVED", "Q2_NOT_ESTABLISHED"),
          "repetition without assay link: unresolved")
    check(g(rec(Q1={"answer": "not_stated", "evidence": []})) == ("DOC_UNRESOLVED", "Q1_NOT_STATED"), "Q1 not stated")
    s1only = rec(Q2={"answer": "established", "evidence": E("S1_i_Investigation.txt", "serum from 100 patients.")},
                 Q3={"count": 100, "scope": "selected_assay", "evidence": E("S1_i_Investigation.txt", "serum from 100 patients.")})
    check(g(s1only) == ("DOC_UNRESOLVED", "S1_ONLY_EVIDENCE"), "S1-only evidence never recovers (A29)")
    check(g(rec(Q4={"answer": "explicit_repeated", "evidence": E("S1_i_Investigation.txt", "visit baseline, 6 months")})) ==
          ("DOC_CONTRADICTED", "REPEATED_SAMPLING"), "S1 may supply contradiction evidence")
    s3q2 = rec(Q2={"answer": "established", "evidence": E("S3_PMC2.view.txt", "Serum data from 100 patients.")})
    check(g(s3q2) == ("DOC_UNRESOLVED", "Q2_S3_WITHOUT_ACCESSION"), "S3-only Q2 needs the accession in the quote")
    check(g(rec(Q2={"answer": "established", "evidence": E("S3_PMC2.view.txt", "We reanalysed MTBLS123 serum data from 100 patients.")})) ==
          ("DOC_RECOVERED_STRICT", "PROVISIONAL"), "S3-only Q2 with the accession in the quote is admissible")
    vol = rec(Q3={"count": 100, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "Plasma volume was 100 µL per sample.")})
    V["S2_PMC1.view.txt"] += " Plasma volume was 100 µL per sample."
    check(grade_record(vol, 100, U, V)[0] == "INVALID", "re-review #3: '100 µL' does not support a count of 100")
    mixed = rec(Q3={"count": 60, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "Serum samples from 100 patients"),
                    "count_derivation": [{"text": "100 patients", "value": 60}]})
    check(grade_record(mixed, 60, U, V)[0] == "INVALID", "re-review #3: derivation text mixing digits and words rejected")
    check(words_to_int("One hundred and twenty-two") == 122 and words_to_int("sixty") == 60 and words_to_int("many") is None,
          "number words parsed")
    check(g(dict(s3q2, Q4=rep["Q4"])) == ("DOC_CONTRADICTED", "REPEATED_SAMPLING"),
          "review #4: S3 restriction sits at rule 8, after the contradiction rules")
    pooled = rec(Q4={"answer": "not_stated", "evidence": E("S2_PMC1.view.txt", "Each participant provided a single fasting sample.")},
                 Q7={"answer": "one_to_one", "evidence": E("S1_i_Investigation.txt", "serum from 100 patients.")})
    check(g(pooled) == ("DOC_RECOVERED_COUNT", "SENSITIVITY_ONLY"), "review #3: Q4/Q7 evidence not pooled; S1-only Q7 gives COUNT")
    check(g(rec(Q4={"answer": "explicit_one", "evidence": E("S1_i_Investigation.txt", "serum from 100 patients.")})) ==
          ("DOC_RECOVERED_COUNT", "SENSITIVITY_ONLY"), "one-sample evidence only from S1: COUNT")
    for bad, why in [(rec(Q1={"answer": "yes", "evidence": E("S2_PMC1.view.txt", "Serum samples from 999 patients")}), "fabricated quote"),
                     (rec(Q1={"answer": "yes", "evidence": E("S3_PMC3.view.txt", "irrelevant text here")}), "unusable S3 source"),
                     (rec(Q1={"answer": "yes", "evidence": E("S2_PMC1.xml", "Serum samples from 100")}), "full document (not view)"),
                     (rec(Q1={"answer": "yes", "evidence": E("S2_PMC1.view.txt", "S")}), "one-character quote"),
                     (rec(Q1={"answer": "yes", "evidence": E("S2_PMC1.view.txt", "[SECTION: Methods] Serum samples")}), "quote spans a marker"),
                     (rec(Q3={"count": 100, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "Serum samples from 100 patients")[:0]
                          + E("S2_PMC1.view.txt", "Each participant provided a single fasting sample.")}), "count absent from its quote"),
                     (rec(Q3={"count": "100", "scope": "selected_assay", "evidence": []}), "count type"),
                     (rec(Q3={"count": True, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "100 patients were analysed")}), "bool count"),
                     (rec(Q4={"answer": "explicit_one", "evidence": []}), "missing evidence"),
                     (rec(Q1="yes"), "review #9: Q1 not an object"),
                     (rec(Q1={"answer": "yes", "evidence": [{"source": "S2_PMC1.view.txt", "locator": "M", "quote": 7}]}), "review #9: integer quote"),
                     (rec(Q6={"flags": [], "confirmed_overlap": True, "evidence": []}), "review #9: boolean confirmed_overlap"),
                     (rec(Q6={"flags": [], "confirmed_overlap": ["MTBLS123"], "evidence": []}), "self-overlap"),
                     (_rec_base(tranche="A35", **{k: v for k, v in rec().items() if k.startswith("Q")}), "wrong tranche field"),
                     ({k: v for k, v in rec().items() if k != "extracted_utc"}, "missing extracted_utc"),
                     ("not an object", "record not an object")]:
        out = grade_record(bad, 100, U, V)
        check(out[0] == "INVALID", f"INVALID: {why} -> {out[:2]}")
    sumrec = rec(Q3={"count": 100, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "We included 60 cases and 40 controls."),
                     "count_derivation": [{"text": "60", "value": 60}, {"text": "40", "value": 40}]})
    check(g(sumrec) == ("DOC_RECOVERED_STRICT", "PROVISIONAL"), "count as a sum of stated disjoint groups")
    check(grade_record(dict(sumrec, Q3=dict(sumrec["Q3"], count_derivation=[{"text": "60", "value": 60}, {"text": "40", "value": 41}])), 101, U, V)[0]
          == "INVALID", "count_derivation must match its text")
    words = rec(Q3={"count": 100, "scope": "selected_assay", "evidence": E("S2_PMC1.view.txt", "One hundred patients were enrolled."),
                    "count_derivation": [{"text": "One hundred", "value": 100}]})
    gw = grade_record(words, 100, U, V)
    check(gw[0] == "DOC_RECOVERED_STRICT" and "NUMBER_WORDS" in gw[3], "number words accepted with a NUMBER_WORDS flag")
    pdfrec = rec(Q2={"answer": "established", "evidence": E("S2M_paper.view.txt", "Serum samples from 100 patients were analysed by LC-MS.")})
    gp = grade_record(pdfrec, 100, U, V)
    check(gp[0] == "DOC_RECOVERED_STRICT" and "VIEW_HEURISTIC_PDF" in gp[3], "review #13: owner-supplied PDF evidence flagged")
    # V-grading
    vrec = _rec_base(tranche="A35", V1={"answer": "pass", "evidence": E("S2_PMC1.view.txt", "Serum samples from 100 patients")},
                     V2={"answer": "pass", "evidence": E("S2_PMC1.view.txt", "100 patients were analysed")},
                     V3={"answer": "pass", "count": 100, "evidence": E("S2_PMC1.view.txt", "100 patients were analysed")},
                     V4={"answer": "not_established", "evidence": []})
    check(grade_v_record(vrec, U, V)[0] == "UNVERIFIABLE", "V: not established -> UNVERIFIABLE")
    vrec["V4"] = {"answer": "pass", "evidence": E("S2_PMC1.view.txt", "Each participant provided a single fasting sample.")}
    gv = grade_v_record(vrec, U, V)
    check((gv[0], gv[2]) == ("VERIFIED", 100), "V: all pass -> VERIFIED with count")
    vrec["V1"] = {"answer": "fail", "evidence": E("S1_i_Investigation.txt", "serum from 100 patients.")}
    check(grade_v_record(vrec, U, V)[0] == "FAILED", "V: any fail -> FAILED")
    check(grade_v_record(dict(vrec, V2="pass"), U, V)[0] == "INVALID", "V: malformed -> INVALID, no crash")

    # 7. bounds and futility
    check(abs(cp_upper(5, 40) * 371 - 90.9) < 0.1 and abs(cp_upper(6, 40) * 371 - 101.9) < 0.1, "CP upper pools 90.9 / 101.9")
    check(abs(cp_lower(16, 40) * 371 - 99.9) < 0.1 and abs(cp_lower(25, 40) * 371 - 179.1) < 0.1, "CP lower pools 99.9 / 179.1")
    check(hyper_upper(5, 371, 40) == 88 and hyper_upper(6, 371, 40) == 99, "hypergeometric upper K 88 / 99")
    check(hyper_lower(15, 371, 40) == 94 and hyper_lower(25, 371, 40) == 182, "hypergeometric lower K 94 / 182")
    for opt in ("a", "b"):
        check(futility(5, salt_option=opt)["DR7_futility_reached"] and not futility(6, salt_option=opt)["DR7_futility_reached"],
              f"DR-7 threshold x <= 5 under option {opt}")
    check(components(["A", "B", "C"], [("A", "B"), ("B", "Z")]) == 2, "confirmed overlap collapses")

    # 8. authorisation lines
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, LOCK), "w").write("{}")
        lsha = sha_file(os.path.join(td, LOCK))
        lkb = {"salt_option": "b", "beacon_url": BEACON + "1"}

        def line(salt, lk_sha=lsha):
            return (f"DOCREC FETCH AUTHORISED 2026-10-30 tranche=T1 candidates={CANDIDATES_SHA} protocol={PROTOCOL} "
                    f"lock={lk_sha} salt={salt} owner=Anas\n")
        open(os.path.join(td, AUTH_FILE), "w").write(line("stageA-docrec-v1"))
        check(refused(lambda: check_fetch_auth(lkb, "T1", base=td)), "option (b) refuses the fixed literal")
        open(os.path.join(td, AUTH_FILE), "a").write(line(good_b))
        check(not refused(lambda: check_fetch_auth(lkb, "T1", base=td)), "the last matching line governs (corrected salt)")
        open(os.path.join(td, AUTH_FILE), "a").write(line(good_b, "0" * 64))
        check(refused(lambda: check_fetch_auth(lkb, "T1", base=td)), "authorisation must name the signed lock's hash")
        open(os.path.join(td, AUTH_FILE), "a").write("DOCREC SEAL SIGNED 2026-10-30 tranche=T1 universe=" + "a" * 64 + " owner=Anas\n")
        check(seal_signed("T1", base=td) and not seal_signed("A35", base=td), "seal-signed detection by tranche")
        check(refused(lambda: check_seal_auth("T1", "b" * 64, base=td)) and not refused(lambda: check_seal_auth("T1", "a" * 64, base=td)),
              "seal line must name the universe")

    # 9. integration: plan -> fetch (scripted network) -> seal -> grade
    _integration(check)

    print(f"selftest: {n[0] - len(fails)}/{n[0]} checks passed")
    for f in fails:
        print("  FAIL:", f)
    return 0 if not fails else 1


def _integration(check):
    """End-to-end run of the real commands against a scripted mirror / Europe PMC (no network)."""
    import contextlib
    import tempfile
    global CANDIDATES_SHA, _TEST_OPENER, _TEST_SLEEP
    saved = (CANDIDATES_SHA, _TEST_OPENER, _TEST_SLEEP, os.getcwd())

    def jats(acc, n, rep, amp=False):
        extra = " Samples were taken at two visits." if rep else ""
        amp_s = " AT&T kits." if amp else ""
        return (f"<article><front><article-meta><title-group><article-title>Study {acc}</article-title></title-group></article-meta></front>"
                f"<body><sec><title>Methods</title><p>Plasma from {n} participants was analysed; each participant provided one sample."
                f"{extra}{amp_s} Data: {acc}.</p></sec>"
                f"<sec><title>Results</title><p>Glucose was higher (p=0.001).</p></sec></body></article>").encode()

    class Net:
        def __init__(self):
            self.urls, self.epmc_down, self.idx404, self.bad_json_once = [], False, set(), set()

        def open(self, req, timeout=None):
            u = req.full_url
            self.urls.append(u)
            guard(u, allow_sa=True)
            if u.startswith(EPMC) and self.epmc_down:
                raise urllib.error.HTTPError(u, 503, "down", {}, None)
            m = re.match(re.escape(MTBLS) + r"(MTBLS\d+)/(.*)$", u)
            if m:
                acc, name = m.group(1), urllib.parse.unquote(m.group(2))
                num = int(acc[5:])
                if name == "" and acc in self.idx404:
                    raise urllib.error.HTTPError(u, 404, "NF", {}, None)
                if name == "":
                    return _FakeResp(200, (f'<a href="i_Investigation.txt">i</a><a href="s_{acc}.txt">s</a><a href="a_{acc}.txt">a</a>'
                                           f'<a href="m_{acc}.tsv">m</a><a href="FILES/">F</a>').encode())
                if name.startswith("i_"):
                    doi = f"https://doi.org/10.9999/{num}" if num % 5 == 0 else ""
                    return _FakeResp(200, (f'"Study PubMed ID"\t"{num}"\n"Study Publication DOI"\t"{doi}"\n"Study File Name"\t"s_{acc}.txt"\n'
                                           f'"Study Assay File Name"\t"a_{acc}.txt"\n').encode())
                return _FakeResp(200, f"Sample Name\tCharacteristics[Organism]\n{acc}_1\thomo sapiens\n".encode())
            q = urllib.parse.parse_qs(urllib.parse.urlparse(u).query)
            if u.startswith(EPMC + "search"):
                query = q["query"][0]
                mm = re.match(r"EXT_ID:(\d+) AND SRC:MED", query)
                if mm:
                    k = mm.group(1)
                    if k in self.bad_json_once:
                        self.bad_json_once.discard(k)
                        return _FakeResp(200, b'{"hitCount": 1, "resultL')  # truncated, no Content-Length mismatch
                    return _FakeResp(200, json.dumps({"hitCount": 1, "resultList": {"result": [
                        {"id": k, "source": "MED", "pmid": k, "pmcid": f"PMC{k}", "doi": f"10.9999/{k}", "isOpenAccess": "Y",
                         "firstPublicationDate": "2020-01-01"}]}}).encode())
                md = re.match(r'DOI:"10\.9999/(\d+)"$', query)
                if md:
                    k = md.group(1)
                    return _FakeResp(200, json.dumps({"hitCount": 1, "resultList": {"result": [
                        {"id": k, "source": "MED", "pmid": k, "pmcid": f"PMC{k}", "doi": f"10.9999/{k}", "isOpenAccess": "Y"}]}}).encode())
                k = query.strip('"')[5:]
                if q["cursorMark"][0] == "*":
                    return _FakeResp(200, json.dumps({"hitCount": 3, "nextCursorMark": "AoJ/x+1", "resultList": {"result": [
                        {"id": "7000" + k, "source": "PPR", "doi": f"10.9999/7000{k}", "firstPublicationDate": "2019-01-01"}]}}).encode())
                return _FakeResp(200, json.dumps({"hitCount": 3, "nextCursorMark": "AoJ/x+1", "resultList": {"result": [
                    {"id": k, "source": "MED", "pmid": k, "pmcid": f"PMC{k}", "doi": f"10.9999/{k}"},
                    {"id": "8000" + k, "source": "MED", "pmid": "8000" + k, "pmcid": f"PMC8000{k}", "doi": f"10.9999/8000{k}",
                     "firstPublicationDate": "2021-01-01"}]}}).encode())
            mm = re.match(re.escape(EPMC) + r"PMC(\d+)/fullTextXML$", u)
            if mm:
                k = mm.group(1)
                if k.startswith("8000") and len(k) > 4:  # S3 paper: accession only for odd numbers; contains a bare '&'
                    return _FakeResp(200, jats("MTBLS" + k[4:] if int(k[4:]) % 2 else "MTBLS0", 60, False, amp=True))
                return _FakeResp(200, jats("MTBLS" + k, 60, int(k) % 10 == 3))
            raise urllib.error.HTTPError(u, 404, "NF", {}, None)

    def fields(acc, srcs):
        return {"accession": acc, "tranche": "T1", "extractor": "test", "extracted_utc": "2026-10-30T10:00:00Z", "sources_consulted": srcs}

    try:
        with tempfile.TemporaryDirectory() as td, contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            os.chdir(td)
            rows = [("MTBLS%d" % i, "M5_SOURCE_UNIQUE" if i <= 60 else "A35_ALLQC_DIAGNOSTIC", 60) for i in range(1, 64)]
            buf = "accession,origin,census_grade,n_sources,n_rows,platform_class,matrix,census_version,title\n" + "".join(
                f"{a},{o},,0,{n},MS,plasma,0.1.7,t\n" for a, o, n in rows)
            open(CANDIDATES, "w").write(buf)
            CANDIDATES_SHA = sha_bytes(buf.encode())
            open(PROTOCOL, "w").write("protocol text\n")
            _TEST_SLEEP = lambda s: None  # noqa: E731
            net = Net()
            _TEST_OPENER = net
            t1 = [o["accession"] for o in tranche_members(order_rows(load_candidates(CANDIDATES), "stageA-docrec-v1"), "T1")]
            check(main(["fetch", "--tranche", "T1"]) == 1, "fetch refused without a lock")
            check(main(["plan", "--salt-option", "b", "--beacon-not-before", "2020-01-01T00:00:00Z", "--abstract", "yes", "--manual", "no"]) == 1,
                  "review #11: plan refuses a beacon time in the past")
            main(["plan", "--salt-option", "a", "--abstract", "yes", "--manual", "no"])
            check(main(["fetch", "--tranche", "T1"]) == 1, "fetch refused with an unsigned lock")
            lk = json.load(open(LOCK))
            lk["owner_signoff"] = "Test 2026-10-29T00:00:00Z"
            json.dump(lk, open(LOCK, "w"))
            check(main(["fetch", "--tranche", "T1"]) == 1, "fetch refused without an authorisation line")

            def auth(tr):
                return (f"DOCREC FETCH AUTHORISED 2026-10-29 tranche={tr} candidates={CANDIDATES_SHA} protocol={PROTOCOL} "
                        f"lock={sha_file(LOCK)} salt=stageA-docrec-v1 owner=Test\n")
            open(AUTH_FILE, "w").write(auth("T1"))
            # run 1: Europe PMC down -> per-host outage stop, partial manifests incomplete
            net.epmc_down = True
            check(main(["fetch", "--tranche", "T1"]) == 2, "review #1: Europe PMC outage stops the run (exit 2)")
            store = os.path.join(OUTDIR, "store")
            mans = {a: json.load(open(os.path.join(store, a, "manifest.json"))) for a in t1 if os.path.exists(os.path.join(store, a, "manifest.json"))}
            check(mans and not any(m["complete"] for m in mans.values()), f"review #1: no study marked complete during the outage ({len(mans)})")
            check(main(["seal", "--tranche", "T1"]) == 1, "review #2: seal refuses an incomplete tranche")
            for _ in range(3):
                check(main(["fetch", "--tranche", "T1"]) == 2, "outage persists across reruns: exit 2 each time")
            mans = {a: json.load(open(os.path.join(store, a, "manifest.json"))) for a in t1 if os.path.exists(os.path.join(store, a, "manifest.json"))}
            check(not any(m["complete"] for m in mans.values()) and all(m["attempt"] == 0 for m in mans.values()),
                  "re-review #1: outage-stopped runs consume no attempt and finalise nothing")
            # run 2: recovered; one index 404 (definite), one truncated JSON (indeterminate once)
            net.epmc_down = False
            net.idx404.add(t1[7])
            net.bad_json_once.add(t1[9][5:])
            code2 = main(["fetch", "--tranche", "T1"])
            check(code2 == 3, f"truncated JSON leaves one study incomplete: exit 3 ({code2})")
            m9 = json.load(open(os.path.join(store, t1[9], "manifest.json")))
            check(not m9["complete"] and any("BAD_JSON" in x for x in m9["notes"]), "truncated JSON is indeterminate")
            check(main(["fetch", "--tranche", "T1"]) == 0, "rerun completes the tranche")
            m7 = json.load(open(os.path.join(store, t1[7], "manifest.json")))
            check(m7["complete"] and any(x.startswith("S1_INDEX 404") for x in m7["notes"]), "definite S1 failure: manifest written, complete")
            check(all((lambda m: m["complete"] and m["attempt"] == 1 and len(m["s2"]) == 1 and m["s2"][0]["fulltext"] == "OK")(
                json.load(open(os.path.join(store, a, "manifest.json")))) for a in mans if a != t1[7]),
                "studies interrupted by the outage were fully re-fetched after recovery (S2 present)")
            check(not any(("/m_" in u) or ("FILES" in u) or ("/s_" in u) or ("/a_" in u) for u in net.urls),
                  "T1 never requests m_, FILES/, s_ or a_")
            n_before = len(net.urls)
            check(main(["fetch", "--tranche", "T1"]) == 0 and len(net.urls) == n_before, "resume skips complete studies")
            acc0 = next(a for a in t1 if a != t1[7])
            man = json.load(open(os.path.join(store, acc0, "manifest.json")))
            check(len(man["s2"]) == 1 and len(man["s3"]) == 2 and man["s2"][0]["fulltext"] == "OK", f"S2 and S3 resolved {man['s3']}")
            check(any("query=DOI%3A%2210.9999%2F" in u for u in net.urls), "DOI search issued (URL-encoded slash passes the guard)")
            check(any("cursorMark=AoJ%2Fx%2B1" in u for u in net.urls), "S3 pagination follows nextCursorMark")
            a5 = next(a for a in t1 if int(a[5:]) % 5 == 0 and a != t1[7])
            m5 = json.load(open(os.path.join(store, a5, "manifest.json")))
            check(len(m5["investigation"]["dois"]) == 1 and len(m5["s2"]) == 1, "PMID and DOI of one paper de-duplicated in S2")
            view = open(os.path.join(store, acc0, f"S2_PMC{acc0[5:]}.view.txt")).read()
            check("each participant provided one sample" in view and "Glucose" not in view, "stored view is value-blind")
            check(main(["seal", "--tranche", "T1"]) == 0, "seal")
            up = os.path.join(OUTDIR, "docrec_T1_evidence_universe.json")
            u = json.load(open(up))
            check(sorted(u["studies"]) == sorted(t1) and t1[7] in u["studies"], "review #2: universe = exactly the 40 tranche members")
            open(AUTH_FILE, "a").write(f"DOCREC SEAL SIGNED 2026-10-29 tranche=T1 universe={sha_file(up)} owner=Test\n")
            check(main(["seal", "--tranche", "T1"]) == 1, "no resealing after the seal is signed")
            check(main(["fetch", "--tranche", "T1"]) == 1, "no fetch after the seal is signed")
            open(AUTH_FILE, "a").write(auth("A35"))
            check(main(["fetch", "--tranche", "A35"]) == 0, "A35 fetch completes after T1 is sealed")
            check(any("/s_MTBLS61.txt" in u for u in net.urls) and not any("search?query=%22MTBLS61" in u for u in net.urls),
                  "A35 fetches s_/a_ and runs no S3 search")
            check(os.path.exists(order_path("T1")) and os.path.exists(order_path("A35")), "re-review #2: one order file per tranche")
            recs, strict = [], []
            for acc in t1:
                k = acc[5:]
                src = f"S2_PMC{k}.view.txt"
                rep = int(k) % 10 == 3
                if acc != t1[7] and (int(k) % 4 == 0 or rep):
                    r = dict(fields(acc, [src]),
                             Q1={"answer": "yes", "evidence": [{"source": src, "locator": "Methods", "quote": "Plasma from 60 participants was analysed"}]},
                             Q2={"answer": "established", "evidence": [{"source": src, "locator": "Methods", "quote": "Plasma from 60 participants was analysed"}]},
                             Q3={"count": 60, "scope": "selected_assay", "evidence": [{"source": src, "locator": "Methods", "quote": "Plasma from 60 participants"}]},
                             Q4={"answer": "explicit_repeated" if rep else "explicit_one", "evidence": [{"source": src, "locator": "Methods",
                                 "quote": "Samples were taken at two visits." if rep else "each participant provided one sample."}]},
                             Q6={"flags": [], "confirmed_overlap": [], "evidence": []}, Q7={"answer": "none", "evidence": []})
                    if not rep:
                        strict.append(acc)
                else:
                    r = dict(fields(acc, []), Q1={"answer": "not_stated", "evidence": []}, Q2={"answer": "not_established", "evidence": []},
                             Q3={"count": None, "scope": "none", "evidence": []}, Q4={"answer": "not_stated", "evidence": []},
                             Q6={"flags": [], "confirmed_overlap": [], "evidence": []}, Q7={"answer": "none", "evidence": []})
                recs.append(r)
            ep = os.path.join(OUTDIR, "docrec_T1_evidence.jsonl")
            open(ep, "w").write("".join(json.dumps(r) + "\n" for r in recs[:-1]))
            main(["grade", "--tranche", "T1"])
            rp = os.path.join(OUTDIR, "docrec_T1_report.json")
            rep_ = json.load(open(rp))
            check(rep_["missing_records"] == [t1[-1]] and str(rep_["futility"]).startswith("NOT COMPUTED"), "missing record blocks futility")
            open(ep, "w").write("".join(json.dumps(r) + "\n" for r in recs))
            main(["grade", "--tranche", "T1"])
            rep_ = json.load(open(rp))
            n_rep = sum(int(a[5:]) % 10 == 3 for a in t1 if a != t1[7])
            check(str(rep_["futility"]).startswith("NOT COMPUTED") and len(rep_["pending_adjudication"]) == len(strict) + n_rep,
                  "futility withheld while adjudication is pending")
            res = list(csv.DictReader(open(os.path.join(OUTDIR, "docrec_T1_results.csv"))))
            adj = "﻿accession,provisional_grade,record_sha256,confirmed,note\n" + "".join(
                f"{r['accession']},{r['provisional_grade']},{r['record_sha256']},yes,ok\n" for r in res if r["reason"] == "PROVISIONAL")
            open(os.path.join(OUTDIR, "docrec_T1_adjudication.csv"), "w", encoding="utf-8").write(adj)
            main(["grade", "--tranche", "T1"])
            rep_ = json.load(open(rp))
            fx = rep_["futility"]
            check(isinstance(fx, dict) and fx["x"] == len(strict) and fx["DR7_futility_reached"] == (len(strict) <= 5),
                  f"report x = {len(strict)} adjudicated STRICT and DR-7 applied ({fx}); BOM-prefixed adjudication CSV read")
            check(any(k.startswith("DOC_CONTRADICTED") for k in rep_["grade_counts"]), "repeated sampling graded CONTRADICTED")
            # review #10: a re-extracted record does not inherit an earlier adjudication
            changed = [dict(r) for r in recs]
            i_s = t1.index(strict[0])
            changed[i_s]["notes"] = "re-extracted"
            open(ep, "w").write("".join(json.dumps(r) + "\n" for r in changed))
            main(["grade", "--tranche", "T1"])
            rep_ = json.load(open(rp))
            check(strict[0] in rep_["pending_adjudication"] and strict[0] in rep_["stale_adjudication_rows"],
                  "review #10: adjudication bound to the record hash; stale rows reported")
            with open(os.path.join(store, acc0, f"S2_PMC{acc0[5:]}.view.txt"), "a") as f:
                f.write("tampered\n")
            check(main(["grade", "--tranche", "T1"]) == 1, "grade refuses a store changed after sealing")
    finally:
        os.chdir(saved[3])
        CANDIDATES_SHA, _TEST_OPENER, _TEST_SLEEP = saved[0], saved[1], saved[2]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    p = sub.add_parser("plan")
    p.add_argument("--salt-option", choices=["a", "b"], required=True)
    p.add_argument("--beacon-not-before")
    p.add_argument("--abstract", choices=["yes", "no"], required=True)
    p.add_argument("--manual", choices=["yes", "no"], required=True)
    p = sub.add_parser("order")
    p.add_argument("--salt", required=True)
    for c in ("fetch", "register", "seal", "grade"):
        p = sub.add_parser(c)
        p.add_argument("--tranche", choices=["T1", "A35"], required=True)
    a = ap.parse_args(argv)
    try:
        return {"selftest": selftest, "plan": lambda: cmd_plan(a), "order": lambda: cmd_order(a),
                "fetch": lambda: cmd_fetch(a), "register": lambda: cmd_register(a), "seal": lambda: cmd_seal(a),
                "grade": lambda: cmd_grade(a)}[a.cmd]()
    except Refused as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
