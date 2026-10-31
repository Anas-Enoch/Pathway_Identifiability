#!/usr/bin/env python3
"""
stageA_census.py  v0.1.7   Stage A — A18: Route M (MetaboLights) value-blind cohort census over the A17 universe.

Rules implemented verbatim from EXTERNAL_DATA_FEASIBILITY_PROTOCOL §2 (D1–D4), §3.3–3.6, §5, §6, via stagea_isa (frozen
recogniser and validators). Nothing here reads a metabolite-assignment file, a raw/derived file or any measured value:
the only requests per study are the mirror folder index (names only) and the i_/s_/a_ metadata files it lists.

Universe: stageA_universe_manifest.csv from A17 (verdict BOUNDED/ESTABLISHED; rows with status UNIVERSE only; a
PROVISIONAL_ manifest is refused), SHA-256 verified against the value hashed into the lock; exemplars (A7) excluded.

Per study (mirror only; REST never used by the census):
  FTPDIR -> class; not INDEX_ACCESSIBLE -> status INDEX_<class>, not screened
  FETCH every i_/s_/a_ the index lists (guards: name class, host, redirect, size, binary); validate per class
  i_: Study Title/Description/Design Type/Factor Name/Assay File Name<->Technology Type/Person Affiliation
  §3.4 assay choice: MS or NMR a_ with most rows (tie: lowest name); none -> EXCLUDED (no MS/NMR assay)
  §3.5 M1 human -> M2 QC removal -> M3 identity -> M4 matrix -> M5 D1 (R1 / R2 via stagea_isa.subject_evidence)
       -> M6 temporal (two-level) -> M7 n>=50 (inside M5) -> M8 FEATURE_COUNT_UNKNOWN (A3) -> M9 touched
Failure policy (fixed before execution; owner decision A20): a transport failure, refusal or validator failure on ONE
study marks that study (TRANSPORT_FAILED / FORMAT_REFUSED / INDEX_*) and the census continues; nothing is retried
within a run; a second pass over TRANSPORT_FAILED studies happens only as a separate, disclosed `--retry-pass` run on the
same lock. Progress is checkpointed per study (JSONL); a restart skips completed accessions, so a run is resumable
without re-requesting anything already recorded.
Outputs (--out): stageA_M_cohorts.csv, stageA_M_exclusions.csv, stageA_M_units.csv, stageA_M_summary.json,
stageA_M_files.csv (every fetched file: accession, name, bytes, sha256), stageA_census_fetch_log.txt,
stageA_census_checkpoint.jsonl, SHA256SUMS.txt.
"""
from __future__ import annotations

import argparse, csv, hashlib, http.client, json, re, socket, sys, time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, unquote
from urllib.request import HTTPRedirectHandler, Request, build_opener

here = Path(__file__).resolve().parent
sys.path.insert(0, str(here))
import stagea_isa as isa  # noqa: E402

VERSION = "0.1.7"
PROTOCOL = "EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0.7.7.md"
AUTH_FILE = "stageA_census_authorisation.txt"
AUTH_LINE = re.compile(r"LIVE STAGE A CENSUS AUTHORISED (\d{4}-\d{2}-\d{2}) universe=([0-9a-f]{64}) protocol=(EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0\.\d+(?:\.\d+)?\.md) owner=(\S+)(?: pilot=(\d+))?")
MANIFEST = "stageA_universe_manifest.csv"
FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"
EXEMPLARS = ("MTBLS1", "MTBLS2")            # A7: excluded from every census count
MIN_N = 50                                  # A5
PAUSE = 0.5                                 # seconds between requests (A21)
MAX_BYTES, TIMEOUT = 20_000_000, 300
OUTAGE_K = 5     # v0.1.7 (D-full-1): K consecutive studies ending TRANSPORT_FAILED -> stop the run (exit 2); resume later
ACC = re.compile(r"^MTBLS\d+$")
META_NAME = re.compile(r"^(i_[^/\\]*\.txt|s_[^/\\]*\.txt|a_[^/\\]*\.txt)$")
FORBIDDEN_NAME = re.compile(r"(^m_|/FILES/|\.tsv$|\.raw$|\.mzml$|\.zip$|\.\.)", re.I)
HREF = re.compile(r'href="([^"?#]+)"', re.I)
BINARY_CT = re.compile(r"zip|gzip|octet-stream|x-tar|compress", re.I)
BINARY_SIG = (b"PK\x03\x04", b"\x1f\x8b", b"PK\x05\x06")
# v0.1.1 (A22, owner): two tiers. STRONG tokens are analytical by nature and remove the row wherever they occur in the
# designated columns. AMBIGUOUS tokens (reference, standard, pool, ...) can name a biological group: they remove the row
# only in a sample-describing column (Sample Name, Material type, Sample type); in a Factor Value / Comment they make the
# study UNRESOLVED (M2 ambiguous QC label), never a silent removal and never a silent inclusion.
QC_STRONG = re.compile(r"(?<![a-z0-9])(qc|qcs|quality control|blank|blanks|pooled qc|qc pool|calibrant|calibration|solvent|wash|spike|spiked|internal standard|system suitability)(?![a-z0-9])")
QC_AMBIG = re.compile(r"(?<![a-z0-9])(reference|standard|standards|std|pool|pooled|pools)(?![a-z0-9])")
SAMPLE_DESCRIBING = ("Sample Name", "Characteristics[Material type]", "Characteristics[Material Type]", "Characteristics[Sample type]", "Characteristics[Sample Type]")
GENERIC_ID = re.compile(r"^(sample|src|source|s|p|sub|subject|patient|participant|id|individual|donor|ctrl|control|case)?[ _\-.]?\d+$")
# C3 (owner): a constant-prefix SHA-256 is deterministic PSEUDONYMISATION, not anonymisation — a low-entropy identifier
# ('P0001') is recoverable by dictionary. Hashes therefore live only in the checkpoint (a working file, listed as
# NOT_FOR_PUBLICATION in the summary); no published table carries them.
D3_SALT = "stageA-D3-v1|"
NOT_FOR_PUBLICATION = ["stageA_census_checkpoint.jsonl"]
MS_NMR = re.compile(r"mass spectrometry|nmr", re.I)
POOL_SPLIT = re.compile(r"(?<![a-z0-9])(pool|pooled|pooling|split|splitting|aliquot)(?![a-z0-9])")
HUMAN = re.compile(r"^(homo sapiens|human)$")


class Halt(Exception):
    pass


class Refused(Exception):
    pass


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


class Redirect:
    def __init__(self, code, location):
        self.code, self.location = code, location


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def url_for(kind: str, acc: str, name: str = "") -> str:
    if not ACC.match(acc):
        raise Refused(f"REFUSED: {acc!r} is not an accession")
    if kind == "FTPDIR":
        return f"{FTP}/{acc}/"
    if kind == "FETCH":
        if FORBIDDEN_NAME.search(name) or not META_NAME.match(name):
            raise Refused(f"REFUSED: {name!r} is not an i_/s_/a_ metadata file")
        return f"{FTP}/{acc}/{quote(name, safe='')}"
    raise Refused(f"REFUSED: unknown request kind {kind!r}")


class Transport:
    """v0.1.4: every request is journaled to the fetch log AS IT OCCURS (append + flush), with its elapsed seconds
    (monotonic clock), for successes and failures alike — an abrupt termination loses at most the request in flight."""
    def __init__(self, offline: dict | None = None, log=None, pause: float = PAUSE, journal: Path | None = None):
        self.offline, self.log, self.pause = offline, (log if log is not None else []), pause
        self.opener = build_opener(_NoRedirect())
        self._jf = journal.open("a") if journal else None
        self._t0 = None

    def _rec(self, entry):
        url, status, nbytes, note = entry
        el = f"{time.monotonic() - self._t0:.2f}" if self._t0 is not None else "0.00"
        self.log.append((url, status, nbytes, note, el))
        if self._jf:
            self._jf.write(f"{url}\t{status}\t{nbytes}\t{note}\telapsed_s={el}\n"); self._jf.flush()

    def close(self):
        if self._jf:
            self._jf.close(); self._jf = None

    def get(self, kind: str, acc: str, name: str = "") -> bytes:
        url = url_for(kind, acc, name)
        self._t0 = time.monotonic()
        if self.offline is not None:
            v = self.offline.get(url)
            if callable(v):
                v = v()
            if v is None:
                v = HTTPError(url, 404, "mock 404", None, None)
            if isinstance(v, Redirect):
                self._rec((url, v.code, 0, f"REDIRECT-NOT-FOLLOWED location={v.location}")); raise Halt(f"{url}: HTTP {v.code} redirect not followed")
            if isinstance(v, HTTPError):
                self._rec((url, v.code, 0, "")); raise Halt(f"{url}: HTTP {v.code}")
            if isinstance(v, Exception):
                self._rec((url, "ERR", 0, f"TRANSPORT-FAILURE {type(v).__name__}")); raise Halt(f"{url}: transport failure ({type(v).__name__})")
            b = v
        else:
            if self.pause:
                time.sleep(self.pause)
            try:
                with self.opener.open(Request(url, headers={"User-Agent": f"stageA_census/{VERSION}"}), timeout=TIMEOUT) as r:
                    if r.status != 200:
                        raise Halt(f"{url}: HTTP {r.status}")
                    ct = str(r.headers.get("Content-Type", "")); cl = r.headers.get("Content-Length")
                    if kind == "FETCH" and (BINARY_CT.search(ct) or (cl and cl.isdigit() and int(cl) > MAX_BYTES)):
                        self._rec((url, 200, 0, f"BINARY-OR-OVERSIZE-REFUSED-UNREAD content-type={ct[:60]}")); raise Halt(f"{url}: non-text or oversize refused unread")
                    b = r.read(MAX_BYTES + 1)
            except HTTPError as e:
                if e.code in (301, 302, 303, 307, 308):
                    loc = str(e.headers.get("Location", ""))[:200] if e.headers else ""
                    self._rec((url, e.code, 0, f"REDIRECT-NOT-FOLLOWED location={loc}")); raise Halt(f"{url}: HTTP {e.code} redirect not followed")
                self._rec((url, e.code, 0, "")); raise Halt(f"{url}: HTTP {e.code}")
            except (URLError, http.client.HTTPException, socket.timeout, OSError) as e:
                self._rec((url, "ERR", 0, f"TRANSPORT-FAILURE {type(e).__name__}: {str(e)[:120]}")); raise Halt(f"{url}: transport failure ({type(e).__name__})")
        if len(b) > MAX_BYTES:
            raise Halt(f"{url}: exceeds {MAX_BYTES} bytes")
        if kind == "FETCH" and any(b.startswith(s) for s in BINARY_SIG):
            self._rec((url, 200, len(b), "BINARY-SIGNATURE-REFUSED-NOT-PARSED")); raise Halt(f"{url}: binary body refused")
        self._rec((url, 200, len(b), sha_bytes(b)))
        return b


# ---------------------------------------------------------------- mirror index
def dir_index(t: Transport, acc: str) -> tuple[str, list[str]]:
    try:
        body = t.get("FTPDIR", acc)
    except Halt as e:
        m = re.search(r"HTTP (\d{3})", str(e))
        if "redirect" in str(e):
            return "INDEX_REDIRECT", []
        if "transport" in str(e):
            return "TRANSPORT_FAILED", []
        return {"404": "INDEX_404", "403": "INDEX_403"}.get(m.group(1) if m else "", "INDEX_OTHER"), []
    html = body.decode("utf-8", "replace")
    # v0.1.7 (D-full-2): index hrefs are percent-encoded ('metabolite%20profiling'); names are DECODED here and encoded
    # exactly once in url_for. v0.1.6 encoded the already-encoded href again ('%2520') -> HTTP 404 on 9 studies.
    names = sorted({unquote(h.rsplit("/", 1)[-1]) for h in HREF.findall(html) if not h.endswith("/") and not h.startswith("/")})
    if not names:
        return "INDEX_UNRECOGNISED", []
    return "INDEX_ACCESSIBLE", names


# ---------------------------------------------------------------- ISA parsing helpers (positional: duplicate headers are real)
def table_rows(text: str) -> tuple[list[str], list[list[str]]]:
    """Positional re-split of a table that stagea_isa has already validated (its own rows are dicts, which collapse the
    duplicate 'Term Source REF' / 'Term Accession Number' headers; the census needs the positions)."""
    lines = [ln for ln in text.splitlines() if ln.strip()]
    cols = [c.strip().strip('"') for c in lines[0].split("\t")]
    return cols, [[c.strip().strip('"') for c in ln.split("\t")] for ln in lines[1:]]


def pairs(cols: list[str], row: list[str]) -> list[tuple[str, str]]:
    return list(zip(cols, row))


def first(pr: list[tuple[str, str]], col: str) -> str | None:
    for c, v in pr:
        if c == col:
            return v
    return None


def accession_after(pr: list[tuple[str, str]], col: str) -> str:
    for i, (c, v) in enumerate(pr):
        if c == col:
            for c2, v2 in pr[i + 1:i + 3]:
                if c2 == "Term Accession Number":
                    return v2
            return ""
    return ""


def row_dict(pr: list[tuple[str, str]]) -> dict:
    """First occurrence wins; keys normalised (what stagea_isa.subject_evidence expects)."""
    d = {}
    for c, v in pr:
        k = isa.norm(c)
        if k not in d:
            d[k] = v
    return d


def inv_fields(v: dict) -> dict:
    f = v["fields"]
    col = lambda label: (f.get(label) or [[]])[0]
    assay_files = col("Study Assay File Name"); techs = col("Study Assay Technology Type"); meas = col("Study Assay Measurement Type")
    return {"title": (col("Study Title") or [""])[0], "description": (col("Study Description") or [""])[0],
            "design_types": [x for x in col("Study Design Type") if x], "factor_names": [x for x in col("Study Factor Name") if x],
            "assays": [{"file": a, "technology": techs[i] if i < len(techs) else "", "measurement": meas[i] if i < len(meas) else ""} for i, a in enumerate(assay_files) if a],
            "protocol_names": [x for x in col("Study Protocol Name") if x], "protocol_types": [x for x in col("Study Protocol Type") if x],
            "affiliations": [x for x in col("Study Person Affiliation") if x], "submission_date": (col("Study Submission Date") or [""])[0],
            "n_labels": len(f)}


# ---------------------------------------------------------------- qualification (§3.5; first failure recorded)
ORGANISM_QC = re.compile(r"(?<![a-z0-9])(reference compound|reference material|reference standard|blank|qc|quality control|standard mixture|calibrant)(?![a-z0-9])")
TAXON_ACC = re.compile(r"ncbitaxon|taxonomy|txid|obo/NCBITaxon", re.I)


def organism_class(pr: list[tuple[str, str]]) -> str:
    """A30 (owner, 2026-10-27, evidentiary clarification of v0.1.6). Classes of one row's Characteristics[Organism]:
       HUMAN                Homo sapiens / human, or a 9606 accession beside a non-empty organism name
       ANALYTICAL           the cell EXPLICITLY names an analytical reference (reference compound/material/standard, blank,
                            qc, quality control, standard mixture, calibrant) and no other metadata conflicts with that
                            reading: the row is QC. The evidence is the explicit designation itself. A missing taxonomy
                            accession is NOT corroboration; it is merely the absence of a conflicting annotation.
       ANALYTICAL_CONFLICT  an explicit analytical label that carries a taxonomy accession (conflicting annotation) and
                            no other QC evidence in the designated columns: study UNRESOLVED, never QC by default
       MISSING              organism cell empty or sentinel: study UNRESOLVED unless other QC evidence removes the row
       NONHUMAN             any other explicit organism: EXCLUDED under M1
    The cell's text takes precedence over a (possibly template-copied) accession. Raw values are preserved."""
    o = isa.norm(first(pr, "Characteristics[Organism]") or "")
    acc = accession_after(pr, "Characteristics[Organism]")
    if ORGANISM_QC.search(o):
        other_qc = any((c in SAMPLE_DESCRIBING or c.startswith("Factor Value[") or c.startswith("Comment[")) and (QC_STRONG.search(isa.norm(v)) or QC_AMBIG.search(isa.norm(v))) for c, v in pr)
        conflicting_taxonomy = bool(TAXON_ACC.search(acc))
        return "ANALYTICAL" if (not conflicting_taxonomy or other_qc) else "ANALYTICAL_CONFLICT"
    if isa.is_missing(o):
        return "MISSING"                           # an accession without an organism name is not accepted as human
    if HUMAN.match(o) or "9606" in acc:
        return "HUMAN"
    return "NONHUMAN"


def qc_class(pr: list[tuple[str, str]]) -> str:
    """'QC' (remove), 'AMBIGUOUS' (study -> UNRESOLVED), or '' (biological)."""
    amb = False
    for c, v in pr:
        if c in SAMPLE_DESCRIBING or c.startswith("Factor Value[") or c.startswith("Comment["):
            n = isa.norm(v)
            if QC_STRONG.search(n):
                return "QC"
            if QC_AMBIG.search(n):
                if c in SAMPLE_DESCRIBING:
                    return "QC"
                amb = True
    if organism_class(pr) == "ANALYTICAL":                # A30: explicit analytical designation, no conflicting annotation
        return "QC"
    return "AMBIGUOUS" if amb else ""


def id_hash(v: str) -> str:
    return hashlib.sha256((D3_SALT + isa.norm(v)).encode()).hexdigest()[:32]


def reconcile_assays(inv_assays: list[dict], mirror_a: list[str]) -> tuple[dict, str]:
    """v0.1.1 (owner's blocking finding): investigation-listed assay names and mirror assay files may differ in name.
    Deterministic rule: exact name match; else, if exactly one investigation assay and exactly one mirror a_ remain
    unmatched, pair them (recorded as 'unique-unmatched pairing'); anything else stays unmatched with technology ''.
    Returns ({mirror_name: {technology, source, inv_name}}, note)."""
    out, note = {}, []
    inv_by = {a["file"]: a for a in inv_assays}
    for n in mirror_a:
        if n in inv_by:
            out[n] = {"technology": inv_by[n]["technology"], "source": "investigation (exact name)", "inv_name": n}
    un_m = [n for n in mirror_a if n not in out]
    un_i = [a for a in inv_assays if a["file"] not in out]
    if len(un_m) == 1 and len(un_i) == 1:
        out[un_m[0]] = {"technology": un_i[0]["technology"], "source": "unique-unmatched pairing", "inv_name": un_i[0]["file"]}
        note.append(f"paired mirror {un_m[0]} with investigation {un_i[0]['file']}")
    for n in un_m:
        if n not in out:
            out[n] = {"technology": "", "source": "unmatched", "inv_name": ""}
            note.append(f"mirror {n} unmatched")
    for a in un_i:
        if a["file"] not in {o["inv_name"] for o in out.values()}:
            note.append(f"investigation {a['file']} absent from mirror")
    return out, "; ".join(note)


def qualify(acc: str, inv: dict, sample: dict, assays: list[dict]) -> dict:
    """inv: inv_fields(); sample: validated s_ (columns, rows); assays: [{file, technology, columns, rows}] validated."""
    rec = {"accession": acc, "title": inv["title"][:200], "design_types": "|".join(inv["design_types"]), "factor_names": "|".join(inv["factor_names"]),
           "status": "UNRESOLVED", "rule": "", "reason": "", "assay_file": "", "platform_class": "", "n_assay_rows": 0,
           "n_sample_rows": sample["n_rows"], "n_qc_rows": 0, "n_biological_rows": 0, "n_sources": 0, "matrix": "", "organism_values": "",
           "design_evidence": "", "subject_evidence": "", "subject_grade": "", "temporal": "", "L1": "NOT_ESTABLISHED", "L1_strict": 0, "L1_sensitivity": 0,
           "feature_count": "FEATURE_COUNT_UNKNOWN", "L2_preliminary": "UNRESOLVED", "L2_confirmed_status": "NOT_ASSESSED", "affiliation": (inv["affiliations"] or [""])[0][:120]}
    # §3.4 analysis choice
    cands = [a for a in assays if MS_NMR.search(a["technology"] or "")]
    if not cands:
        if any(not a["technology"] for a in assays):
            rec.update(status="UNRESOLVED_STRUCTURE", rule="3.4", reason="assay sheet(s) with no reconcilable technology: " + "|".join(a["file"] for a in assays if not a["technology"])[:120]); return rec
        rec.update(status="EXCLUDED", rule="3.4", reason="no mass-spectrometry or NMR assay sheet"); return rec
    cands.sort(key=lambda a: (-a["n_rows"], a["file"]))
    A = cands[0]; rec.update(assay_file=A["file"], n_assay_rows=A["n_rows"], platform_class="NMR" if re.search("nmr", A["technology"], re.I) else "MS")
    scols = sample["columns"]
    srows_all = [pairs(scols, r) for r in sample["rows"]]
    srows = [pr for pr in srows_all if any(v.strip() for c, v in pr)]        # (blank lines are already dropped by the validator)
    rec["n_empty_rows"] = len(srows_all) - len(srows)
    # v0.1.3: a row with an empty Sample Name (MTBLS3 pattern: 64 sheet rows, 63 assay rows, one row with no name) is
    # counted separately; it still makes the study UNRESOLVED under M2 (it is neither a QC row nor an assay sample), because
    # the protocol treats an unidentified row as missing information, not as noise to discard.
    rec["n_rows_without_sample_name"] = sum(1 for pr in srows if not (first(pr, "Sample Name") or "").strip())
    rec["numeric_fraction_s"] = sample.get("numeric_fraction", "")
    a_names = {isa.norm(first(pairs(A["columns"], r), "Sample Name") or "") for r in A["rows"]}
    other_names = {isa.norm(first(pairs(a["columns"], r), "Sample Name") or "") for a in assays if a is not A for r in a["rows"]}
    # M2 QC/blank/pool removal; unidentified rows -> UNRESOLVED
    bio, qc, amb = [], 0, 0
    for pr in srows:
        sn = isa.norm(first(pr, "Sample Name") or "")
        qcls = qc_class(pr)
        if qcls == "QC":
            qc += 1; continue
        if qcls == "AMBIGUOUS":
            amb += 1
        if sn in a_names:
            bio.append(pr)
        elif sn in other_names:
            continue
        else:
            n_unid = sum(1 for q in srows if not qc_class(q) and isa.norm(first(q, "Sample Name") or "") not in a_names | other_names)
            rec.update(status="UNRESOLVED", rule="M2", reason=f"{n_unid} sample-sheet row(s) in no assay and not QC-matched", n_qc_rows=qc, n_rows_unidentified=n_unid); return rec
    sample_bio_names = {isa.norm(first(pr, "Sample Name") or "") for pr in bio}
    missing_in_sheet = [n for n in a_names if n not in sample_bio_names and n not in {isa.norm(first(pr, "Sample Name") or "") for pr in srows}]
    rec.update(n_qc_rows=qc, n_biological_rows=len(bio), n_ambiguous_qc_rows=amb)
    if amb:
        rec.update(status="UNRESOLVED", rule="M2", reason=f"{amb} row(s) carry an ambiguous QC-like label (reference/standard/pool) in a Factor Value or Comment"); return rec
    if missing_in_sheet:
        rec.update(status="UNRESOLVED", rule="M2", reason=f"{len(missing_in_sheet)} assay sample(s) absent from the sample sheet"); return rec
    if not bio:
        rec.update(status="UNRESOLVED", rule="M2", reason="no biological rows after QC removal"); return rec
    # M1 human on every biological row (A30: raw values preserved; analytical-reference rows were removed as QC above;
    # missing or unsupported-analytical organism cells make the study UNRESOLVED, never EXCLUDED and never QC by default)
    rec["organism_values_raw"] = "|".join(sorted({isa.norm(first(pr, "Characteristics[Organism]") or "") for pr in srows}))[:160]
    rec["n_rows_organism_qc_removed"] = sum(1 for pr in srows if organism_class(pr) == "ANALYTICAL" and not any(
        (c in SAMPLE_DESCRIBING or c.startswith("Factor Value[") or c.startswith("Comment[")) and QC_STRONG.search(isa.norm(v)) for c, v in pr))
    ocls = [organism_class(pr) for pr in bio]
    rec["organism_values"] = "|".join(sorted({isa.norm(first(pr, "Characteristics[Organism]") or "") for pr in bio}))[:120]
    rec["n_rows_organism_missing"] = ocls.count("MISSING"); rec["n_rows_organism_analytical_conflict"] = ocls.count("ANALYTICAL_CONFLICT")
    if any(c == "NONHUMAN" for c in ocls):
        rec.update(status="EXCLUDED", rule="M1", reason=f"organism not Homo sapiens on every biological row: {rec['organism_values']}"); return rec
    if any(c in ("MISSING", "ANALYTICAL_CONFLICT") for c in ocls):
        rec.update(status="UNRESOLVED", rule="M1", reason=f"organism missing on {ocls.count('MISSING')} row(s) and/or analytical-reference label with a conflicting taxonomy annotation on {ocls.count('ANALYTICAL_CONFLICT')} row(s)"); return rec
    # M3 identity
    for pr in bio:
        if isa.is_missing(first(pr, "Sample Name") or "") or isa.is_missing(first(pr, "Source Name") or ""):
            rec.update(status="UNRESOLVED", rule="M3", reason="missing or sentinel Sample Name / Source Name"); return rec
    names = [isa.norm(first(pr, "Sample Name") or "") for pr in bio]
    if len(set(names)) != len(names):
        rec.update(status="UNRESOLVED", rule="M3", reason="duplicate Sample Name among biological rows"); return rec
    # M4 matrix
    mcol = "Characteristics[Organism part]" if any(c == "Characteristics[Organism part]" for c in scols) else ("Characteristics[Sample type]" if any(c == "Characteristics[Sample type]" for c in scols) else None)
    if mcol is None:
        rec.update(status="UNRESOLVED", rule="M4", reason="no Organism part / Sample type column"); return rec
    mats = [isa.norm(first(pr, mcol) or "") for pr in bio]
    if any(isa.is_missing(m) for m in mats):
        rec.update(status="UNRESOLVED", rule="M4", reason=f"{mcol} missing on a biological row"); return rec
    if len(set(mats)) > 1:
        rec.update(status="EXCLUDED", rule="M4", reason=f"> 1 matrix: {sorted(set(mats))[:4]}", matrix="|".join(sorted(set(mats)))[:80]); return rec
    rec["matrix"] = mats[0][:80]
    # M5 pooling/splitting protocol on the path
    prot_vals = {isa.norm(v) for pr in bio for c, v in pr if c == "Protocol REF"}
    pooled = sorted(p for p in prot_vals if POOL_SPLIT.search(p)) + sorted(p for p in inv["protocol_names"] if POOL_SPLIT.search(isa.norm(p)) and isa.norm(p) in prot_vals)
    if pooled:
        rec.update(status="EXCLUDED", rule="M5", reason=f"pooling/splitting Protocol REF on the path: {pooled[:3]}"); return rec
    # M5 D1 via the frozen recogniser (includes M7 n >= 50 and R1/R2 grades)
    # v0.1.7 (reporting only): person-key-like columns present on the sheet, whatever the D1 path; lets the owner see
    # how often R1 was never evaluated because the Source Name test came first (A27 / multiplicity).
    rec["person_key_columns"] = "|".join(sorted({isa.norm(c) for c in scols if isa.PERSON_KEY.match(isa.norm(c))}))
    ev = isa.subject_evidence([row_dict(pr) for pr in bio], inv["design_types"], inv["description"], min_n=MIN_N)
    rec.update(n_sources=ev["n_sources"], subject_grade=ev["grade"], subject_evidence=ev["reason"][:200],
               design_evidence=(f"layer1 pass={ev['layer1']['pass']} +{ev['layer1']['positive']} -{ev['layer1']['negative']}" if "layer1" in ev else ("layer1 not reached (R1)" if ev["grade"] == "R1" else "layer1 not reached"))[:200])
    if ev["grade"] == "EXCLUDED" and "Source Name carries" in ev["reason"] and ev["n_sources"] == 1 and len(bio) > 1:
        # v0.1.3 (A27, proposed after pilot 25: MTBLS19/25/28 carry ONE Source Name for every sample — the source node is
        # a study-level placeholder, not a subject). A degenerate source node is uninformative: UNRESOLVED, not EXCLUDED.
        # n_sources > 1 with multiplicity stays EXCLUDED (demonstrated > 1 sample per source).
        rec.update(status="UNRESOLVED", rule="M5", reason=f"source node degenerate (one Source Name for {len(bio)} samples): subject identity not evidenced", subject_grade="SOURCE_DEGENERATE"); return rec
    if ev["grade"] == "EXCLUDED":
        rec.update(status="EXCLUDED", rule="M5/M7", reason=ev["reason"]); return rec
    if ev["grade"] in ("UNRESOLVED", "SOURCE_UNIQUE"):
        rec.update(status="UNRESOLVED", rule="M5", reason=ev["reason"]); return rec
    # M6 temporal (two levels)
    fcols = [c for c in scols if c.startswith("Factor Value[")]
    tclasses = {}
    for c in fcols:
        tc = isa.temporal_class(c[13:-1])
        if tc != "none":
            levels = {isa.norm(first(pr, c) or "") for pr in bio}
            tclasses[c] = (tc, len(levels))
    for fn in inv["factor_names"]:
        tc = isa.temporal_class(fn)
        if tc != "none" and f"Factor Value[{fn}]" not in tclasses:
            tclasses[f"factor:{fn}"] = (tc, -1)
    note = ""
    strong = [(c, n) for c, (tc, n) in tclasses.items() if tc == "strong"]
    ambig = [(c, n) for c, (tc, n) in tclasses.items() if tc == "ambiguous"]
    if any(n != 1 for c, n in strong):
        # per-source multiplicity already 1 (subject_evidence) -> independent groups at different times
        note = "CONFIRMED-with-note: strong temporal factor with > 1 level, one sample per source (independent groups)"
    elif ambig:
        if ev["grade"] == "R1" or ev.get("layer1", {}).get("pass"):
            note = "CONFIRMED-with-note: ambiguous temporal factor; non-temporal use established by R1 / Layer 1 design"
        else:
            rec.update(status="UNRESOLVED", rule="M6", reason=f"ambiguous temporal factor without design evidence: {[c for c, n in ambig][:3]}"); return rec
    rec["temporal"] = (note or "none") + (f" {sorted(tclasses)[:3]}" if tclasses else "")
    # M7 inside subject_evidence; M8 A3 -> UNKNOWN; M9 touched
    if acc in EXEMPLARS:
        rec.update(status="EXCLUDED", rule="M9", reason="touched (schema exemplar, A7)"); return rec
    pk = sorted(k for k in row_dict(bio[0]) if isa.PERSON_KEY.match(k))
    ids = {"source": sorted({id_hash(first(pr, "Source Name") or "") for pr in bio}), "sample": sorted({id_hash(first(pr, "Sample Name") or "") for pr in bio}),
           "person": sorted({id_hash(row_dict(pr).get(pk[0], "")) for pr in bio}) if pk else [],
           "specific_source": sorted({id_hash(first(pr, "Source Name") or "") for pr in bio if not GENERIC_ID.match(isa.norm(first(pr, "Source Name") or ""))}),
           "specific_person": sorted({id_hash(row_dict(pr).get(pk[0], "")) for pr in bio if not GENERIC_ID.match(isa.norm(row_dict(pr).get(pk[0], "")))}) if pk else []}
    rec["_ids"] = ids
    rec.update(status="ELIGIBLE", rule="D2", reason=f"L1 {ev['grade']} n={ev['n_sources']}" + ("" if not note else f"; {note}"),
               L1="CONFIRMED" if ev["grade"] in ("R1", "R2_T12") else "CONFIRMED_SENSITIVITY_ONLY",
               L1_strict=int(ev["grade"] in ("R1", "R2_T12")), L1_sensitivity=1, L2_preliminary="UNRESOLVED (FEATURE_COUNT_UNKNOWN)")
    return rec


# ---------------------------------------------------------------- per-study pipeline
def screen_study(t: Transport, acc: str) -> tuple[dict, list[dict]]:
    files: list[dict] = []
    base = {"accession": acc, "status": "", "rule": "", "reason": "", "index_class": "", "n_files_fetched": 0}
    cls, names = dir_index(t, acc)
    base["index_class"] = cls
    if cls != "INDEX_ACCESSIBLE":
        base.update(status=cls, rule="3.3", reason="mirror folder index not usable; not screened"); return base, files
    meta = [n for n in names if META_NAME.match(n) and not FORBIDDEN_NAME.search(n)]
    base["n_forbidden_names_seen_never_fetched"] = sum(1 for n in names if FORBIDDEN_NAME.search(n))
    if not any(n.startswith("i_") for n in meta) or not any(n.startswith("s_") for n in meta):
        base.update(status="FORMAT_REFUSED", rule="3.3", reason="mirror lists no investigation or no sample sheet"); return base, files
    texts = {}
    for n in sorted(meta):
        try:
            b = t.get("FETCH", acc, n)
        except Halt as e:
            m = re.search(r"HTTP (\d{3})$", str(e))
            # v0.1.7 (D-full-3): an HTTP error on a file the index lists is FETCH_HTTP_<code>; ACCESS_REFUSED is kept for
            # what the protocol meant by it (redirect, binary or oversize response, never followed or parsed).
            st = "TRANSPORT_FAILED" if "transport" in str(e) else (f"FETCH_HTTP_{m.group(1)}" if (m and "redirect" not in str(e)) else "ACCESS_REFUSED")
            base.update(status=st, rule="3.3", reason=str(e)[:160], n_files_fetched=len(files)); return base, files
        files.append({"accession": acc, "name": n, "bytes": len(b), "sha256": sha_bytes(b)})
        texts[n] = b.decode("utf-8", "replace")
    base["n_files_fetched"] = len(files)
    invs = [n for n in texts if n.startswith("i_")]
    if len(invs) != 1:
        base.update(status="FORMAT_REFUSED", rule="3.3", reason=f"{len(invs)} investigation files"); return base, files
    v = isa.validate_investigation(texts[invs[0]])
    if not v["ok"]:
        base.update(status="FORMAT_REFUSED", rule="3.3", reason=f"i_: {v['errors'][:2]}"); return base, files
    inv = inv_fields(v)
    ss = [n for n in texts if n.startswith("s_")]
    if len(ss) != 1:
        base.update(status="UNRESOLVED_STRUCTURE", rule="A23", reason=f"{len(ss)} sample sheets (parser precondition: exactly one); not a biological exclusion"); return base, files
    sv = isa.validate_sample_sheet(texts[ss[0]])
    if not sv["ok"]:
        base.update(status="FORMAT_REFUSED", rule="3.3", reason=f"s_: {sv['errors'][:2]}"); return base, files
    sv["columns"], sv["rows"] = table_rows(texts[ss[0]])
    assays = []
    mirror_a = sorted(k for k in texts if k.startswith("a_"))
    recon, recon_note = reconcile_assays(inv["assays"], mirror_a)
    for n in mirror_a:
        av = isa.validate_assay_sheet(texts[n])
        if not av["ok"]:
            base.update(status="FORMAT_REFUSED", rule="3.3", reason=f"{n}: {av['errors'][:2]}"); return base, files
        acols, arows = table_rows(texts[n])
        assays.append({"file": n, "technology": recon[n]["technology"], "technology_source": recon[n]["source"], "investigation_name": recon[n]["inv_name"],
                       "columns": acols, "rows": arows, "n_rows": av["n_rows"]})
    rec = qualify(acc, inv, sv, assays)
    chosen = next((a for a in assays if a["file"] == rec.get("assay_file")), None)
    rec.update(index_class=cls, n_files_fetched=len(files), assay_reconciliation=recon_note[:200],
               assay_technology_source=(chosen or {}).get("technology_source", ""), assay_investigation_name=(chosen or {}).get("investigation_name", ""),
               n_forbidden_names_seen_never_fetched=base["n_forbidden_names_seen_never_fetched"])
    return rec, files


# ---------------------------------------------------------------- D3 units and strata (§5)
def units(cohorts: list[dict]) -> tuple[list[dict], list[dict]]:
    """D3 (§5) with the owner's safeguards (A24/C2): CONFIRMED merge ONLY on overlap of SPECIFIC person-key hashes with
    the same affiliation. Same affiliation + title-token Jaccard >= 0.6, and any other identifier-hash overlap (generic
    ids, sample names, different affiliations), are SUSPECTED_OVERLAP review flags, never merges."""
    el = [c for c in cohorts if c["status"] == "ELIGIBLE"]
    flags = []
    parent = {c["accession"]: c["accession"] for c in el}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def tokens(s):
        return {w for w in re.findall(r"[a-z0-9]+", isa.norm(s)) if len(w) > 3}
    reasons = {}
    for i, a in enumerate(el):
        for b in el[i + 1:]:
            why = None
            ia, ib = a.get("_ids") or {}, b.get("_ids") or {}
            same_aff = bool(a["affiliation"]) and isa.norm(a["affiliation"]) == isa.norm(b["affiliation"])
            sp = set(ia.get("specific_person", [])) & set(ib.get("specific_person", []))
            if sp and same_aff:
                why = f"specific person-identifier overlap ({len(sp)}) with the same affiliation"
            if why:
                parent[find(a["accession"])] = find(b["accession"]); reasons[(a["accession"], b["accession"])] = why
            else:
                if same_aff:
                    ta, tb = tokens(a["title"]), tokens(b["title"])
                    if ta and tb and len(ta & tb) / len(ta | tb) >= 0.6:
                        flags.append({"accession_a": a["accession"], "accession_b": b["accession"], "identifier_kind": "title", "n_overlap": len(ta & tb),
                                      "n_specific_overlap": 0, "same_affiliation": 1, "status": "SUSPECTED_OVERLAP (same affiliation, title similarity >= 0.6; review flag; not merged)"})
                for kind in ("person", "source", "sample"):
                    ov = set(ia.get(kind, [])) & set(ib.get(kind, []))
                    if ov:
                        spec = len(set(ia.get(f"specific_{kind}", [])) & set(ib.get(f"specific_{kind}", []))) if kind != "sample" else 0
                        flags.append({"accession_a": a["accession"], "accession_b": b["accession"], "identifier_kind": kind, "n_overlap": len(ov),
                                      "n_specific_overlap": spec, "same_affiliation": int(same_aff), "status": "SUSPECTED_OVERLAP (review flag; not merged)"}); break
    groups: dict[str, list[dict]] = {}
    for c in el:
        groups.setdefault(find(c["accession"]), []).append(c)
    out = []
    for root, members in sorted(groups.items()):
        strata = {(m["platform_class"], m["matrix"]) for m in members}
        strict = any(m["L1_strict"] for m in members)
        out.append({"unit": root, "accessions": "|".join(sorted(m["accession"] for m in members)), "n_accessions": len(members),
                    "repository": "MetaboLights", "modality": "metabolomics", "platform_class": "|".join(sorted({p for p, _ in strata})),
                    "matrix": "|".join(sorted({x for _, x in strata}))[:80], "L1_strict": int(strict), "L1_sensitivity": 1,
                    "merge_reason": "; ".join(sorted({w for (x, y), w in reasons.items() if find(x) == root}))[:160]})
    return out, flags


# ---------------------------------------------------------------- manifest / lock / run
def load_universe(path: Path, expected_sha: str) -> list[str]:
    b = path.read_bytes()
    if sha_bytes(b) != expected_sha:
        raise Refused(f"REFUSED: manifest SHA-256 {sha_bytes(b)[:12]}… != authorised {expected_sha[:12]}…")
    rows = list(csv.DictReader(b.decode().splitlines()))
    if any(r["status"].startswith("PROVISIONAL") for r in rows):
        raise Refused("REFUSED: provisional manifest (A17 verdict was UNRESOLVED)")
    accs = [r["accession"] for r in rows if r["status"] == "UNIVERSE" and ACC.match(r["accession"])]
    return [a for a in sorted(accs, key=lambda a: int(a[5:])) if a not in EXEMPLARS]


def auth_status(p: Path | None) -> dict:
    if not p or not p.exists():
        return {"present": False, "valid": False}
    m = AUTH_LINE.fullmatch(p.read_text().strip())
    return {"present": True, "valid": bool(m), "line": p.read_text().strip()[:220], "date": m.group(1) if m else None,
            "universe_sha256": m.group(2) if m else None, "protocol": m.group(3) if m else None, "owner": m.group(4) if m else None,
            "pilot": int(m.group(5)) if m and m.group(5) else None}


def sha_file(p: Path | None):
    return sha_bytes(p.read_bytes()) if p and p.exists() else None


def lock_manifest(protocol: Path | None, auth: Path | None, manifest: Path | None) -> dict:
    me = Path(__file__).resolve()
    return {"version": VERSION, "protocol": PROTOCOL, "isa_version": isa.VERSION,
            "requests_per_study": ["FTPDIR names-only index", "FETCH of every listed i_/s_/a_ (guards: name class, host, redirect, size, binary)"],
            "never": ["REST", "m_ files", "FILES/", "raw/derived", "any value", "redirect destinations", "retries within a run"],
            "rules": "protocol §2 D1–D4, §3.3–3.6, §5, §6; stagea_isa frozen recogniser; MIN_N 50; exemplars MTBLS1, MTBLS2 excluded (A7)",
            "failure_policy": "per-study TRANSPORT_FAILED / FORMAT_REFUSED / INDEX_*; census continues; separate disclosed --retry-pass only",
            "pause_s": PAUSE, "hashes": {me.name: sha_file(me), "stagea_isa.py": sha_file(here / "stagea_isa.py"), PROTOCOL: sha_file(protocol),
                                         AUTH_FILE: sha_file(auth), MANIFEST: sha_file(manifest)},
            "authorisation_file": auth_status(auth), "owner_signoff": None, "created_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}


def cmd_plan(a):
    m = lock_manifest(Path(a.protocol), Path(a.auth), Path(a.manifest))
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    (out / "stageA_census_lock.json").write_text(json.dumps(m, indent=2) + "\n"); print(json.dumps(m, indent=2)); print("[plan] written; owner adds owner_signoff before any run")


def retry_eligible(r: dict) -> bool:
    """v0.1.7: the disclosed retry pass re-requests (i) TRANSPORT_FAILED, (ii) FETCH_HTTP_* and (iii) legacy v0.1.6
    ACCESS_REFUSED records whose reason shows the double-encoding defect ('%25', not a redirect/binary/oversize refusal).
    Nothing else is revisited."""
    s_, why = r.get("status", ""), r.get("reason", "")
    # the v0.1.6 reason is truncated at 160 characters, so 'HTTP 404' may be cut off: the signature is '%25' (a
    # double-encoded '%') in an ACCESS_REFUSED reason that is not a redirect / binary / oversize refusal
    return s_ == "TRANSPORT_FAILED" or s_.startswith("FETCH_HTTP_") or (
        s_ == "ACCESS_REFUSED" and "%25" in why and not any(k in why for k in ("redirect", "binary", "oversize", "non-text")))


def run_census(t: Transport, accs: list[str], out: Path, retry_pass: bool = False, progress=None) -> dict:
    out.mkdir(parents=True, exist_ok=True)
    ck = out / "stageA_census_checkpoint.jsonl"
    done: dict[str, dict] = {}
    if ck.exists():
        for ln in ck.read_text().splitlines():
            if ln.strip():
                r = json.loads(ln); done[r["accession"]] = r
    todo = [a for a in accs if a not in done or (retry_pass and retry_eligible(done[a]))]
    jpath = out / "stageA_census_fetch_log.txt"
    with jpath.open("a") as jf:
        jf.write(f"# run {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} census v{VERSION} todo={len(todo)} retry_pass={int(retry_pass)}\n")
    if t._jf is None:
        t._jf = jpath.open("a")
    files_path = out / "stageA_M_files.csv"
    sample_ids: dict[str, set] = {}
    with ck.open("a") as ckf, files_path.open("a", newline="") as ff:
        fw = csv.writer(ff)
        if files_path.stat().st_size == 0:
            fw.writerow(["accession", "name", "bytes", "sha256"])
        streak = 0
        for i, acc in enumerate(todo):
            try:
                rec, files = screen_study(t, acc)
            except Refused as e:                      # A20: a guard refusal is a protocol violation -> global halt
                summ = finalize([done[a] for a in accs if a in done], out, t.log)
                summ["halted_protocol_violation"] = f"{acc}: {e}"; (out / "stageA_M_summary.json").write_text(json.dumps(summ, indent=2) + "\n")
                return summ
            rec["run_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()); rec["census_version"] = VERSION
            for f in files:
                fw.writerow([f["accession"], f["name"], f["bytes"], f["sha256"]])
            ckf.write(json.dumps(rec) + "\n"); ckf.flush(); ff.flush()
            done[acc] = rec
            if progress and (i % 25 == 0 or i == len(todo) - 1):
                progress(i + 1, len(todo), acc, rec["status"])
            streak = streak + 1 if rec["status"] == "TRANSPORT_FAILED" else 0
            if streak >= OUTAGE_K:                    # D-full-1: a network outage must not "screen" the rest of the list
                t.close()
                summ = finalize([done[a] for a in accs if a in done], out, t.log)
                summ["halted_outage"] = f"{OUTAGE_K} consecutive studies ended TRANSPORT_FAILED (last {acc}); run stopped; resume with the same command once the network is back"
                (out / "stageA_M_summary.json").write_text(json.dumps(summ, indent=2) + "\n")
                return summ
    recs = [done[a] for a in accs if a in done]
    t.close()
    return finalize(recs, out, t.log)


FIELDS = ["accession", "status", "rule", "reason", "index_class", "n_files_fetched", "assay_file", "assay_technology_source", "assay_investigation_name", "assay_reconciliation",
          "platform_class", "n_assay_rows", "n_sample_rows", "n_empty_rows", "n_rows_without_sample_name", "n_rows_unidentified", "numeric_fraction_s", "n_qc_rows", "n_ambiguous_qc_rows",
          "n_biological_rows", "n_sources", "organism_values", "organism_values_raw", "n_rows_organism_qc_removed", "n_rows_organism_missing",
          "n_rows_organism_analytical_conflict", "matrix", "design_evidence", "subject_evidence", "subject_grade", "temporal", "L1", "L1_strict",
          "L1_sensitivity", "feature_count", "L2_preliminary", "L2_confirmed_status", "title", "design_types", "factor_names", "affiliation",
          "n_forbidden_names_seen_never_fetched", "person_key_columns", "census_version", "run_utc"]


def finalize(recs: list[dict], out: Path, log: list) -> dict:
    with (out / "stageA_M_cohorts.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore"); w.writeheader()
        for r in recs:
            w.writerow({k: r.get(k, "") for k in FIELDS})
    with (out / "stageA_M_exclusions.csv").open("w", newline="") as f:
        w = csv.writer(f); w.writerow(["accession", "status", "rule", "reason"])
        for r in recs:
            if r["status"] != "ELIGIBLE":
                w.writerow([r["accession"], r["status"], r.get("rule", ""), r.get("reason", "")])
    us, flags = units(recs)
    with (out / "stageA_M_overlap_flags.csv").open("w", newline="") as f:
        cols = ["accession_a", "accession_b", "identifier_kind", "n_overlap", "n_specific_overlap", "same_affiliation", "status"]
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); [w.writerow(x) for x in flags]
    with (out / "stageA_M_units.csv").open("w", newline="") as f:
        cols = ["unit", "accessions", "n_accessions", "repository", "modality", "platform_class", "matrix", "L1_strict", "L1_sensitivity", "merge_reason"]
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); [w.writerow(u) for u in us]
    by = lambda key: {k: sum(1 for r in recs if key(r) == k) for k in sorted({key(r) for r in recs})}
    strata = {}
    for u in us:
        k = f"{u['repository']}×{u['modality']}×{u['platform_class']}×{u['matrix']}"
        strata.setdefault(k, {"units": 0, "units_strict": 0}); strata[k]["units"] += 1; strata[k]["units_strict"] += u["L1_strict"]
    summ = {"version": VERSION, "isa_version": isa.VERSION, "n_screened": len(recs), "status_counts": by(lambda r: r["status"]),
            "rule_counts": by(lambda r: f"{r['status']}:{r.get('rule', '')}"),
            "L1_strict": sum(r.get("L1_strict", 0) for r in recs), "L1_sensitivity": sum(r.get("L1_sensitivity", 0) for r in recs),
            "L1_unresolved": sum(1 for r in recs if r["status"] == "UNRESOLVED"),
            "subject_grades_among_screened": by(lambda r: r.get("subject_grade", "") or "-"),
            "L2_preliminary": {"confirmed": 0, "unresolved_FEATURE_COUNT_UNKNOWN": sum(1 for r in recs if r["status"] == "ELIGIBLE"), "excluded": 0},
            "L2_confirmed_status": "NOT_ASSESSED", "D3_metadata_units": len(us), "D3_metadata_units_strict": sum(u["L1_strict"] for u in us),
            "D3_metadata_units_by_stratum": strata, "L3_confirmed_status": "NOT_ASSESSED",
            "artefacts_not_for_publication": NOT_FOR_PUBLICATION, "identifier_hashing": "deterministic pseudonymisation (constant prefix); not anonymisation; hashes only in the checkpoint",
            "transport_failed": sum(1 for r in recs if r["status"] == "TRANSPORT_FAILED"),
            "structure_unresolved": sum(1 for r in recs if r["status"] == "UNRESOLVED_STRUCTURE"), "access_refused": sum(1 for r in recs if r["status"] == "ACCESS_REFUSED"),
            "format_refused": sum(1 for r in recs if r["status"] == "FORMAT_REFUSED"), "D3_suspected_overlap_flags": len(flags), "n_requests_this_invocation": len(log),
            "run_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    # v0.1.4: the fetch log is a journal written by the transport as requests occur (run header at invocation start);
    # finalize only counts it. (v0.1.3 appended it here, after the fact — fine for a completed run, lost on interruption.)
    jp = out / "stageA_census_fetch_log.txt"
    summ["n_requests_cumulative"] = sum(1 for ln in jp.read_text().splitlines() if ln and not ln.startswith("#")) if jp.exists() else 0
    (out / "stageA_M_summary.json").write_text(json.dumps(summ, indent=2) + "\n")
    names = ["stageA_M_cohorts.csv", "stageA_M_exclusions.csv", "stageA_M_units.csv", "stageA_M_overlap_flags.csv", "stageA_M_summary.json", "stageA_M_files.csv", "stageA_census_fetch_log.txt", "stageA_census_checkpoint.jsonl"]
    (out / "SHA256SUMS.txt").write_text("".join(f"{sha_bytes((out / n).read_bytes())}  {n}\n" for n in names if (out / n).exists()))
    return summ


def cmd_run(a):
    out = Path(a.out); lock = out / "stageA_census_lock.json"
    if not lock.exists():
        sys.exit("[run] REFUSED: no stageA_census_lock.json")
    m = json.loads(lock.read_text())
    if not m.get("owner_signoff"):
        sys.exit("[run] REFUSED: lock not signed by owner")
    now = lock_manifest(Path(a.protocol), Path(a.auth), Path(a.manifest))["hashes"]
    if any(v is None for v in now.values()) or now != m.get("hashes"):
        sys.exit(f"[run] REFUSED: hashes differ from the lock or a locked file is missing\n locked={m.get('hashes')}\n now={now}")
    st = auth_status(Path(a.auth))
    if not st["valid"] or st["protocol"] != PROTOCOL or Path(a.protocol).name != PROTOCOL or st["universe_sha256"] != now[MANIFEST]:
        sys.exit(f"[run] REFUSED: authorisation invalid, not naming {PROTOCOL}, or its universe hash differs from the manifest present")
    try:
        accs = load_universe(Path(a.manifest), st["universe_sha256"])
    except Refused as e:
        sys.exit(f"[run] {e}")
    # A21: a pilot is authorised by `pilot=N` in the authorisation line; selection rule fixed: the first N accessions of the
    # manifest in numeric order. A full run is refused while the line carries pilot=N; a --pilot flag is refused unless it
    # equals the authorised N.
    if st["pilot"] is not None:
        if a.pilot != st["pilot"]:
            sys.exit(f"[run] REFUSED: authorisation is for pilot={st['pilot']}; run with --pilot {st['pilot']}")
        accs = accs[:st["pilot"]]
    elif a.pilot:
        sys.exit("[run] REFUSED: --pilot given but the authorisation line carries no pilot=N")
    t = Transport()
    summ = run_census(t, accs, out, retry_pass=a.retry_pass, progress=lambda i, n, acc, s: print(f"[{i}/{n}] {acc} {s}", flush=True))
    finish_run(summ)


def run_exit_code(summ: dict) -> int:
    """C1 (owner): a protocol violation is a failed run (1). v0.1.7: an outage stop is 2. 0 only when the list was finished."""
    return 1 if summ.get("halted_protocol_violation") else (2 if summ.get("halted_outage") else 0)


def finish_run(summ: dict):
    print(json.dumps(summ, indent=2)); sys.exit(run_exit_code(summ))


# ---------------------------------------------------------------- selftest (synthetic studies; offline)
def _inv(acc, title, desc, design, factors, assays, prot=("Sample collection", "Extraction"), affil="Inst A"):
    lines = ["ONTOLOGY SOURCE REFERENCE", "Term Source Name\tOBI", "INVESTIGATION", f"Investigation Identifier\t{acc}", "STUDY", f"Study Identifier\t{acc}",
             f"Study Title\t{title}", f"Study Description\t{desc}", f"Study File Name\ts_{acc}.txt", "STUDY DESIGN DESCRIPTORS",
             "Study Design Type\t" + "\t".join(design), "STUDY FACTORS", "Study Factor Name\t" + "\t".join(factors), "STUDY ASSAYS",
             "Study Assay File Name\t" + "\t".join(f for f, _ in assays), "Study Assay Measurement Type\t" + "\t".join("metabolite profiling" for _ in assays),
             "Study Assay Technology Type\t" + "\t".join(t for _, t in assays), "STUDY PROTOCOLS", "Study Protocol Name\t" + "\t".join(prot),
             "Study Protocol Type\t" + "\t".join(prot), "STUDY CONTACTS", "Study Person Last Name\tX", f"Study Person Affiliation\t{affil}"]
    return "\n".join(lines) + "\n"


def _sheet(n, organism="Homo sapiens", acc="NCBITAXON:9606", part="urine", person=False, qc=0, src_fn=None, extra_factor=None, parts=None, protocol="Sample collection"):
    cols = ["Source Name", "Characteristics[Organism]", "Term Source REF", "Term Accession Number", "Characteristics[Organism part]"] + (["Characteristics[Subject ID]"] if person else []) + \
           ["Protocol REF", "Sample Name", "Factor Value[Group]"] + ([f"Factor Value[{extra_factor}]"] if extra_factor else [])
    rows = []
    for i in range(n):
        src = (src_fn(i) if src_fn else f"SRC{i:04d}")
        rows.append([src, organism, "NCBITAXON", acc, (parts[i % len(parts)] if parts else part)] + ([f"P{i:04d}"] if person else []) + [protocol, f"S{i:04d}", "case" if i % 2 else "control"] + ([f"visit {i % 2 + 1}"] if extra_factor else []))
    for j in range(qc):
        rows.append([f"QC{j}", organism, "NCBITAXON", acc, part] + (["QC"] if person else []) + [protocol, f"QC_sample_{j}", "QC"] + (["-"] if extra_factor else []))
    return "\n".join(["\t".join(cols)] + ["\t".join(r) for r in rows]) + "\n"


def _assay(n, qc=0, tech="nmr"):
    cols = ["Sample Name", "Protocol REF", "NMR Assay Name" if tech == "nmr" else "MS Assay Name", "Metabolite Assignment File"]
    rows = [[f"S{i:04d}", "Extraction", f"A{i}", "m_x.tsv"] for i in range(n)] + [[f"QC_sample_{j}", "Extraction", f"Q{j}", "m_x.tsv"] for j in range(qc)]
    return "\n".join(["\t".join(cols)] + ["\t".join(r) for r in rows]) + "\n"


def _study(off, acc, inv, sheet, assays: dict, extra_names=("m_x.tsv",)):
    names = ["i_Investigation.txt", f"s_{acc}.txt"] + list(assays) + list(extra_names)
    off[url_for("FTPDIR", acc)] = ("<html>" + "".join(f'<a href="{n}">{n}</a>' for n in names) + '<a href="/pub/">up</a></html>').encode()
    off[url_for("FETCH", acc, "i_Investigation.txt")] = inv.encode(); off[url_for("FETCH", acc, f"s_{acc}.txt")] = sheet.encode()
    for n, txt in assays.items():
        off[url_for("FETCH", acc, n)] = txt.encode()


def _mock():
    off = {}
    cs = ["cross sectional"]
    _study(off, "MTBLS10", _inv("MTBLS10", "Urine NMR cohort A", "Urine was collected from each participant.", cs, ["Group"], [("a_MTBLS10_nmr.txt", "NMR spectroscopy")]),
           _sheet(60, person=True, qc=3), {"a_MTBLS10_nmr.txt": _assay(60, qc=3)})                                   # R1, n=60
    _study(off, "MTBLS11", _inv("MTBLS11", "Plasma MS cohort B", "One blood sample was collected from each participant at enrolment.", cs, ["Group"], [("a_MTBLS11_ms.txt", "mass spectrometry")]),
           _sheet(55, part="plasma"), {"a_MTBLS11_ms.txt": _assay(55, tech="ms")})                                     # R2_T12, n=55
    _study(off, "MTBLS12", _inv("MTBLS12", "Serum MS cohort C", "Samples were collected from 52 unique individuals.", cs, ["Group"], [("a_MTBLS12_ms.txt", "mass spectrometry")]),
           _sheet(52, part="serum"), {"a_MTBLS12_ms.txt": _assay(52, tech="ms")})                                      # R2_COUNT_SUPPORTED, n=52
    _study(off, "MTBLS13", _inv("MTBLS13", "Urine NMR cohort D", "Participants were recruited in two hospitals.", cs, ["Group"], [("a_MTBLS13_nmr.txt", "NMR spectroscopy")]),
           _sheet(70), {"a_MTBLS13_nmr.txt": _assay(70)})                                                             # SOURCE_UNIQUE -> UNRESOLVED
    _study(off, "MTBLS14", _inv("MTBLS14", "Arabidopsis leaves", "Leaves were harvested.", ["Plant Metabolomics"], ["genotype"], [("a_MTBLS14_ms.txt", "mass spectrometry")]),
           _sheet(16, organism="Arabidopsis thaliana", acc="NCBITAXON:3702", part="leaf"), {"a_MTBLS14_ms.txt": _assay(16, tech="ms")})  # M1 EXCLUDED
    _study(off, "MTBLS15", _inv("MTBLS15", "Repeated plasma", "One sample was collected from each participant.", cs, ["Group"], [("a_MTBLS15_ms.txt", "mass spectrometry")]),
           _sheet(60, person=True, src_fn=lambda i: f"SRC{i // 2:04d}"), {"a_MTBLS15_ms.txt": _assay(60, tech="ms")})  # source with 2 samples -> EXCLUDED
    _study(off, "MTBLS16", _inv("MTBLS16", "Two matrices", "One sample was collected from each participant.", cs, ["Group"], [("a_MTBLS16_ms.txt", "mass spectrometry")]),
           _sheet(60, person=True, parts=["plasma", "urine"]), {"a_MTBLS16_ms.txt": _assay(60, tech="ms")})             # M4 EXCLUDED
    _study(off, "MTBLS17", _inv("MTBLS17", "Small", "One sample was collected from each participant.", cs, ["Group"], [("a_MTBLS17_ms.txt", "mass spectrometry")]),
           _sheet(49, person=True), {"a_MTBLS17_ms.txt": _assay(49, tech="ms")})                                       # n=49 EXCLUDED
    _study(off, "MTBLS18", _inv("MTBLS18", "Longitudinal", "One sample was collected from each participant.", ["longitudinal", "cross sectional"], ["Group"], [("a_MTBLS18_ms.txt", "mass spectrometry")]),
           _sheet(60), {"a_MTBLS18_ms.txt": _assay(60, tech="ms")})                                                    # Layer 1 negative -> UNRESOLVED
    _study(off, "MTBLS19", _inv("MTBLS19", "Visits", "Urine was collected from each participant.", cs, ["Group", "Visit"], [("a_MTBLS19_nmr.txt", "NMR spectroscopy")]),
           _sheet(60, person=True, extra_factor="Visit"), {"a_MTBLS19_nmr.txt": _assay(60)})                           # strong temporal, 1 sample/source -> with-note
    _study(off, "MTBLS20", _inv("MTBLS20", "Pooled", "One sample was collected from each participant.", cs, ["Group"], [("a_MTBLS20_ms.txt", "mass spectrometry")], prot=("Pooling", "Extraction")),
           _sheet(60, person=True, protocol="Pooling"), {"a_MTBLS20_ms.txt": _assay(60, tech="ms")})                   # pooling protocol -> EXCLUDED
    _study(off, "MTBLS21", _inv("MTBLS21", "No MS", "x", cs, ["Group"], [("a_MTBLS21_arr.txt", "DNA microarray")]), _sheet(60, person=True), {"a_MTBLS21_arr.txt": _assay(60, tech="ms")})  # 3.4 EXCLUDED
    _study(off, "MTBLS22", _inv("MTBLS22", "Urine NMR cohort A", "Urine was collected from each participant.", cs, ["Group"], [("a_MTBLS22_nmr.txt", "NMR spectroscopy")]),
           _sheet(60, person=True), {"a_MTBLS22_nmr.txt": _assay(60)})                                                # duplicate of MTBLS10 (same inst, same title) -> one unit
    # exemplar-shaped fixture: the MTBLS1 column layout observed in run 3 (values synthetic), no person key, template description
    cols1 = ["Source Name", "Characteristics[Organism]", "Term Source REF", "Term Accession Number", "Characteristics[Organism part]", "Term Source REF", "Term Accession Number",
             "Characteristics[Variant]", "Term Source REF", "Term Accession Number", "Characteristics[Sample type]", "Term Source REF", "Term Accession Number", "Protocol REF", "Sample Name",
             "Factor Value[Gender]", "Term Source REF", "Term Accession Number", "Factor Value[Metabolic syndrome]", "Term Source REF", "Term Accession Number"]
    rows1 = [[f"SRC{i}", "Homo sapiens", "NCBITAXON", "http://purl.obolibrary.org/obo/NCBITaxon_9606", "urine", "BTO", "x", "-", "", "", "experimental sample", "OBI", "x", "Sample collection", f"ADG{i}",
              "Male" if i % 2 else "Female", "NCIT", "x", "diabetes mellitus" if i % 3 else "Control Group", "EFO", "x"] for i in range(132)]
    s1 = "\n".join(["\t".join(cols1)] + ["\t".join(r) for r in rows1]) + "\n"
    a1 = "\n".join(["Sample Name\tProtocol REF\tParameter Value[Extraction Method]\tExtract Name\tNMR Assay Name\tNormalization Name\tData Transformation Name\tMetabolite Assignment File"] +
                    [f"ADG{i}\tExtraction\tnone\tE{i}\tA{i}\tN{i}\tD{i}\tm_MTBLS1.tsv" for i in range(132)]) + "\n"
    _study(off, "MTBLS27", _inv("MTBLS27", "A metabolomic study of urinary changes in type 2 diabetes in human compared to the control group",
                                "Type 2 diabetes mellitus is the result of many conditions. Urine samples were analysed by NMR.", ["diabetes mellitus", "Human Study Subject", "untargeted metabolites"],
                                ["Gender", "Metabolic syndrome"], [("a_MTBLS27_nmr.txt", "NMR spectroscopy")]), s1, {"a_MTBLS27_nmr.txt": a1}, extra_names=("m_MTBLS1.tsv",))
    # v0.1.1 regression fixtures (owner's blocking finding): investigation lists the REST-style assay name, mirror holds the
    # public name — the actual MTBLS1 / MTBLS2 pairs observed in runs 2–3
    _study(off, "MTBLS28", _inv("MTBLS28", "Urine NMR cohort E", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS1_NMR_metabolite_profiling_NMR_spectroscopy.txt", "NMR spectroscopy")]),
           _sheet(60, person=True), {"a_MTBLS1_metabolite_profiling_NMR_spectroscopy.txt": _assay(60)})
    _study(off, "MTBLS29", _inv("MTBLS29", "Plasma MS cohort F", "One blood sample was collected from each participant.", cs, ["Group"], [("a_MTBLS2_LC-MS_metabolite_profiling_mass_spectrometry.txt", "mass spectrometry")]),
           _sheet(60, person=True, part="plasma"), {"a_MTBLS2_metabolite_profiling_mass_spectrometry.txt": _assay(60, tech="ms")})
    # two investigation assays, two differently named mirror files -> no unique pairing -> UNRESOLVED_STRUCTURE, not EXCLUDED
    _study(off, "MTBLS30", _inv("MTBLS30", "Two assays", "One blood sample was collected from each participant.", cs, ["Group"], [("a_X_pos.txt", "mass spectrometry"), ("a_X_neg.txt", "mass spectrometry")]),
           _sheet(60, person=True, part="plasma"), {"a_MTBLS30_pos.txt": _assay(60, tech="ms"), "a_MTBLS30_neg.txt": _assay(60, tech="ms")})
    # ambiguous QC label in a Factor Value -> UNRESOLVED; ambiguous label in Sample Name -> removed as QC
    sh = _sheet(60, person=True).replace("\tcase\n", "\treference group\n", 1)
    _study(off, "MTBLS31", _inv("MTBLS31", "Ambiguous factor", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS31_nmr.txt", "NMR spectroscopy")]), sh, {"a_MTBLS31_nmr.txt": _assay(60)})
    sh2 = _sheet(60, person=True, qc=2).replace("QC_sample_0", "pool_1").replace("QC_sample_1", "standard_2").replace("\tQC\tSample collection", "\tx\tSample collection")
    sh2 = "\n".join(ln.replace("\tQC\n", "\tx\n") for ln in sh2.splitlines(True))
    a2 = _assay(60, qc=2).replace("QC_sample_0", "pool_1").replace("QC_sample_1", "standard_2")
    _study(off, "MTBLS32", _inv("MTBLS32", "Pool names", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS32_nmr.txt", "NMR spectroscopy")]), sh2, {"a_MTBLS32_nmr.txt": a2})
    # two sample sheets -> UNRESOLVED_STRUCTURE
    _study(off, "MTBLS33", _inv("MTBLS33", "Two sheets", "x", cs, ["Group"], [("a_MTBLS33_nmr.txt", "NMR spectroscopy")]), _sheet(60, person=True), {"a_MTBLS33_nmr.txt": _assay(60)})
    off[url_for("FTPDIR", "MTBLS33")] = off[url_for("FTPDIR", "MTBLS33")].replace(b"</html>", b'<a href="s_MTBLS33_extra.txt">s2</a></html>'); off[url_for("FETCH", "MTBLS33", "s_MTBLS33_extra.txt")] = _sheet(5).encode()
    # D3: specific person ids overlapping with same affiliation -> CONFIRMED merge; generic ids elsewhere -> flag only
    def sp_sheet(n, prefix):
        return _sheet(n, person=True).replace("Characteristics[Subject ID]", "Characteristics[Participant]")
    s34 = _sheet(60, person=True); s35 = _sheet(60, person=True)
    s34 = "\n".join(ln.replace(f"\tP{i:04d}\t", f"\tHOSP-{i + 1000:05d}-K\t") if i < 60 else ln for i, ln in enumerate(s34.split("\n")[1:], 0)); s34 = _sheet(0, person=True).strip("\n") + "\n" + s34
    s35 = "\n".join(ln.replace(f"\tP{i:04d}\t", f"\tHOSP-{i + 1000:05d}-K\t") if i < 60 else ln for i, ln in enumerate(s35.split("\n")[1:], 0)); s35 = _sheet(0, person=True).strip("\n") + "\n" + s35
    _study(off, "MTBLS34", _inv("MTBLS34", "Hospital cohort serum", "One serum sample was collected from each participant.", cs, ["Group"], [("a_MTBLS34_ms.txt", "mass spectrometry")], affil="Inst B"), s34, {"a_MTBLS34_ms.txt": _assay(60, tech="ms")})
    _study(off, "MTBLS35", _inv("MTBLS35", "Hospital cohort urine", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS35_nmr.txt", "NMR spectroscopy")], affil="Inst B"), s35, {"a_MTBLS35_nmr.txt": _assay(60)})
    deg = _sheet(60, src_fn=lambda i: "STUDY")                                   # one Source Name for all rows (MTBLS28 pattern)
    _study(off, "MTBLS36", _inv("MTBLS36", "Urine cohort degenerate source", "Urine samples from lung cancer patients and controls.", cs, ["Group"], [("a_MTBLS36_ms.txt", "mass spectrometry")]), deg, {"a_MTBLS36_ms.txt": _assay(60, tech="ms")})
    emp = _sheet(60, person=True) + "SRC9999\tHomo sapiens\tNCBITAXON\tNCBITAXON:9606\turine\tP9999\tSample collection\t\tcontrol\n"   # row without Sample Name (MTBLS3 pattern)
    _study(off, "MTBLS37", _inv("MTBLS37", "Empty row", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS37_nmr.txt", "NMR spectroscopy")]), emp, {"a_MTBLS37_nmr.txt": _assay(60)})
    # A30 fixtures
    def org_rows(spec):   # spec: list of (organism, accession, sample_suffix)
        cols = ["Source Name", "Characteristics[Organism]", "Term Source REF", "Term Accession Number", "Characteristics[Organism part]", "Characteristics[Subject ID]", "Protocol REF", "Sample Name", "Factor Value[Group]"]
        rows = [[f"SRC{i:04d}", o, "NCBITAXON" if a else "", a, "urine", f"P{i:04d}", "Sample collection", f"S{i:04d}", "case"] for i, (o, a) in enumerate(spec)]
        return "\n".join(["\t".join(cols)] + ["\t".join(r) for r in rows]) + "\n"
    spec27 = [("Homo sapiens", "NCBITAXON:9606")] * 60 + [("blank", "")] * 3 + [("reference compound", "")] * 2
    _study(off, "MTBLS38", _inv("MTBLS38", "Human urine with blanks", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS38_nmr.txt", "NMR spectroscopy")]), org_rows(spec27), {"a_MTBLS38_nmr.txt": _assay(65)})
    spec_missing = [("Homo sapiens", "NCBITAXON:9606")] * 60 + [("", "")] * 2
    _study(off, "MTBLS39", _inv("MTBLS39", "Missing organism", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS39_nmr.txt", "NMR spectroscopy")]), org_rows(spec_missing), {"a_MTBLS39_nmr.txt": _assay(62)})
    spec_unsup = [("Homo sapiens", "NCBITAXON:9606")] * 60 + [("reference compound", "NCBITAXON:9606")] * 2
    _study(off, "MTBLS40", _inv("MTBLS40", "Unsupported reference", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS40_nmr.txt", "NMR spectroscopy")]), org_rows(spec_unsup), {"a_MTBLS40_nmr.txt": _assay(62)})
    spec_c = [("Homo sapiens", "NCBITAXON:9606")] * 60 + [("reference compound", "NCBITAXON:9606")] * 2
    sh_c = org_rows(spec_c).replace("\tS0060\tcase", "\tS0060\tQC").replace("\tS0061\tcase", "\tS0061\tQC")
    _study(off, "MTBLS41", _inv("MTBLS41", "Conflict resolved", "One urine sample was collected from each participant.", cs, ["Group"], [("a_MTBLS41_nmr.txt", "NMR spectroscopy")]), sh_c, {"a_MTBLS41_nmr.txt": _assay(62)})
    # v0.1.7 fixtures: percent-encoded hrefs (MTBLS242 pattern) and an index-listed file that 404s
    sp42 = "s_MTBLS42 obesity.txt"; ap42 = "a_MTBLS42_metabolite profiling_mass spectrometry.txt"
    _study(off, "MTBLS42", _inv("MTBLS42", "Encoded names", "One blood sample was collected from each participant.", cs, ["Group"], [(ap42, "mass spectrometry")]).replace("s_MTBLS42.txt", sp42),
           _sheet(60, person=True, part="plasma"), {ap42: _assay(60, tech="ms")})
    off[url_for("FETCH", "MTBLS42", sp42)] = off.pop(url_for("FETCH", "MTBLS42", "s_MTBLS42.txt"))
    off[url_for("FTPDIR", "MTBLS42")] = ('<html><a href="i_Investigation.txt">i</a><a href="' + quote(sp42) + '">s</a><a href="' + quote(ap42) + '">a</a></html>').encode()
    _study(off, "MTBLS43", _inv("MTBLS43", "Listed but missing", "x", cs, ["Group"], [("a_MTBLS43_ms.txt", "mass spectrometry")]), _sheet(60, person=True), {"a_MTBLS43_ms.txt": _assay(60, tech="ms")})
    off.pop(url_for("FETCH", "MTBLS43", "a_MTBLS43_ms.txt"))
    off[url_for("FTPDIR", "MTBLS23")] = HTTPError(url_for("FTPDIR", "MTBLS23"), 404, "x", None, None)                 # INDEX_404
    off[url_for("FTPDIR", "MTBLS24")] = http.client.RemoteDisconnected("x")                                           # TRANSPORT_FAILED
    _study(off, "MTBLS25", _inv("MTBLS25", "Bad sheet", "x", cs, ["Group"], [("a_MTBLS25_ms.txt", "mass spectrometry")]), "Source Name\tSample Name\tIntensity\nA\tS\t3.2\n", {"a_MTBLS25_ms.txt": _assay(1, tech="ms")})  # FORMAT_REFUSED
    off[url_for("FTPDIR", "MTBLS26")] = b'<html><a href="i_Investigation.txt">i</a><a href="s_MTBLS26.txt">s</a></html>'
    off[url_for("FETCH", "MTBLS26", "i_Investigation.txt")] = Redirect(302, FTP + "/MTBLS26/m_x.tsv")                  # redirect -> FORMAT_REFUSED (never followed)
    return off


def cmd_selftest(a):
    ok = True

    def check(c, msg):
        nonlocal ok
        print(("PASS " if c else "FAIL ") + msg); ok &= bool(c)
    import tempfile
    off = _mock(); accs = [f"MTBLS{i}" for i in range(10, 44)]
    d = Path(tempfile.mkdtemp()); t = Transport(offline=off, pause=0)
    summ = run_census(t, accs, d)
    recs = {json.loads(l)["accession"]: json.loads(l) for l in (d / "stageA_census_checkpoint.jsonl").read_text().splitlines()}
    g = lambda a, k: recs[a].get(k)
    check(g("MTBLS10", "status") == "ELIGIBLE" and g("MTBLS10", "subject_grade") == "R1" and g("MTBLS10", "n_sources") == 60 and g("MTBLS10", "n_qc_rows") == 3 and g("MTBLS10", "L1_strict") == 1, "R1 cohort n=60 with 3 QC rows removed -> ELIGIBLE strict")
    check(g("MTBLS11", "status") == "ELIGIBLE" and g("MTBLS11", "subject_grade") == "R2_T12" and g("MTBLS11", "L1_strict") == 1 and g("MTBLS11", "platform_class") == "MS", "R2_T12 cohort -> ELIGIBLE strict, MS")
    check(g("MTBLS12", "status") == "ELIGIBLE" and g("MTBLS12", "subject_grade") == "R2_COUNT_SUPPORTED" and g("MTBLS12", "L1_strict") == 0 and g("MTBLS12", "L1_sensitivity") == 1, "T3 count cohort -> sensitivity only")
    check(g("MTBLS13", "status") == "UNRESOLVED" and g("MTBLS13", "subject_grade") == "SOURCE_UNIQUE", "SOURCE_UNIQUE -> UNRESOLVED, never confirmed")
    check(g("MTBLS14", "status") == "EXCLUDED" and g("MTBLS14", "rule") == "M1", "non-human -> EXCLUDED at M1")
    check(g("MTBLS15", "status") == "EXCLUDED" and "Source Name carries" in g("MTBLS15", "reason"), "source with 2 samples -> EXCLUDED")
    check(g("MTBLS16", "status") == "EXCLUDED" and g("MTBLS16", "rule") == "M4", "two matrices -> EXCLUDED at M4")
    check(g("MTBLS17", "status") == "EXCLUDED" and "49" in g("MTBLS17", "reason"), "n=49 -> EXCLUDED (A5)")
    check(g("MTBLS18", "status") == "UNRESOLVED" and "Layer 1" in g("MTBLS18", "subject_evidence"), "longitudinal descriptor blocks R2 -> UNRESOLVED")
    check(g("MTBLS19", "status") == "ELIGIBLE" and "with-note" in g("MTBLS19", "temporal"), "strong temporal factor, one sample per source -> CONFIRMED-with-note")
    check(g("MTBLS20", "status") == "EXCLUDED" and g("MTBLS20", "rule") == "M5" and "Pooling" in g("MTBLS20", "reason"), "pooling Protocol REF on the path -> EXCLUDED")
    check(g("MTBLS21", "status") == "EXCLUDED" and g("MTBLS21", "rule") == "3.4", "no MS/NMR assay -> EXCLUDED at 3.4")
    check(g("MTBLS23", "status") == "INDEX_404" and g("MTBLS24", "status") == "TRANSPORT_FAILED" and g("MTBLS25", "status") == "FORMAT_REFUSED" and g("MTBLS26", "status") == "ACCESS_REFUSED",
          "index 404 / transport failure / bad sheet / redirect -> per-study statuses (A20 classes); census continued")
    check(g("MTBLS27", "status") == "UNRESOLVED" and g("MTBLS27", "subject_grade") == "SOURCE_UNIQUE" and g("MTBLS27", "n_sources") == 132 and g("MTBLS27", "matrix") == "urine"
          and "Layer 1 fail" in g("MTBLS27", "subject_evidence"), "exemplar-shaped fixture (MTBLS1 layout, no person key, no design descriptor) -> UNRESOLVED SOURCE_UNIQUE, n=132, matrix urine")
    check(g("MTBLS28", "status") == "ELIGIBLE" and g("MTBLS28", "assay_file") == "a_MTBLS1_metabolite_profiling_NMR_spectroscopy.txt" and g("MTBLS28", "assay_technology_source") == "unique-unmatched pairing"
          and g("MTBLS28", "assay_investigation_name") == "a_MTBLS1_NMR_metabolite_profiling_NMR_spectroscopy.txt" and g("MTBLS28", "platform_class") == "NMR",
          "reconciliation: MTBLS1 pair (investigation NMR_ name vs mirror name) -> paired, NMR, both names preserved")
    check(g("MTBLS29", "status") == "ELIGIBLE" and g("MTBLS29", "assay_technology_source") == "unique-unmatched pairing" and g("MTBLS29", "platform_class") == "MS", "reconciliation: MTBLS2 pair (LC-MS_ name) -> paired, MS")
    check(g("MTBLS30", "status") == "UNRESOLVED_STRUCTURE" and g("MTBLS30", "rule") == "3.4" and "unmatched" in g("MTBLS30", "assay_reconciliation"), "two unmatched assays -> no pairing invented -> UNRESOLVED_STRUCTURE, not EXCLUDED")
    check(g("MTBLS31", "status") == "UNRESOLVED" and g("MTBLS31", "rule") == "M2" and g("MTBLS31", "n_ambiguous_qc_rows") == 1, "ambiguous QC label ('reference group') in a Factor Value -> UNRESOLVED, not removed")
    check(g("MTBLS32", "status") == "ELIGIBLE" and g("MTBLS32", "n_qc_rows") == 2 and g("MTBLS32", "n_sources") == 60, "ambiguous tokens in Sample Name (pool_1, standard_2) -> analytical context -> removed as QC")
    check(g("MTBLS33", "status") == "UNRESOLVED_STRUCTURE" and g("MTBLS33", "rule") == "A23", "two sample sheets -> UNRESOLVED_STRUCTURE (A23), not a biological exclusion")
    us2 = list(csv.DictReader((d / "stageA_M_units.csv").open())); fl = list(csv.DictReader((d / "stageA_M_overlap_flags.csv").open()))
    check(any(u["accessions"] == "MTBLS34|MTBLS35" and "specific person-identifier overlap" in u["merge_reason"] for u in us2), "D3: specific person-id hash overlap + same affiliation -> CONFIRMED merge")
    check(any(f["accession_a"] == "MTBLS10" and f["identifier_kind"] in ("person", "source", "sample") and "SUSPECTED" in f["status"] for f in fl) and not any(u["accessions"] == "MTBLS10|MTBLS11" for u in us2),
          "D3: generic identifier overlap (P0000…, SRC0000…) across different cohorts -> review flag only, no merge")
    check(all("_ids" not in r for r in csv.DictReader((d / "stageA_M_cohorts.csv").open())) and not any("SRC0000" in (d / "stageA_census_checkpoint.jsonl").read_text() for _ in [0]), "no raw identifier in cohorts.csv or checkpoint (hashes only)")
    hs = run_census(Transport(offline=off, pause=0), ["MTBLS10/../MTBLS11"], Path(tempfile.mkdtemp()))
    check("halted_protocol_violation" in hs and "REFUSED" in hs["halted_protocol_violation"], "A20: guard refusal -> global halt recorded, census stops")
    codes = {}
    for label, sm in (("halted", hs), ("completed", summ)):
        try:
            import io, contextlib
            with contextlib.redirect_stdout(io.StringIO()):
                finish_run(sm)
        except SystemExit as e:
            codes[label] = e.code
    check(codes == {"halted": 1, "completed": 0} and run_exit_code(hs) == 1, "C1: the run's terminal exit path returns 1 on a protocol-violation halt and 0 on completion")
    check(any(f["identifier_kind"] == "title" and {f["accession_a"], f["accession_b"]} == {"MTBLS10", "MTBLS22"} for f in fl), "C2: title-similarity pair is a SUSPECTED_OVERLAP flag")
    fl_txt = (d / "stageA_census_fetch_log.txt").read_text(); n_before = sum(1 for ln in fl_txt.splitlines() if ln and not ln.startswith("#"))
    check(n_before > 0 and summ["n_requests_cumulative"] == n_before and fl_txt.startswith("# run ") and all("elapsed_s=" in ln for ln in fl_txt.splitlines() if ln and not ln.startswith("#")),
          "v0.1.4: fetch log journaled with a run header, cumulative count, elapsed seconds on every line")
    check(any("\tERR\t" in ln and "elapsed_s=" in ln for ln in fl_txt.splitlines()) and any("\t404\t" in ln for ln in fl_txt.splitlines()), "C6: failed requests carry elapsed time too")
    check("'" not in g("MTBLS37", "reason") and g("MTBLS37", "n_rows_unidentified") == 1, "C5: M2 diagnostic is a count, no identifier in the reason")
    check("D3_metadata_units" in summ and summ["L3_confirmed_status"] == "NOT_ASSESSED" and "L3_units" not in summ, "C7: D3 metadata units named as such; L3_confirmed_status NOT_ASSESSED")
    class Boom(Exception):
        pass
    offb = dict(off)
    def boom():
        raise Boom("simulated process termination")
    offb[url_for("FTPDIR", "MTBLS11")] = boom
    db = Path(tempfile.mkdtemp()); tb = Transport(offline=offb, pause=0)
    try:
        run_census(tb, ["MTBLS10", "MTBLS11"], db); check(False, "interruption fixture raised")
    except Boom:
        jl = (db / "stageA_census_fetch_log.txt").read_text().splitlines()
        check(sum(1 for ln in jl if ln and not ln.startswith("#")) == 4 and (db / "stageA_census_checkpoint.jsonl").read_text().count("\n") == 1,
              "v0.1.4: after an abrupt termination the journal holds every completed request (4) and the checkpoint the completed study (1)")
    check(summ["artefacts_not_for_publication"] == ["stageA_census_checkpoint.jsonl"] and "pseudonymisation" in summ["identifier_hashing"], "C3: pseudonymisation stated; checkpoint marked not for publication")
    check(g("MTBLS36", "status") == "UNRESOLVED" and g("MTBLS36", "subject_grade") == "SOURCE_DEGENERATE" and g("MTBLS36", "n_sources") == 1, "A27: one Source Name for all samples -> UNRESOLVED (source node uninformative), not EXCLUDED")
    check(g("MTBLS15", "status") == "EXCLUDED", "A27: > 1 source with multiplicity remains EXCLUDED")
    check(g("MTBLS37", "status") == "UNRESOLVED" and g("MTBLS37", "rule") == "M2" and g("MTBLS37", "n_rows_without_sample_name") == 1, "v0.1.3: a row without Sample Name -> UNRESOLVED (M2) and counted separately")
    check(g("MTBLS38", "status") == "ELIGIBLE" and g("MTBLS38", "n_qc_rows") == 5 and g("MTBLS38", "n_rows_organism_qc_removed") == 5 and g("MTBLS38", "n_sources") == 60
          and "blank" in g("MTBLS38", "organism_values_raw") and g("MTBLS38", "organism_values") == "homo sapiens", "A30: explicit blank / reference-compound organism labels (no conflicting taxonomy annotation) -> QC; raw values preserved; human rows evaluated")
    check(g("MTBLS39", "status") == "UNRESOLVED" and g("MTBLS39", "rule") == "M1" and g("MTBLS39", "n_rows_organism_missing") == 2, "A30: missing organism without QC evidence -> UNRESOLVED, not QC, not EXCLUDED")
    check(g("MTBLS40", "status") == "UNRESOLVED" and g("MTBLS40", "n_rows_organism_analytical_conflict") == 2, "A30: 'reference compound' with a conflicting taxonomy accession and no other QC evidence -> UNRESOLVED")
    check(g("MTBLS41", "status") == "ELIGIBLE" and g("MTBLS41", "n_qc_rows") == 2 and g("MTBLS41", "n_rows_organism_analytical_conflict") == 0,
          "A30: conflicting taxonomy accession resolved by other QC evidence (Factor Value 'QC') -> QC rows removed")
    check(g("MTBLS14", "status") == "EXCLUDED" and g("MTBLS14", "rule") == "M1", "A30: explicit non-human organism still EXCLUDED under M1")
    check(g("MTBLS42", "status") == "ELIGIBLE" and any("MTBLS42/a_MTBLS42_metabolite%20profiling_mass%20spectrometry.txt" in l[0] for l in t.log)
          and not any("%25" in l[0] for l in t.log), "v0.1.7 D-full-2: percent-encoded index hrefs decoded once, re-encoded once (no %2520)")
    check(g("MTBLS43", "status") == "FETCH_HTTP_404" and g("MTBLS26", "status") == "ACCESS_REFUSED", "v0.1.7 D-full-3: index-listed file returning 404 -> FETCH_HTTP_404; redirect stays ACCESS_REFUSED")
    check(g("MTBLS10", "census_version") == VERSION and g("MTBLS10", "person_key_columns") == "characteristics[subject id]" and g("MTBLS36", "person_key_columns") == "",
          "v0.1.7: census_version and person_key_columns recorded on every screened record")
    check(retry_eligible({"status": "ACCESS_REFUSED", "reason": "https://x/MTBLS242/s_a%2520b.txt: HTTP 404"}) and not retry_eligible({"status": "ACCESS_REFUSED", "reason": "x: HTTP 302 redirect not followed"})
          and retry_eligible({"status": "FETCH_HTTP_404", "reason": ""}) and not retry_eligible({"status": "EXCLUDED", "reason": ""}), "v0.1.7: retry scope = TRANSPORT_FAILED, FETCH_HTTP_*, legacy double-encoding ACCESS_REFUSED only")
    offo = dict(off)
    for k in range(50, 57):
        offo[url_for("FTPDIR", f"MTBLS{k}")] = URLError("[Errno 8] nodename nor servname provided, or not known")
    do = Path(tempfile.mkdtemp()); so = run_census(Transport(offline=offo, pause=0), [f"MTBLS{k}" for k in range(50, 57)], do)
    nck = sum(1 for ln in (do / "stageA_census_checkpoint.jsonl").read_text().splitlines() if ln.strip())
    check("halted_outage" in so and run_exit_code(so) == 2 and nck == OUTAGE_K, f"v0.1.7 D-full-1: {OUTAGE_K} consecutive transport failures stop the run (exit 2); the rest of the list is not touched")
    check(not any("m_x.tsv" in l[0] for l in t.log) and not any("/MTBLS26/m_x.tsv" in l[0] for l in t.log) and all(l[0].startswith(FTP) for l in t.log), "no MAF, no redirect destination, no non-mirror host requested")
    check(summ["L1_strict"] == 12 and summ["L1_sensitivity"] == 13 and summ["L1_unresolved"] == 8 and summ["L2_confirmed_status"] == "NOT_ASSESSED"
          and summ["L2_preliminary"]["unresolved_FEATURE_COUNT_UNKNOWN"] == 13, "summary: strict 12, sensitivity 13, unresolved 8; L2 NOT_ASSESSED; FEATURE_COUNT_UNKNOWN")
    us = list(csv.DictReader((d / "stageA_M_units.csv").open()))
    check(not any({"MTBLS10", "MTBLS22"} <= set(u["accessions"].split("|")) for u in us), "C2: same affiliation + similar title is NOT merged")
    n_log = len(t.log); t2 = Transport(offline=off, pause=0); run_census(t2, accs, d)
    check(len(t2.log) == 0, "resume: completed accessions are not re-requested")
    fl2 = (d / "stageA_census_fetch_log.txt").read_text(); n_after = sum(1 for ln in fl2.splitlines() if ln and not ln.startswith("#"))
    s2 = json.loads((d / "stageA_M_summary.json").read_text())
    check(n_after == n_before and fl2.count("# run ") == 2 and s2["n_requests_cumulative"] == n_before and s2["n_requests_this_invocation"] == 0,
          "v0.1.3: resume appends an empty run block; the original fetch log survives; cumulative count unchanged")
    t3 = Transport(offline=off, pause=0); run_census(t3, accs, d, retry_pass=True)
    rq = {l[0].split("/public/")[1].split("/")[0] for l in t3.log}
    check(rq == {"MTBLS24", "MTBLS43"}, "retry pass re-requests only TRANSPORT_FAILED (MTBLS24) and FETCH_HTTP_* (MTBLS43) studies")
    man = d / MANIFEST; man.write_text("accession,in_rest,in_mirror,status\nMTBLS1,1,1,UNIVERSE\nMTBLS2,1,1,UNIVERSE\nMTBLS10,1,1,UNIVERSE\nMTBLS9,0,1,UNRESOLVED\n")
    check(load_universe(man, sha_file(man)) == ["MTBLS10"], "universe: UNIVERSE rows only, exemplars excluded, UNRESOLVED never screened")
    try:
        load_universe(man, "0" * 64); check(False, "manifest hash mismatch refused")
    except Refused:
        check(True, "manifest hash mismatch refused")
    man.write_text("accession,in_rest,in_mirror,status\nMTBLS10,1,1,PROVISIONAL_UNIVERSE\n")
    try:
        load_universe(man, sha_file(man)); check(False, "provisional manifest refused")
    except Refused:
        check(True, "provisional manifest refused")
    for kind, acc, name in (("FETCH", "MTBLS10", "m_x.tsv"), ("FETCH", "MTBLS10", "FILES/x.raw"), ("FETCH", "MTBLS10/../MTBLS11", "i_Investigation.txt"), ("DOWNLOAD", "MTBLS10", "")):
        try:
            url_for(kind, acc, name); check(False, f"guard refuses {kind} {acc} {name}")
        except Refused:
            check(True, f"guard refuses {kind} {acc} {name}")
    ap = d / AUTH_FILE; ap.write_text(f"LIVE STAGE A CENSUS AUTHORISED 2026-10-25 universe={'a' * 64} protocol={PROTOCOL} owner=TEST\n")
    check(auth_status(ap)["valid"] and auth_status(ap)["universe_sha256"] == "a" * 64 and auth_status(ap)["pilot"] is None, "authorisation line carries the universe hash; no pilot")
    ap.write_text(f"LIVE STAGE A CENSUS AUTHORISED 2026-10-25 universe={'a' * 64} protocol={PROTOCOL} owner=TEST pilot=25\n")
    check(auth_status(ap)["valid"] and auth_status(ap)["pilot"] == 25, "authorisation line with pilot=25 recognised (A21)")
    print("SELFTEST PASSED" if ok else "SELFTEST FAILED"); sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest")
    for name in ("plan", "run"):
        p = sub.add_parser(name)
        p.add_argument("--protocol", default=PROTOCOL); p.add_argument("--auth", default=AUTH_FILE); p.add_argument("--manifest", default=MANIFEST)
        p.add_argument("--out", default="stageA_census_out"); p.add_argument("--pilot", type=int, default=0); p.add_argument("--retry-pass", action="store_true")
    a = ap.parse_args()
    {"selftest": cmd_selftest, "plan": cmd_plan, "run": cmd_run}[a.cmd](a)


if __name__ == "__main__":
    main()
