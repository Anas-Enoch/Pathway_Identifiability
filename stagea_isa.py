#!/usr/bin/env python3
"""
stagea_isa.py  v0.2.2   Stage A — ISA-Tab structural validators and the frozen A10 recognition function.

Offline module. No network. Implements protocol v0.4 §2 D1 (R1, R2 two-layer), §3.3 validators, §3.5 M6 temporal
detector. Everything here is deterministic text processing on files already on disk; it never opens a MAF.

Evidence grades returned by `subject_evidence()`:
  R1                 explicit person identifier on the source node, one sample per person
  R2_T12             Layer 1 passes and a T1/T2 sentence (one sample from each participant) qualifies  -> strict L1
  R2_COUNT_SUPPORTED Layer 1 passes and only a T3/T4/T4b count sentence qualifies, integer == source count -> sensitivity L1
  SOURCE_UNIQUE      unique source labels, one sample each, no affirmative evidence                      -> UNRESOLVED
  UNRESOLVED         anything else (missing fields, sentinels, contradictions, Layer 1 fail)
  EXCLUDED           a source/person with > 1 sample, pooling/splitting protocol, or demonstrated repeated sampling

Usage: python3 -I stagea_isa.py selftest
"""
from __future__ import annotations
import re, sys

VERSION = "0.2.2"
MISSING = {"", "-", "na", "n/a", "none", "null", "unknown", "nan", "not applicable", "not available"}

# ---------------------------------------------------------------- normalisation
def norm(s: str) -> str:
    s = str(s or "").casefold().replace("_", " ").replace("-", " ")
    return re.sub(r"\s+", " ", s).strip()


def is_missing(v) -> bool:
    return norm(v) in MISSING


# ---------------------------------------------------------------- A8 Layer 1: design descriptors (frozen sets, owner 2026-10-25)
CROSS_SECTIONAL = ("cross sectional", "case control", "case control study", "matched case control",
                   "population based cross sectional", "observational cross sectional")
LONGITUDINAL = ("longitudinal", "time series", "repeated measures", "crossover", "cross over", "prospective cohort",
                "follow up", "before and after", "pre post intervention", "paired sampling", "serial sampling")


def _phrase(p: str) -> re.Pattern:
    return re.compile(r"(?<![a-z0-9])" + re.escape(p) + r"(?![a-z0-9])")


_POS = [_phrase(p) for p in CROSS_SECTIONAL]
_NEG = [_phrase(p) for p in LONGITUDINAL]


def layer1(design_types: list[str], description: str) -> dict:
    text = norm(" | ".join(list(design_types) + [description]))
    pos = [p for p, rx in zip(CROSS_SECTIONAL, _POS) if rx.search(text)]
    neg = [p for p, rx in zip(LONGITUDINAL, _NEG) if rx.search(text)]
    return {"pass": bool(pos) and not neg, "positive": pos, "negative": neg}


# ---------------------------------------------------------------- A8 Layer 2: affirmative one-person-one-source statement (A10)
P = r"(participants?|subjects?|patients?|individuals?|donors?|volunteers?|persons?|people|women|men|children|cases|controls)"
ONE = r"(one|a single|single|1)"
SAMPLE = r"(samples?|specimens?|biopsy|biopsies|blood draws?|collections?)"
ADJ = r"( [a-z]+){0,2}"          # up to two qualifying words, e.g. "one BLOOD sample", "a single PLASMA specimen"
T1 = re.compile(ONE + ADJ + r" " + SAMPLE + r" (was |were )?(collected|obtained|taken|drawn|provided|analy[sz]ed) (from|for|by) (each|every) " + P)
T2 = re.compile(r"(each|every) " + P + r" (contributed|provided|donated|gave|supplied) " + ONE + ADJ + r" " + SAMPLE)
T3 = re.compile(SAMPLE + r" (were |was )?(collected|obtained|taken) from (\d+) (unique|distinct|different|independent|separate) " + P)
T4 = re.compile(r"(\d+) (unique|distinct|different|independent|separate) " + P + r" (were |was )?(enrolled|recruited|included|studied|sampled|analy[sz]ed)")
T4b = re.compile(r"(enrolled|recruited|included|studied|sampled|analy[sz]ed) (\d+) (unique|distinct|different|independent|separate) " + P)
TEMPLATES = {"T1": T1, "T2": T2, "T3": T3, "T4": T4, "T4b": T4b}

NEGATION = [_phrase(t) for t in ("not", "no", "never", "n't", "without", "except", "rather than", "instead of", "unless")]
PLURAL = [_phrase(t) for t in ("samples per", "multiple samples", "several samples", "serial", "repeated", "longitudinal",
                               "each visit", "each time point", "each timepoint", "subset of", "some of the", "twice", "three times",
                               "four times", "times a", "per week", "per day", "per month", "weekly", "daily", "monthly", "every week",
                               "every day", "every month")]
# sentences anywhere in the description that demonstrate repeated sampling block R2 regardless of a positive sentence
REPEATED_SAMPLING = [_phrase(t) for t in ("sampled at", "samples were collected at", "two time points", "three time points",
                                          "multiple time points", "before and after", "pre and post", "paired samples",
                                          "serial samples", "repeated sampling", "at each visit", "at baseline and",
                                          "follow up samples", "longitudinal samples")]


def sentences(description: str) -> list[str]:
    return [s.strip() for s in re.split(r"[.;:]\s+|[.;:]$", norm(description)) if s.strip()]


def layer2(description: str, n_sources: int) -> dict:
    """Returns the best grade among qualifying sentences; records every match and every guard trigger."""
    out = {"grade": None, "matches": [], "blocked": [], "repeated_sampling": []}
    for s in sentences(description):
        for rx in REPEATED_SAMPLING:
            if rx.search(s):
                out["repeated_sampling"].append(s)
    for s in sentences(description):
        for name, rx in TEMPLATES.items():
            m = rx.search(s)
            if not m:
                continue
            if any(g.search(s) for g in NEGATION):
                out["blocked"].append((name, "negation", s)); continue
            if any(g.search(s) for g in PLURAL):
                out["blocked"].append((name, "plural", s)); continue
            if name in ("T3", "T4", "T4b"):
                n_doc = int([g for g in m.groups() if g and g.isdigit()][0])
                if n_doc != n_sources:
                    out["blocked"].append((name, f"integer {n_doc} != sources {n_sources}", s)); continue
            out["matches"].append((name, s))
    if out["repeated_sampling"]:
        out["grade"] = None
        return out
    names = {n for n, _ in out["matches"]}
    if names & {"T1", "T2"}:
        out["grade"] = "R2_T12"
    elif names & {"T3", "T4", "T4b"}:
        out["grade"] = "R2_COUNT_SUPPORTED"
    return out


# ---------------------------------------------------------------- R1 person identifier on the source node
PERSON_KEY = re.compile(r"^characteristics\[(subject|participant|individual|donor|patient|volunteer|person)( ?(id|identifier|number|no|code))?\]$")


def subject_evidence(sheet_rows: list[dict], design_types: list[str], description: str, min_n: int = 50) -> dict:
    """sheet_rows: biological rows (QC removed) of the chosen assay's samples, as dicts with normalised column names."""
    rec = {"grade": "UNRESOLVED", "reason": "", "n_sources": 0, "n_samples": len(sheet_rows)}
    if not sheet_rows:
        rec["reason"] = "no biological rows"; return rec
    for r in sheet_rows:
        if is_missing(r.get("sample name")) or is_missing(r.get("source name")):
            rec["reason"] = "missing or sentinel Sample Name / Source Name"; return rec
    samples = [norm(r["sample name"]) for r in sheet_rows]
    if len(set(samples)) != len(samples):
        rec["reason"] = "duplicate Sample Name"; return rec
    sources = [norm(r["source name"]) for r in sheet_rows]
    per_source = {}
    for s in sources:
        per_source[s] = per_source.get(s, 0) + 1
    rec["n_sources"] = len(per_source)
    if max(per_source.values()) > 1:
        rec.update(grade="EXCLUDED", reason="a Source Name carries > 1 sample"); return rec
    # R1
    pkeys = sorted(k for k in sheet_rows[0] if PERSON_KEY.match(norm(k)))
    if pkeys:
        for k in pkeys:
            vals = [norm(r.get(k)) for r in sheet_rows]
            if any(v in MISSING for v in vals):
                rec["reason"] = f"person key '{k}' has missing values"; return rec
            if len(set(vals)) != len(vals):
                rec.update(grade="EXCLUDED", reason=f"person key '{k}': > 1 sample per person"); return rec
        if rec["n_sources"] < min_n:
            rec.update(grade="EXCLUDED", reason=f"n={rec['n_sources']} < {min_n}"); return rec
        rec.update(grade="R1", reason=f"person key(s) {pkeys}, one sample each, n={rec['n_sources']}"); return rec
    # R2
    l1 = layer1(design_types, description)
    rec["layer1"] = l1
    if not l1["pass"]:
        rec.update(grade="SOURCE_UNIQUE", reason=f"Layer 1 fail: positive={l1['positive']} negative={l1['negative']}"); return rec
    l2 = layer2(description, rec["n_sources"])
    rec["layer2"] = {k: v for k, v in l2.items()}
    if l2["repeated_sampling"]:
        rec.update(grade="EXCLUDED", reason="description demonstrates repeated sampling"); return rec
    if l2["grade"] is None:
        rec.update(grade="SOURCE_UNIQUE", reason="Layer 2: no qualifying sentence"); return rec
    if rec["n_sources"] < min_n:
        rec.update(grade="EXCLUDED", reason=f"n={rec['n_sources']} < {min_n}"); return rec
    rec.update(grade=l2["grade"], reason=f"Layer 1 {l1['positive']}; Layer 2 {[n for n, _ in l2['matches']]}")
    return rec


# ---------------------------------------------------------------- M6 temporal detector (two levels; underscore/camel-case aware)
def _flex(term: str) -> re.Pattern:
    words = term.split()
    return re.compile(r"(?<![a-z0-9])" + r"[ _]?".join(map(re.escape, words)) + r"(?![a-z0-9])")


STRONG_TIME = [_flex(t) for t in ("timepoint", "time point", "visit", "follow up", "collection time", "sampling day",
                                   "sampling week", "sampling month", "baseline", "endpoint", "day", "week", "month", "hour")]
AMBIG_TIME = [_flex(t) for t in ("stage", "phase", "pre", "post", "collection", "cycle", "duration", "time")]


def split_camel(k: str) -> str:
    return re.sub(r"([a-z])([A-Z])", r"\1 \2", str(k))


def temporal_class(factor_name: str) -> str:
    n = norm(split_camel(factor_name))
    if any(rx.search(n) for rx in STRONG_TIME):
        return "strong"
    if any(rx.search(n) for rx in AMBIG_TIME):
        return "ambiguous"
    return "none"


# ---------------------------------------------------------------- ISA structural validators (per file class)
I_SECTIONS = {"ONTOLOGY SOURCE REFERENCE", "INVESTIGATION", "INVESTIGATION PUBLICATIONS", "INVESTIGATION CONTACTS", "STUDY",
              "STUDY DESIGN DESCRIPTORS", "STUDY PUBLICATIONS", "STUDY FACTORS", "STUDY ASSAYS", "STUDY PROTOCOLS", "STUDY CONTACTS"}
_TERM = r"( Term Accession Number| Term Source REF)?"
_PERSON = r"(Last Name|First Name|Mid Initials|Email|Phone|Fax|Address|Affiliation|Roles" + _TERM + r")"
_PUB = r"(PubMed ID|Publication DOI|Publication Author List|Publication Title|Publication Status" + _TERM + r")"
I_FIELDS = {  # per section: regex over the full label (ISA-Tab 1.0 investigation file); Comment[...] allowed everywhere
    "ONTOLOGY SOURCE REFERENCE": r"Term Source (Name|File|Version|Description)",
    "INVESTIGATION": r"Investigation (Identifier|Title|Description|Submission Date|Public Release Date)",
    "INVESTIGATION PUBLICATIONS": r"Investigation " + _PUB,
    "INVESTIGATION CONTACTS": r"Investigation Person " + _PERSON,
    "STUDY": r"Study (Identifier|Title|Description|Submission Date|Public Release Date|File Name)",
    "STUDY DESIGN DESCRIPTORS": r"Study Design Type" + _TERM,
    "STUDY PUBLICATIONS": r"Study " + _PUB,
    "STUDY FACTORS": r"Study Factor (Name|Type" + _TERM + r")",
    "STUDY ASSAYS": r"Study Assay (File Name|Measurement Type" + _TERM + r"|Technology Type" + _TERM + r"|Technology Platform)",
    "STUDY PROTOCOLS": r"Study Protocol (Name|Type" + _TERM + r"|Description|URI|Version|Parameters Name" + _TERM + r"|Components Name|Components Type" + _TERM + r")",
    "STUDY CONTACTS": r"Study Person " + _PERSON,
}
I_FIELD_RX = {sec: re.compile(r"^(" + rx + r"|Comment\[.+\])$") for sec, rx in I_FIELDS.items()}
# v0.2.1 (A-4): the sample-table allowlist is likewise completed to the ISA-Tab 1.0 vocabulary (Material Type,
# Description, Provider, Performer, Date) from the specification, before the s_ file of any exemplar has been seen.
S_HEADERS = re.compile(r"^(Source Name|Sample Name|Characteristics\[.+\]|Factor Value\[.+\]|Protocol REF|Parameter Value\[.+\]|"
                       r"Material Type|Description|Provider|Performer|Date|Term Source REF|Term Accession Number|Unit|Comment\[.+\])$")
# v0.2.1 (A-4, after run 3 attempt 2 FORMAT-REFUSED the MTBLS1 assay sheet on 'Normalization Name' and
# 'Data Transformation Name'): the assay allowlist is completed to the ISA-Tab 1.0 assay-table node, file and
# attribute column names (process nodes Normalization Name / Data Transformation Name / Scan Name / Hybridization Assay
# Name; Array Design REF; Array Data File, Array Data Matrix File, Derived Array Data File, Derived Array Data Matrix
# File, Image File, Spot Picking File, Protein Assignment File, Peptide Assignment File, Post Translational Modification
# Assignment File; material attributes Characteristics[...] / Factor Value[...] / Material Type). Source: the ISA-Tab
# 1.0 specification's column vocabulary, not the observed file. Values are still never read beyond header and width
# checks; the numeric guard is unchanged (a numeric Factor Value column refuses the file, which is disclosed, not
# bypassed).
A_HEADERS = re.compile(r"^(Sample Name|Protocol REF|Parameter Value\[.+\]|Extract Name|Labeled Extract Name|Label|Assay Name|MS Assay Name|NMR Assay Name|"
                       r"Hybridization Assay Name|Scan Name|Normalization Name|Data Transformation Name|Array Design REF|"
                       r"Raw Spectral Data File|Derived Spectral Data File|Raw Data File|Derived Data File|Metabolite Assignment File|"
                       r"Array Data File|Array Data Matrix File|Derived Array Data File|Derived Array Data Matrix File|Image File|Spot Picking File|"
                       r"Protein Assignment File|Peptide Assignment File|Post Translational Modification Assignment File|"
                       r"Free Induction Decay Data File|Acquisition Parameter Data File|Characteristics\[.+\]|Factor Value\[.+\]|Material Type|"
                       r"Term Source REF|Term Accession Number|Unit|Comment\[.+\]|Date|Performer)$")
NUMERIC = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
# v0.2.2 (A26, after pilot 25: 17/25 real sample sheets refused at 6–22 % numeric cells — ages, doses, replicate numbers,
# time values in Factor Value / Characteristics columns are ordinary ISA metadata). The guard's purpose is to refuse a
# VALUE MATRIX masquerading as a sheet (such a table is almost entirely numeric). Factor Value[...] joins the permitted
# set; the refusal threshold is 50 % of the remaining cells. The header allowlist stays the primary defence: a value
# matrix fails it on its column names before the guard is consulted. The guard still refuses a sheet whose non-permitted
# cells are mostly numbers, and the observed fraction is reported for every accepted sheet.
NUM_OK = re.compile(r"^(Characteristics\[Age.*\]|Unit|Parameter Value\[.+\]|Factor Value\[.+\]|Comment\[.+\]|Term Accession Number)$")
NUM_GUARD_FRACTION = 0.50


def validate_investigation(text: str) -> dict:
    section, out = None, {"ok": True, "errors": [], "fields": {}}
    for ln in text.splitlines():
        if not ln.strip():
            continue
        if ln.strip() in I_SECTIONS:
            section = ln.strip(); continue
        if ln.startswith("#"):
            continue
        parts = ln.split("\t")
        if section is None:
            out["ok"] = False; out["errors"].append(f"line outside any section: {ln[:40]!r}"); break
        label = parts[0].strip().strip('"')
        if not I_FIELD_RX[section].match(label):
            out["ok"] = False; out["errors"].append(f"{section}: unknown field label {label!r}"); break
        if len(parts) > 51:
            out["ok"] = False; out["errors"].append(f"{section}: > 50 value columns"); break
        out["fields"].setdefault(label, []).append([p.strip().strip('"') for p in parts[1:]])
    for req in ("Study Title", "Study File Name", "Study Assay File Name"):
        if req not in out["fields"]:
            out["ok"] = False; out["errors"].append(f"missing {req}")
    return out


def _validate_table(text: str, header_rx: re.Pattern, mandatory: tuple[str, ...]) -> dict:
    lines = [ln for ln in text.splitlines() if ln.strip()]
    out = {"ok": True, "errors": [], "columns": [], "n_rows": 0, "rows": []}
    if not lines:
        out["ok"] = False; out["errors"].append("empty"); return out
    cols = [c.strip().strip('"') for c in lines[0].split("\t")]
    out["columns"] = cols
    bad = [c for c in cols if not header_rx.match(c)]
    if bad:
        out["ok"] = False; out["errors"].append(f"unexpected columns {bad[:5]}")
    for m in mandatory:
        if m not in cols:
            out["ok"] = False; out["errors"].append(f"missing {m}")
    if not out["ok"]:
        return out
    rows = [[c.strip().strip('"') for c in ln.split("\t")] for ln in lines[1:]]
    for i, r in enumerate(rows, start=2):
        if len(r) != len(cols):
            out["ok"] = False; out["errors"].append(f"row width mismatch at line {i}: {len(r)} fields vs {len(cols)} columns"); return out
    out["n_rows"] = len(rows)
    checked = guard = 0
    for r in rows:
        for c, v in zip(cols, r):
            if NUM_OK.match(c):
                continue
            checked += 1
            if NUMERIC.match(v):
                guard += 1
    out["numeric_fraction"] = round(guard / checked, 4) if checked else 0.0
    if checked and guard / checked > NUM_GUARD_FRACTION:
        out["ok"] = False; out["errors"].append(f"numeric guard: {guard}/{checked} numeric cells outside permitted columns (> {NUM_GUARD_FRACTION:.0%})")
        return out
    out["rows"] = [{norm(c): v for c, v in zip(cols, r)} for r in rows]
    return out


def validate_sample_sheet(text: str) -> dict:
    return _validate_table(text, S_HEADERS, ("Source Name", "Sample Name"))


def validate_assay_sheet(text: str) -> dict:
    return _validate_table(text, A_HEADERS, ("Sample Name",))


# ---------------------------------------------------------------- selftest
def selftest() -> bool:
    ok = True

    def check(c, msg):
        nonlocal ok
        print(("PASS " if c else "FAIL ") + msg); ok &= bool(c)

    rows = lambda n, person=False: [{"source name": f"SRC{i:03d}", "sample name": f"S{i:03d}", "characteristics[sex]": "F",
                                     **({"characteristics[subject id]": f"P{i:03d}"} if person else {})} for i in range(n)]
    cs = ["case-control"]
    d_t1 = "A case-control study. One blood sample was collected from each participant at enrolment."
    d_t2 = "Cross-sectional design. Each subject provided a single plasma specimen."
    d_t3 = "Samples were collected from 60 distinct patients with type 2 diabetes."
    d_t4 = "60 unique participants were recruited. Serum was analysed by LC-MS."
    check(subject_evidence(rows(60, True), cs, "")["grade"] == "R1", "R1: person key, one sample each")
    check(subject_evidence(rows(60), cs, d_t1)["grade"] == "R2_T12", "R2 T1 sentence -> R2_T12")
    check(subject_evidence(rows(60), cs, d_t2)["grade"] == "R2_T12", "R2 T2 sentence -> R2_T12")
    check(subject_evidence(rows(60), cs, d_t3)["grade"] == "R2_COUNT_SUPPORTED", "T3 with matching integer -> COUNT_SUPPORTED")
    check(subject_evidence(rows(60), cs, d_t4)["grade"] == "R2_COUNT_SUPPORTED", "T4 with matching integer -> COUNT_SUPPORTED")
    check(subject_evidence(rows(61), cs, d_t3)["grade"] == "SOURCE_UNIQUE", "T3 integer mismatch (60 vs 61 sources) -> SOURCE_UNIQUE")
    check(layer2("patients were not sampled once; one sample was collected from each participant twice a week", 60)["grade"] is None
          and layer2("one blood sample was not collected from each participant", 60)["grade"] is None, "negation and frequency guards block")
    check(layer2("the study enrolled 100 distinct patients", 100)["grade"] == "R2_COUNT_SUPPORTED"
          and layer2("the study enrolled 100 distinct patients", 99)["grade"] is None, "T4b (verb-first count) with exact integer only")
    check(layer2("a subset of participants provided a single sample", 60)["grade"] is None, "subset wording blocks")
    check(layer2("one sample was collected from each participant. serial samples were collected at each visit", 60)["grade"] is None,
          "positive sentence does not override repeated-sampling sentence")
    check(layer2("the study enrolled 60 distinct patients; normal controls were included", 60)["grade"] == "R2_COUNT_SUPPORTED"
          and not any(g.search("normal") for g in NEGATION), "'no' does not match inside 'normal'")
    check(subject_evidence(rows(60), ["longitudinal cohort"], d_t1)["grade"] == "SOURCE_UNIQUE", "Layer 1 negative descriptor blocks R2")
    check(subject_evidence(rows(60), ["cohort study"], "One blood sample was collected from each participant.")["grade"] == "SOURCE_UNIQUE",
          "'cohort' alone is not positive even with a T1 sentence")
    check(subject_evidence(rows(49), cs, d_t1)["grade"] == "EXCLUDED", "49 subjects -> EXCLUDED")
    dup = rows(60); dup[5]["source name"] = "SRC004"
    check(subject_evidence(dup, cs, d_t1)["grade"] == "EXCLUDED", "source with 2 samples -> EXCLUDED")
    sen = rows(60); sen[7]["source name"] = "NA"
    check(subject_evidence(sen, cs, d_t1)["grade"] == "UNRESOLVED", "sentinel source -> UNRESOLVED")
    pm = rows(60, True); pm[3]["characteristics[subject id]"] = "unknown"
    check(subject_evidence(pm, cs, "")["grade"] == "UNRESOLVED", "missing person-key value -> UNRESOLVED")
    check(subject_evidence(rows(60), cs, "Fifty-nine women and one man took part.")["grade"] == "SOURCE_UNIQUE", "no template -> SOURCE_UNIQUE")
    check(temporal_class("collectionWeek") == "strong" and temporal_class("AFTER_MEAL_TIME") == "ambiguous" and temporal_class("Visit") == "strong"
          and temporal_class("Cancer stage") == "ambiguous" and temporal_class("Diagnosis") == "none", "temporal detector handles camel/underscore")
    inv = "INVESTIGATION\nInvestigation Identifier\tMTBLSX\nSTUDY\nStudy Title\tT\nStudy File Name\ts_X.txt\nSTUDY ASSAYS\nStudy Assay File Name\ta_X.txt\n"
    check(validate_investigation(inv)["ok"], "investigation validator accepts vertical sections")
    check(not validate_investigation("Study Title\tT\n")["ok"], "investigation validator refuses line outside sections")
    s_ok = "Source Name\tCharacteristics[Organism]\tProtocol REF\tSample Name\n" + "\n".join(f"SRC{i}\tHomo sapiens\tSample collection\tS{i}" for i in range(5))
    check(validate_sample_sheet(s_ok)["ok"], "sample validator accepts ISA headers")
    check(not validate_sample_sheet("Source Name\tSample Name\tIntensity\nA\tS\t3.2")["ok"], "sample validator refuses unknown header")
    check(validate_sample_sheet("Source Name\tMaterial Type\tDescription\tProtocol REF\tSample Name\nSRC\twhole organism\tx\tP\tS1")["ok"], "A-4: ISA-Tab sample-table attribute columns accepted")
    s_num = "Source Name\tSample Name\t" + "\t".join(f"Characteristics[m{j}]" for j in range(10)) + "\n" + "\n".join(f"SRC{i}\tS{i}\t" + "\t".join(str(3.1 + i + j) for j in range(10)) for i in range(5))
    v = validate_sample_sheet(s_num); check(not v["ok"] and "numeric guard" in v["errors"][0], "A26: value-matrix-like sheet (83 % numeric outside permitted columns) refused")
    s_meta = "Source Name\tSample Name\tCharacteristics[Height]\tFactor Value[dose]\tCharacteristics[Organism]\n" + "\n".join(f"SRC{i}\tS{i}\t{170 + i}\t{10 * i}\tHomo sapiens" for i in range(5))
    v = validate_sample_sheet(s_meta); check(v["ok"] and v["numeric_fraction"] == 0.25, "A26: ordinary numeric metadata (height, dose) accepted; numeric fraction reported (0.25)")
    a_ok = "Sample Name\tProtocol REF\tMS Assay Name\tMetabolite Assignment File\nS1\tExtraction\tA1\tm_X_maf.tsv"
    check(validate_assay_sheet(a_ok)["ok"], "assay validator accepts ISA headers incl. MAF filename column")
    a_r3 = "Sample Name\tProtocol REF\tNMR Assay Name\tNormalization Name\tData Transformation Name\tMetabolite Assignment File\nS1\tNMR\tA1\tN1\tDT1\tm_X_maf.tsv"
    check(validate_assay_sheet(a_r3)["ok"], "A-4: ISA-Tab process-node columns (Normalization Name, Data Transformation Name) accepted")
    check(not validate_assay_sheet("Sample Name\tIntensity\nS1\t3.2")["ok"] and not validate_assay_sheet("Sample Name\tNormalisation Name\nS1\tN1")["ok"],
          "A-4: non-ISA assay column still refused (incl. a near-miss spelling)")
    a_fv = "Sample Name\tFactor Value[dose]\n" + "\n".join(f"S{i}\t{10 * i}" for i in range(5))
    check(validate_assay_sheet(a_fv)["ok"], "A26: numeric Factor Value column in an assay sheet accepted (permitted set)")
    a_mx = "Sample Name\t" + "\t".join(f"Characteristics[m{j}]" for j in range(10)) + "\n" + "\n".join(f"S{i}\t" + "\t".join(str(i + j) for j in range(10)) for i in range(5))
    check(not validate_assay_sheet(a_mx)["ok"], "A26: value-matrix-like assay table still refused")
    inv_bad = "INVESTIGATION\nInvestigation Identifier\tX\nSTUDY\nStudy Title\tT\nUnexpected Measurement Data\t1\t2\nStudy File Name\ts.txt\nSTUDY ASSAYS\nStudy Assay File Name\ta.txt\n"
    v = validate_investigation(inv_bad)
    check(not v["ok"] and "unknown field label" in v["errors"][0], "investigation validator refuses unknown field label inside a known section")
    inv_term = "STUDY\nStudy Title\tT\nStudy File Name\ts.txt\nSTUDY DESIGN DESCRIPTORS\nStudy Design Type\tcase-control\nStudy Design Type Term Accession Number\tOBI:1\nStudy Design Type Term Source REF\tOBI\nComment[x]\ty\nSTUDY ASSAYS\nStudy Assay File Name\ta.txt\nStudy Assay Technology Type\tmass spectrometry\n"
    check(validate_investigation(inv_term)["ok"], "investigation validator accepts ontology variants and Comment[]")
    s_short = "Source Name\tSample Name\tCharacteristics[Organism]\nSRC1\tS1\tHomo sapiens\nSRC2\tS2\n"
    v = validate_sample_sheet(s_short)
    check(not v["ok"] and "row width mismatch at line 3" in v["errors"][0], "sample validator refuses short row (zip truncation closed)")
    s_long = "Source Name\tSample Name\nSRC1\tS1\tEXTRA\n"
    check(not validate_sample_sheet(s_long)["ok"], "sample validator refuses long row")
    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "selftest":
        sys.exit(0 if selftest() else 1)
    sys.exit(__doc__)
