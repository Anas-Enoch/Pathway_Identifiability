# Stage A — External-data feasibility protocol (value-blind cohort recovery)

**Status: v0.7.8 — census COMPLETE to a terminal status for 3557/3582 accessions; strict and sensitivity L1 = 0;
25 operational failures await a final retry invocation under the v0.7.7 lock; diagnostic rescreen of 230 "all rows QC"
studies proposed (A34); nothing authorised by this version.** 2026-10-28. Supersedes v0.7.7
(`f37e94c749f64072ba9c883f3ee7e7534d1d55dd5b7b51de849d7ab2657b9922`). Record `STAGEA_CENSUS_RESULT_v0.9_PROVISIONAL.md`.
Change, census v0.1.8 (64 offline checks; `stagea_isa` v0.2.2 unchanged; classification identical to v0.1.7 — reporting
and scope only): organism values, person-key-like columns and `qc_triggers` (designated column : vocabulary token, never a
cell value) are recorded for every study **before** any M2/M3/M4 return; `--rescreen-allqc` re-screens exactly the records
with status UNRESOLVED, rule M2, reason "no biological rows after QC removal", under the same rules, journaled under its
own header (`rescreen_allqc=1`). No D1, M1–M9, A22, A26, A27 or A30 rule changes; any rule change suggested by the
triggers is a separate dated amendment (A35) decided after the diagnostic, not before.

v0.7.7 had superseded v0.7.6 (`dc4e313c83a76ebaf273742d7ffba44b67dcaac078cc2fcd17116b954384b117`).
Record `STAGEA_CENSUS_RUN1_INTERIM_v1.0.md`. Changes, census v0.1.7 (60 offline checks; `stagea_isa` v0.2.2 unchanged):
**D-full-1** outage stop — 5 consecutive studies ending `TRANSPORT_FAILED` stop the run (exit 2, `halted_outage`);
**D-full-2** index hrefs decoded once and encoded once (v0.1.6 double-encoded → 9 false 404s); **D-full-3** an HTTP error
on an index-listed file is `FETCH_HTTP_<code>`, not `ACCESS_REFUSED`; reporting fields `census_version`,
`person_key_columns`; `design_evidence` label corrected. **A20 retry scope (amended):** the retry pass re-requests
`TRANSPORT_FAILED`, `FETCH_HTTP_*`, and v0.1.6 `ACCESS_REFUSED` records carrying the double-encoding signature (`%25`, not a
redirect/binary/oversize refusal); it may span several invocations if an outage stop interrupts it, each journaled under
its own header. **Run-1 completion rule:** the retry pass for run 1 executes under the v0.7.7 lock against the same output
directory; the run-1 lock and authorisation are preserved beside it; v0.1.6 records completed in run 1 stand (v0.1.7
classifies them identically). No rule of D1, M1–M9, A26, A27 or A30 changes.

v0.7.6 had superseded v0.7.5 (`1238401d852f6d05ea6209bf18da4a5d385aba1ab81f098cb89309b4fc26dc98`).
Change (census v0.1.6, 55 offline checks; `stagea_isa` v0.2.2 unchanged): **A30 evidentiary clarification (owner).**
The evidence that a row is an analytical reference is the **explicit designation in the organism cell** (`reference
compound`, `reference material`, `reference standard`, `blank`, `qc`, `quality control`, `standard mixture`,
`calibrant`), accepted as QC provided no other metadata conflicts with it. A missing taxonomy accession is **not**
corroboration; it is the absence of a conflicting annotation. An explicit analytical label that carries a taxonomy
accession is a conflict: the row is QC only if another designated column affirms the analytical role; otherwise the
study is UNRESOLVED (`n_rows_organism_analytical_conflict`). Missing organism without QC evidence stays UNRESOLVED.
Classes renamed accordingly (`ANALYTICAL`, `ANALYTICAL_CONFLICT`); the boolean outcome for MTBLS27-shaped rows is
unchanged, its justification is not. §3.7: `DOC_CONTRADICTED` contradicts the claimed independent-sampling structure of
the chosen assay; whether a separately defined subset remains usable needs its own evidence — it does not declare
every observation unusable. Post-retry pilot 25b record pending the updated archive (not reconstructed from a summary).

v0.7.5 had superseded v0.7.4
(`80c8e3a5c7fb74bc8d85b4042e588ee570fc2f98b0c6416dfbd59c6e6131c53b`). Records: `STAGEA_PILOT25B_RESULT_v1.0.md`
(`cf529b2e…e644`). Changes, census v0.1.5 (54 offline checks; `stagea_isa` v0.2.2 unchanged): **A30** —
`Characteristics[Organism]` enters QC classification with the owner's distinctions (§3.5 M1/M2): an analytical
annotation (`reference compound|reference material|reference standard|blank|qc|quality control|standard mixture|
calibrant`) makes the row QC only when *supported* — another designated column carries a QC token, or the organism
cell carries no taxonomy accession; an analytical annotation that carries a taxonomy accession and no other QC
evidence, or a missing organism without QC evidence, makes the study **UNRESOLVED** (M1), never QC by default and never
EXCLUDED; an explicit non-human organism remains EXCLUDED; the raw organism values and the counts of rows removed /
missing / unsupported are preserved in every record. The cell's text takes precedence over a template-copied accession.
**D1 (owner, 2026-10-27):** unique `Sample Name` values are not evidence of participant independence; `SOURCE_DEGENERATE`
stays UNRESOLVED (an identifiability limitation of the deposited metadata, not a biological exclusion); D1 unchanged.
**§3.7 (new, development only):** a separate documentary-recovery route for degenerate-source studies, reported as
`L1_documentary_recovered`, never added to strict or sensitivity L1. **§1 limitation (new):** the strict census measures
E = 1 (evidence satisfying D1), not I = 1 (true independence); P(E = 1 | I = 1) < 1, so the strict count is a lower bound
on truly eligible cohorts only if the classifier makes no false-positive independence assignment — a goal, not a proven
property; the manuscript claim is "no cohort met the prespecified affirmative independence-evidence criteria among the
evaluated accessions", never "no eligible cohort exists". **F5:** the MTBLS8 retry pass is approved under the v0.7.4
pilot lock (outputs copied first; journal must show the original 111 records plus the retry records).

v0.7.4 had superseded v0.7.3 (`ff1f7308f3421c5a7984704da32d6423132f907850cb03f29ea522f5a1074db3`). Census v0.1.4 (50 offline checks);
`stagea_isa` v0.2.2 unchanged (35). Changes: **provenance** — the fetch log is a journal written by the transport as each
request completes (append + flush, run header per invocation), so an abrupt termination loses at most the request in
flight; tested by a simulated termination (journal holds the four completed requests, checkpoint the one completed
study). **C5** — the M2 diagnostic is a count (`n_rows_unidentified`); no sample identifier is written to any public
table; tested. **C6** — elapsed seconds (monotonic) on every journal line, successes and failures alike; tested.
**C7** — summary keys `D3_metadata_units`, `D3_metadata_units_strict`, `D3_metadata_units_by_stratum`,
`L3_confirmed_status = NOT_ASSESSED`; `L3_units` no longer appears. A26 (numeric guard: `Factor Value[…]` permitted,
50 % threshold, fraction reported; an engineering threshold, not a validated boundary) and A27 (`SOURCE_DEGENERATE` →
UNRESOLVED, outside both L1 counts; multiplicity across several sources stays EXCLUDED) unchanged from v0.7.3. **A28**:
same 25 accessions (MTBLS3 … MTBLS28), new lock, `pilot=25`, output `stageA_pilot25b/`; the record compares the two
pilots study by study (format refusals, exclusions, unresolved incl. `SOURCE_DEGENERATE`, strict/sensitivity L1,
transport failures, request provenance preserved). The freeze record v1.1 binds these bytes, not v0.7.2's.

v0.7.3 had superseded v0.7.2 (`c464a603819be9bae7c33ae8487fd942231bcbaadee6c426e66b09dbcbde0984`).
Record: `STAGEA_PILOT25_RESULT_v1.0.md`. Changes: **P1 (defect, census v0.1.3):** the fetch log is appended per
invocation with a run header, never rewritten (the pilot's resume check had overwritten the 112-request log; per-file
SHA-256s survive in `stageA_M_files.csv`); the summary carries `n_requests_this_invocation` and `n_requests_cumulative`.
**A26 (proposed; `stagea_isa` v0.2.2):** sample-/assay-sheet numeric guard recalibrated — `Factor Value[…]` permitted;
refusal above 50 % numeric cells outside the permitted columns (observed 6–22 % on real sheets; a value-matrix-like table
at 83 % still refused); the fraction is reported per sheet. The header allowlist remains the primary defence. The archived
schema probe v0.3.7 stays bound to isa v0.2.1; isa v0.2.2 is for the census only. **A27 (proposed; census v0.1.3):**
a single `Source Name` for all samples (`n_sources == 1`, `n_samples > 1`) is `UNRESOLVED` (`SOURCE_DEGENERATE`), not
EXCLUDED; `n_sources > 1` with multiplicity remains EXCLUDED. This amends D1's "a source with > 1 sample → EXCLUDED" for
the degenerate case only and is the owner's decision. Also: rows without a `Sample Name` are counted
(`n_rows_without_sample_name`; outcome unchanged, UNRESOLVED under M2). Census v0.1.3, 46 checks; isa v0.2.2, 35 checks.

v0.7.2 had superseded v0.7.1 (`8a13625822c843f8ed4fb9e8419467c0b15e1867e548d8e50a10db7b421ef8eb`) after the
owner's A24 audit (conditional approval; corrections C1–C3; the owner withdrew the "demonstrated defect" classification
of the assay-name reconciliation and recorded it as a prospective robustness improvement). Changes, census v0.1.2
(41 offline checks; `stagea_isa` v0.2.1 unchanged): **C1** the run's terminal exit path returns 1 whenever the summary
carries `halted_protocol_violation` (previously 0); tested on the same function the command calls. **C2** D3 CONFIRMED
merging requires overlap of *specific* person-identifier hashes with the same affiliation; same affiliation + title
similarity ≥ 0.6 is now a `SUSPECTED_OVERLAP` review flag (kind `title`), never a merge; tested. **C3** identifier hashes
are documented as deterministic pseudonymisation (constant prefix; low-entropy identifiers are dictionary-recoverable),
live only in the checkpoint, which the summary lists under `artefacts_not_for_publication`; no published table carries
a hash. No change to A1–A23, the evidence grades, QC, assay selection or the pilot rule. **A24:** N = 25, the first 25
non-exemplar accessions of the manifest in numeric order (MTBLS3 … MTBLS28); operational pilot, not a prevalence
estimate; evaluation criteria in §7.4c. The 7.4 freeze record (`STAGEA_FREEZE_RECORD_v1.0.md`, unsigned draft) binds
this file, census v0.1.2, `stagea_isa` v0.2.1, manifest `3e483c7a…9dc5`, N = 25 and A1–A24.

v0.7.1 had superseded v0.7.0
(`855ed2434b5c612721f71563c4d3c1be5487e5063c4f042800edcc6e76b60cb4`) after the owner's A18 audit (conditional approval).
Changes (§7.4a, §7.4b): (1) **assay-name reconciliation** (owner's blocking finding) — investigation-listed and
mirror-listed assay names are reconciled deterministically: exact name; else a unique unmatched pair is paired and
recorded as such; anything else stays unmatched and a study whose only assays are unmatched is `UNRESOLVED_STRUCTURE`,
never `EXCLUDED`; both names are preserved in the output; regression tests use the actual MTBLS1/MTBLS2 name pairs.
A factual qualification is recorded: in run 3 the exemplars' own investigation files listed the **mirror** names, so the
discrepancy observed live was REST-vs-mirror, not investigation-vs-mirror; the correction guards against the latter
case, which has not been observed but cannot be excluded; (2) **A20** three failure classes — operational
(`TRANSPORT_FAILED`), structural (`FORMAT_REFUSED`, `UNRESOLVED_STRUCTURE`), access anomaly on the server side
(`ACCESS_REFUSED`: redirect, binary, oversize) — all per-study; a **guard refusal (protocol violation on our side) halts
the whole census** and is recorded in the summary; (3) **A21** pilot mechanism: the authorisation line may carry
`pilot=N`; selection rule fixed as the first N accessions of the manifest in numeric order; a full run is refused while
the line carries `pilot=`; (4) **A22** two-tier QC: strong analytical tokens remove the row anywhere in the designated
columns; ambiguous tokens (`reference`, `standard`, `std`, `pool(ed/s)`) remove the row only in a sample-describing
column (`Sample Name`, `Material type`, `Sample type`) and make the study `UNRESOLVED` (M2) when they appear in a
`Factor Value` or `Comment`; (5) **A23** more than one sample sheet → `UNRESOLVED_STRUCTURE` (counted; not a biological
exclusion); (6) **D3** hashed identifiers (SHA-256 of the normalised identifier with a public constant prefix; raw
identifiers never written) retained for eligible cohorts with provenance source/sample/person and a generic-identifier
test; CONFIRMED merge only on specific person-identifier overlap with the same affiliation, or the documentary rule;
every other overlap is a `SUSPECTED_OVERLAP` review flag in `stageA_M_overlap_flags.csv`, never a merge. Census
v0.1.1, 38 offline checks; `stagea_isa` v0.2.1 unchanged.

v0.7.0 had superseded v0.6.9
(`6e6b095bad1efe8ee8eb33b1845abc983a6ca2df40314d1e4a7ab8c51c63b2e6`). Changes: (1) §7.5b records the A17 result as
reported by the owner: LIST and mirror repeats stable; REST 3584 ⊂ mirror 3596; twelve mirror-only accessions
(MTBLS9, 30, 202, 224, 351, 539, 557, 643, 728, 926, 4938, 14881) carried UNRESOLVED, never screened without a dated
amendment; census universe = intersection 3584; screening pool after A7 = **3582**; (2) **§7.4a — A18 census
executable** `stageA_census.py` v0.1.0 (27 offline checks on synthetic studies incl. an MTBLS1-layout fixture; no live
request has been made by it) with the execution decisions A20–A23 the owner must take before any lock; (3) the
authorisation statement for the census carries the **manifest SHA-256** so the lock binds the universe. Probe v0.3.7,
`stagea_isa` v0.2.1 and audit v0.1.1 unchanged.

v0.6.9 had superseded v0.6.8
(`b1320f7688029b1e97ff4a007b95cf579196a2ca252018f64c7eff5f6fa29fc6`) by **Amendment A-5 (audit execution integrity; no
verdict-rule change)** after A17 run 1 (`STAGEA_A17_RUN1_RESULT_v1.0.md`): LIST ×3 byte-identical (3584, no pagination key
or header); mirror root index 721 989 B streamed in ~3 min; the second root read timed out at 120 s → recorded halt →
UNRESOLVED. Two corrections, both prospective: (1) the root index's parent-directory link (`…/studies/`) had been
counted as a non-conforming token — absolute and parent hrefs are navigation, listed and ignored; (2) FTPROOT per-read
timeout 600 s, no retry; plus (3) on a halt after one LIST and one root success the comparison and manifest are written
marked `PROVISIONAL_` (verdict unchanged; the census refuses a provisional manifest) and (4) elapsed time per request is
logged. Provisional evidence, not a verdict: mirror 3596 vs REST 3584, same maximum MTBLS16154 → [Likely] `BOUNDED`,
U = 3584, band [3584, 3596], twelve accessions to be named by run 2. Audit v0.1.1, 26 offline checks. Also recorded:
§7.1f now points to `STAGEA_RUN3_RESULT_v1.1.md` (per-file hashes; R1 person key absent on both exemplars; 86-label
template observation). The schema probe and `stagea_isa` are unchanged.

v0.6.8 had superseded v0.6.7
(`ca15b8e0ce5c9f6f9267fc734c86cb78ef6d48c39eb9005110ed63d221283c7d`). Changes: (1) §7.1f records the owner's close-out
of the schema-inspection phase and the two exemplars' findings as reported (MTBLS2: *Arabidopsis thaliana*, 16 sample /
16 assay rows, genotype and replicate factors → not human; MTBLS1: human status and person identifier not demonstrated;
L2-confirmed NOT_ASSESSED for both; both excluded from the census under A7); (2) **§7.5a — A17 Public-universe
verification**, a separately locked, metadata-only audit with fixed verdict rules, delivered as `stageA_universe_audit.py`
v0.1.0 (23 offline checks; loopback-proved redirect refusal and header capture); (3) the authorisation statement for A17
is its own file `stageA_universe_authorisation.txt` with its own line. The schema probe (v0.3.7) and `stagea_isa`
(v0.2.1) are **not modified** by this version; their run-3 lock stands as archived.

v0.6.7 had superseded v0.6.6
(`d80bebdde10bd02017923063d32909e8909d2205afeac7316f3a69a081ec65f1`) by **Amendment A-4 (validator vocabulary only)**:
the `a_MTBLS1_metabolite_profiling_NMR_spectroscopy.txt` fetched from the mirror (200, 120 888 B, SHA-256
`cd615abd…de7a`) was FORMAT-REFUSED for the ISA-Tab 1.0 process-node columns `Normalization Name` and
`Data Transformation Name`, which the v0.2.0 assay allowlist omitted. A-4 completes the assay-table and sample-table
allowlists of §3.3 to the ISA-Tab 1.0 column vocabulary **from the specification** (assay: `Hybridization Assay Name`,
`Scan Name`, `Normalization Name`, `Data Transformation Name`, `Array Design REF`, array/image/spot-picking/protein/
peptide/PTM file columns, `Characteristics[…]`, `Factor Value[…]`, `Material Type`; sample: `Material Type`,
`Description`, `Provider`, `Performer`, `Date`); non-ISA and near-miss headers still refuse; the numeric guard, the
investigation allowlist, every route, guard and the fetch set are unchanged. `stagea_isa` v0.2.1 (33 checks), probe
v0.3.7 (97 checks). Record: `STAGEA_RUN3_RESULT_v1.0.md`. Attempt-2 outputs in `stageA_out3b/` are preserved; attempt 3
uses `stageA_out3c/`.

v0.6.6 had superseded v0.6.5
(`74adde21c84c9c0b7c2779f3f1ceb4fe7ff4d39feeab650c58ada62f4c98c82d`) by **Amendment A-3 (execution integrity only; no
rule, route, guard or fetch-set change)**: on the first request (LIST) the remote end closed the connection without a
response (`http.client.RemoteDisconnected`); that exception class is not a `URLError`, so the v0.3.5 transport branch did
not convert it into a recorded halt — the process died with a traceback, exit 1 by crash, and **no report and no fetch
log were written**. A-3: every transport-level failure (`URLError`, `http.client.HTTPException`, timeout, `OSError`) is a
recorded Halt — logged as `ERR` with the exception class, `endpoint_reachable = false`, report and fetch log written,
exit 1 — with offline tests for `RemoteDisconnected`, `OSError` and timeout and a loopback proof (server closed the
socket without a response; halt recorded, outputs written). Nothing is retried by the code; a transport failure is
reported, not worked around. Probe v0.3.6, 93 offline checks (29 from `stagea_isa` v0.2.0, unchanged). The v0.6.5 lock
(signed `Anas 2026-10-25T20:50:02Z`) retrieved zero bytes; its directory `stageA_out3/` holds only the signed lock and the
traceback log and is preserved as the record of attempt 1. Attempt 2 uses a new lock and `stageA_out3b/`.
[Guessing] cause: a transient server-side close or a connection-rate limit after several same-day runs; no evidence either
way, and the protocol does not act on the guess.

v0.6.5 had superseded v0.6.4 (`a71e5113a4e863c1e63427923a786fe6fc4eb2019174d5210545a51147f0ee08`) by **Amendment A-2**,
supported by observed interface evidence from run 2 (`STAGEA_RUN2_RESULT_v1.0.md`): on both exemplars the FTP mirror index
was `INDEX_ACCESSIBLE`, `i_`/`s_` were present under the REST names, and the **assay-sheet name differed** (REST
`a_MTBLS1_NMR_metabolite_profiling_NMR_spectroscopy.txt` vs mirror `a_MTBLS1_metabolite_profiling_NMR_spectroscopy.txt`;
REST `a_MTBLS2_LC-MS_…` vs mirror `a_MTBLS2_metabolite_profiling_mass_spectrometry.txt`). A-2: (1) the per-exemplar FETCH
set is the **mirror index's own** `i_/s_/a_` names under the same guards; REST FILES names serve the comparison table and
forbidden-name counting only; a REST name absent from the mirror is recorded (`rest_names_not_fetched`) and never
requested; (2) an index class other than `INDEX_ACCESSIBLE` halts with "fetch set undefined" — no fall-back to REST names;
(3) a mirror index lacking an `i_` or an `s_` halts; (4) mirror entries that are neither metadata nor forbidden are listed
by name, never fetched; (5) the fetch set is held outside the report (no transient key survives a halt); (6) the single
FTP template and every guard of v0.6.4 are unchanged; no second route (A13-b stands). Source of record from run 3 on: the
**public release** on the mirror, not the REST working copy. Probe v0.3.5, 90 offline checks (29 from `stagea_isa` v0.2.0,
unchanged). Expected outcome of run 3: FETCH of six files (3 per exemplar) and per-class validation; any refusal halts.

v0.6.4 had superseded v0.6.3 (`3fbec18a70e412c645498c63f2935dcb1240709a0216afb0ab3930c51a1dad44`) by the owner's **A13
decision: A13-a APPROVED (diagnostic-only run 2); A13-b REJECTED (no unverified second FETCH route)**. Changes from v0.6.3:
(1) the guessed web FETCH template is **removed**; the single FETCH route is the public FTP mirror template of v0.6.2
(hosts restricted to `ftp.ebi.ac.uk` for FETCH); (2) the FTP directory index (`FTPDIR`, names-only) is classified into
**fixed diagnostic classes decided before execution** (§3.3a): `INDEX_ACCESSIBLE` (names compared with the FILES listing),
`DIRECTORY_404` (directory template unresolved; mirror contents not established), `DIRECTORY_403` (index inaccessible;
file existence unresolved), `DIRECTORY_OTHER_<code>` (unresolved), redirect → halt (destination never requested),
unrecognised index structure → halt (no names inferred); an unavailable index is never reported as an absence of the
file; (3) the inspection is **two-phase**: FILES + FTPDIR for both exemplars are recorded before any metadata FETCH, so
a FETCH halt on the first exemplar cannot suppress the second exemplar's diagnosis; (4) the per-study record is attached
to the report before any request for that study, so a halt preserves what was established (in v0.3.3 and earlier a halt
discarded the partial study record — a defect of the reporting, not of the access control; run 1's stage flags were
unaffected). Guards retained unchanged: redirect refusal, pre-read binary/bundle guard, size cap, two-exemplar limit,
forbidden-name counting. Probe v0.3.4, 89 offline checks (29 from `stagea_isa` v0.2.0, unchanged). The v0.6.3 bundle is
superseded and was never locked; the archived v0.6.2 run 1 outputs are unchanged.
Expected outcome of run 2: a FETCH on the FTP template will [Likely] 404 again and halt with exit 1 — that is the
designed result; the deliverable is the FTPDIR class and name comparison for MTBLS1 and MTBLS2.

v0.6.3 had superseded v0.6.2 (`65ebae7e1a8b1dfc21d61d8d156a6bebc27a42f3afafc89f9313660f8766791c`) by **Amendment A-1** after the
owner's signed run 1 (§7.1a): LIST and FILES verified live (200, recognised shapes); FETCH template 1 (public FTP mirror)
returned **404** for the first listed MTBLS1 assay file; halt; exit 1; nothing forbidden requested. A-1 changes only the
metadata-retrieval route and adds guards; no rule, grade, threshold, validator or exemplar rule changes:
(1) an **ordered, fixed list of FETCH candidates** — template 1 unchanged, template 2 the per-file web route
`https://www.ebi.ac.uk/metabolights/{MTBLS}/files/{name}` — tried in that order per file; every attempt logged;
all-fail halts; a redirect, oversize or binary response halts immediately without falling through; the REST
`…/studies/{MTBLS}/download` family is **excluded** as a bundled-download endpoint (§3.3). Template 2 is UNVERIFIED and
was chosen from documentation memory, not observation (A13). (2) a **names-only directory index** of the FTP study folder
(`FTPDIR`; hrefs parsed, nothing fetched) recorded per exemplar to distinguish a wrong path from a wrong file name;
its unavailability is reported, not fatal. (3) a **pre-read binary/bundle guard** on every FETCH: `Content-Length` above
the size cap or a zip/gzip/octet-stream `Content-Type` halts before the body is read (loopback-proved: the server's
oversize write was reset unread); a body with a zip/gzip signature halts after reading and is never parsed. (4) the
authorisation line names this file. Probe v0.3.3, 83 offline checks (incl. 29 from `stagea_isa` v0.2.0, unchanged).
Listing completeness after run 1: **UNRESOLVED** (3584 accessions returned; no total field) — §7.5 census remains
blocked; A7 exemplars (MTBLS1, MTBLS2) are the minimum of what was returned, not of the repository.
v0.6.2 (NOT FROZEN, NOT EXECUTED at the time) had superseded v0.6.1 (`61b6d30a…bb33`) after the owner's
A12 pre-authorisation audit: **HTTP redirects are never followed** — the live opener's redirect handler returns None, so
any 301/302/303/307/308 is a recorded halt (status and Location text logged; destination never requested); the offline
mock models a 302 to an `m_*.tsv` and asserts the destination is absent from the fetch log; a loopback HTTP server test
confirmed the live opener stops at the 302 (server saw one request). A redirect destination enters the allowlist only by a
prospective amendment. Probe v0.3.2, 76 offline checks.
v0.6.1 had superseded v0.6 (`f3248c39…1095`) after the owner's
pre-lock audit of probe v0.3.0: the authorisation file is now validated as **exactly one statement** (full-file match after
stripping surrounding whitespace; any additional line, before or after, rejects), with regression tests; probe v0.3.1,
72 offline checks. Owner's qualification recorded: listing completeness, if UNRESOLVED, is acceptable for the schema
inspection but **blocks authorisation of the Stage A census** until pagination/completeness handling is resolved (§7.5).
v0.6 had superseded v0.5 (`5ef5aca6…a2f9`) after the owner's audit
of probe v0.2.0 (milestone approved; A12 HOLD for two execution-integrity corrections): (1) a halted inspection now writes
its report and fetch log and **exits nonzero** (regression-tested); (2) authorisation is no longer a source constant —
it is a separate owner-written file `stageA_live_authorisation.txt`, hashed **into** the lock as the fourth hash, so no byte
of source or protocol changes between lock and run (the v0.5 sequence would have refused itself); (3) the schema report
separates endpoint reachable / response shape recognised / metadata validated / listing completeness (ESTABLISHED only
when a total field equals the number returned; otherwise PAGINATED_UNRESOLVED or UNRESOLVED, with the consequence for
A7 stated). Probe v0.3.0, 69 offline checks.
v0.5 had superseded v0.4 (`849b94ad…5790`) after the owner's code
audit (offline implementation verified; HOLD on schema execution): the schema probe now implements the **complete**
inspection workflow (LIST → exemplars → FILES → guarded FETCH of i_/s_/a_ → per-class validation → report + fetch log),
tested end-to-end on an offline mock including malformed JSON, unexpected shapes, non-200, oversize, forbidden names in
listings (counted, never fetched), unknown investigation fields, malformed rows and numeric injection; the investigation
validator now allowlists every field label per ISA section (unknown label → FORMAT-REFUSED); table validators refuse any
row whose width differs from the header (zip truncation closed); the lock binds protocol + both sources + endpoint
templates + authorisation conditions; live execution is gated last by `LIVE_AUTHORISED = False`, flipped only by a
dated amendment recorded here. Endpoint templates remain **UNVERIFIED**. The regex recogniser is tested against its
synthetic suite, not proven against arbitrary natural language; that limitation is permanent and stated in §2.
v0.4 had superseded v0.3 (`dd821bb5…2249`) after the owner's A10/A11
decisions: evidence grades separated (R1, R2_T12, R2_COUNT_SUPPORTED, SOURCE_UNIQUE, UNRESOLVED, EXCLUDED); T3/T4 count
statements yield `R2_COUNT_SUPPORTED`, reported as a **sensitivity** L1 count, never the strict one; A11 exact integer equality
with no tolerance; `L2_confirmed_status = NOT_ASSESSED` in Stage A; implementation delivered (`stagea_isa.py` v0.1.0,
`stageA_schema_probe.py` v0.1.0) with offline selftests; schema-lock manifest format fixed.
v0.3 had superseded v0.2 (`9b1430d9…a335`) after the
owner's review (architecture accepted; A8 approved with a stronger documentary requirement; A9 approved with
conditions; L2 split). v0.3 changes: (1) R2 is a **two-layer** rule — design classification with frozen
descriptor sets, then an affirmative one-person-one-source statement matched by exact, negation-guarded
templates tied to the source count (§2, A8); (2) `CONFIRMED-with-note` counts toward L1 only when D1 is already
established by R1/R2 (A9); (3) L2 is split into L2-preliminary (Stage A, repository count) and L2-confirmed
(Stage B/C, usable matrix), with L2-confirmed = 0 by construction in Stage A; (4) design evidence and
individual-subject evidence are reported as separate columns everywhere.
v0.2 had changed from v0.1: separate validators per ISA file class; `SOURCE_UNIQUE` ≠ subject independence; two-level
temporal detector; schema lock before the Stage A freeze; A1–A7; three feasibility levels; Stage B rule for
`FEATURE_COUNT_UNKNOWN`.
Authorised: drafting only. **Not authorised:** any schema fetch, any MetaboLights or untouched mwTab retrieval,
any census execution, any measured-value access, the internal 93-study assessment, changes to P4V v0.12 or to the
six frozen programmes.

**Chronology statement.** Written after B0, the E4 adjudication and the B3 audit showed that Workbench REST
`factors` carries no subject identifier for 124 candidates, and before any MetaboLights study metadata or untouched
mwTab file was retrieved. The rules below are informed by that failure mode, not by any candidate they will screen.

## 1. Objective and the three feasibility levels

Stage A reports, from metadata alone and for cohorts never opened by this programme, counts at three levels, each
as **confirmed / unresolved / excluded**, and each broken down by stratum (repository, modality MS/NMR, platform
class, biological matrix):

| Level | Question | Confirmed means |
|---|---|---|
| L1 subject independence | does the metadata affirmatively establish ≥ 50 independent human individuals, one sample each, in the chosen assay? | D1 satisfied by route R1 or R2 (§2) |
| L2-preliminary task feasibility | does a validated repository-level count indicate ≥ 60 named features in the chosen assay? | D4(a) satisfied; a proxy, not a usable matrix |
| L2-confirmed task feasibility | does the assay yield a usable subject-by-feature matrix for the frozen task (≥ 60 columns after Stage C's eligibility rules; design/held-out split; prediction targets = metabolite columns)? | D4(b); **Stage B/C only — reported as `L2_confirmed_status = NOT_ASSESSED` in Stage A** (not assessed ≠ assessed and none qualified) |
| L3 external-validation executability | across L2-confirmed cohorts, how many independent validation units exist per compatible task stratum? | D3 units within one stratum ≥ the Stage C requirement (not computed in Stage A) |

Stage A selects no design, sets no N, tests nothing, and never chooses a stratum after outcomes exist: strata are
reported in full and any later choice is made in a frozen Stage C protocol before any value is opened.

**1.x Sensitivity limitation (owner, 2026-10-27).** Let I = 1 denote a study that truly contains ≥ 50 independent human
subjects in the chosen assay and E = 1 a study whose deposited metadata satisfy D1. The strict census counts E = 1. Since
P(E = 1 | I = 1) < 1 — the ISA source node is often a study-level placeholder (pilot 25b: 4 of 6 human studies) and the
template carries no person key — N_L1,strict ≤ N_truly-eligible holds only if the classifier makes no false-positive
independence assignment, which regex and documentary evidence cannot guarantee; the inequality is a methodological
goal. The UNRESOLVED band is the identifiability gap, not biological dependence. Any manuscript statement is of the form
"no cohort met the prespecified affirmative independence-evidence criteria among the evaluated accessions".

## 2. Definitions

**D1. Independent biological subject.** A human individual contributing exactly one biological sample to the
chosen assay. Two evidence routes, both affirmative. Evidence of *study design* and evidence of
*individual-subject independence* are recorded in separate columns and never substituted for each other.
- **R1 explicit identifier.** A column on the *source* node identifying a person — `Characteristics[Subject ID]`,
  `[Participant]`, `[Individual]`, `[Donor]`, `[Patient]`, `[Volunteer]` (pattern frozen in code) — present and
  non-sentinel on every biological row, exactly one `Sample Name` per identifier in the chosen assay.
- **R2 documented one-person-one-source design (A8, two layers; all conditions required).**
  *Normalisation* for every text test: case-fold; underscores and hyphens → spaces; collapse whitespace.
  *Layer 1 — design classification* on `Study Design Type` descriptors and `Study Description`:
  positive set `CROSS_SECTIONAL = {cross sectional, case control, case control study, matched case control,
  population based cross sectional, observational cross sectional}`; negative set `LONGITUDINAL = {longitudinal,
  time series, repeated measures, crossover, cross over, prospective cohort, follow up, before and after, pre post
  intervention, paired sampling, serial sampling}`. Layer 1 passes iff ≥ 1 positive term matches and 0 negative
  terms match anywhere in the permitted documentation. `cohort`, `observational`, `case series`, `intervention`,
  `clinical study` alone never classify positively. A negative term does not prove dependence; it blocks R2.
  *Layer 2 — affirmative one-person-one-source statement* in `Study Description` (sentence-split on `. ; :`):
  a sentence qualifies iff it matches one of the frozen templates **and** contains no negation token
  (`not, no, never, n't, without, except, rather than, instead of, unless`) **and** no plural-sampling token
  (`samples per, multiple samples, several samples, serial, repeated, longitudinal, each visit, each time point,
  subset of`). Templates (regex over normalised text; `P` = `participant|subject|patient|individual|donor|volunteer|person|
  woman|man|child|case|control`):
  T1 `(one|a single|1) (sample|specimen|biopsy|blood draw|collection) (was |were )?(collected|obtained|taken|drawn|
  provided) (from|for|by) (each|every) P`;
  T2 `(each|every) P (contributed|provided|donated|gave|supplied) (one|a single|1) (sample|specimen|biopsy)`;
  T3 `(samples|specimens) (were )?(collected|obtained|taken) from (\d+) (unique|distinct|different|independent|
  separate) P(s)?` **with the integer equal to the number of biological sources in the chosen assay**;
  T4 `(\d+) (unique|distinct|different|independent|separate) P(s)? (were |was )?(enrolled|recruited|included|
  studied|sampled)` **with the same integer equality** and one sample per source in the sheet.
  Unmatched wording is UNRESOLVED; nothing is interpreted after inspection. **Limitation (permanent):** a regex
  recogniser is tested against its synthetic suite; it cannot be proven to detect every contradictory sentence in
  arbitrary natural language. Its outputs are evidence grades under a frozen rule, not a reading of the study. ("patients were not sampled once"
  fails the negation guard; "a subset of participants gave one sample" fails the plural-sampling guard.)
  **Evidence grades (owner, A10).** R2 yields one of two grades: `R2_T12` when a T1/T2 sentence qualifies
  (explicit one-sample-per-participant statement consistent with the sheet) — counts toward the **strict** L1;
  `R2_COUNT_SUPPORTED` when only a T3/T4 (or verb-first T4b: "enrolled N distinct patients") count sentence
  qualifies with the integer **exactly** equal to the biological source count (A11: no ±1, no percentage, no
  manual exception) — equal cardinalities do not prove a person-to-source bijection, so this grade is reported in a
  separate **sensitivity** L1 count and never in the strict count. Both require Layer 1, no demonstrated repeated
  sampling anywhere in the description (a positive sentence never overrides a repeated-sampling sentence), one
  `Sample Name` per `Source Name`, and ≥ 50 sources. Guard tokens are matched at word boundaries (`no` does not
  match inside `normal`). Frequency words (`twice`, `weekly`, `per week`, …) are plural-sampling tokens.
  Any failure → UNRESOLVED with the layer recorded.
`SOURCE_UNIQUE` — distinct source labels with one sample each but neither R1 nor R2 — is **UNRESOLVED**, never
confirmed: unique labels are bookkeeping, not evidence of distinct individuals. Age/sex characteristics on the
source node are recorded as supporting descriptors and are not an evidence route. Sentinels → UNRESOLVED; never
imputed.

**D2. Eligible cohort.** One study, one chosen assay, L1 confirmed with ≥ 50 subjects, human, one matrix, no
demonstrated repeated sampling (§3.5 M6), not touched (§6).

**D3. Independent validation unit.** Eligible cohorts merged when: same submitter/institute *and* overlapping
sample or subject identifiers; same subjects under several accessions or assays; one cohort in both repositories;
declared related-study links. Units are reported per stratum. Accessions ≠ units ≠ compatible units.

**D4. Task-feasible cohort — two grades.** (a) *L2-preliminary*: an eligible MS or NMR cohort whose
**validated** repository-level count indicates ≥ 60 named features (P4V's hashed 20/20/20 blocks of 60; [A1]
RefMet naming for the moderator only). A repository count may count identifications, not usable columns; it is a
proxy. `FEATURE_COUNT_UNKNOWN` is neither 0 nor ≥ 60 and places the cohort in L2-preliminary-unresolved, giving
N_prelim,confirmed ≤ N_prelim,true ≤ N_prelim,possible. (b) *L2-confirmed*: the assay yields a usable
subject-by-feature matrix for the frozen task — ≥ 60 columns surviving Stage C's eligibility rules (availability,
duplicate-label and variance guards), a design/held-out subject split on L1-confirmed subjects, and prediction
targets that are themselves metabolite columns (the acquisition task needs no clinical target). (b) requires
value-adjacent or value access and belongs to Stage B/C; **in Stage A it is reported as `NOT_ASSESSED`**, never as a
count, so no Stage A figure is ever described as a demonstrated modelling capability.

## 3. Route M — MetaboLights (primary)

**3.1 Why primary.** Subject (`Source Name`) and sample (`Sample Name`) are structurally distinct in the sample
sheet; study documentation is in the investigation file; measured values are in separate metabolite-assignment
files. The value boundary is a file boundary — conditional on the structural validation of §3.3, since a filename
rule is not proof of content.

**3.2 Universe.** All public MTBLS accessions at the census date from one listing source fixed in the schema
lock (§7.1). No programme branch has used MetaboLights; cross-repository duplicates are handled by D3/§6.

**3.3 Access allowlist and per-class validators (fail-closed; refusal in code before any request).**
- Permitted requests: study listing; per-study file-name listing; **[A-1]** a names-only HTML index of the public FTP study
  folder (`{FTP}/{MTBLS}/`; `href` names recorded); **[A-2]** download by exact name of the `i_*.txt`, `s_*.txt`,
  `a_*.txt` names **that the mirror index itself lists** (REST names are compared and counted, never used as a fetch
  name; an index class other than `INDEX_ACCESSIBLE` leaves the fetch set undefined → halt), through the **single FETCH template** `{FTP}/{MTBLS}/{name}` (name
  percent-encoded; host must be `ftp.ebi.ac.uk`); every attempt logged; a non-200, a redirect, an oversize or a binary
  response halts. **[A13-a]** No second candidate exists; one is added only by a dated amendment supported by observed
  interface evidence (the FTPDIR classes of §3.3a). [A2] the `Metabolite Assignment File` *name* in the assay sheet may be read; the file is never
  fetched.
- Forbidden: `m_*.tsv` and any metabolite-assignment file; `FILES/`; raw, derived, mzML, spectral, image, results
  files; any bundled study download (the public CLI's default includes the MAF and is not used; **[A-1]** the REST
  `…/studies/{MTBLS}/download` endpoint family is therefore never a FETCH candidate); **any HTTP
  redirect** (never followed; a halt with the Location recorded as text); **[A-1]** any FETCH response whose
  `Content-Length` exceeds the size cap or whose `Content-Type` is zip/gzip/octet-stream/tar (halt before the body is
  read) or whose body starts with a zip/gzip signature (halt; never parsed).
- **Investigation validator (`i_`).** Vertical label–value structure: lines are `LABEL<TAB>value…`; accepted
  section headers are the ISA set (`ONTOLOGY SOURCE REFERENCE`, `INVESTIGATION`, `INVESTIGATION PUBLICATIONS`,
  `INVESTIGATION CONTACTS`, `STUDY`, `STUDY DESIGN DESCRIPTORS`, `STUDY PUBLICATIONS`, `STUDY FACTORS`, `STUDY
  ASSAYS`, `STUDY PROTOCOLS`, `STUDY CONTACTS`); any line outside a known section, or any section with > 50
  value columns, refuses the file. Fields read: `Study Title`, `Study Description`, `Study Design Type`, `Study
  Factor Name/Type`, `Study Assay Technology Type`, `Study Assay File Name`, `Study File Name`, `Study Protocol
  Name/Type`, `Study Person Last Name/Affiliation`, `Comment[…]` cross-references.
- **Sample-sheet validator (`s_`).** Tabular; header must consist of `Source Name`, `Sample Name`,
  `Characteristics[…]`, `Factor Value[…]`, `Protocol REF`, `Parameter Value[…]`, **[A-4]** `Material Type`,
  `Description`, `Provider`, `Performer`, `Date`, `Term Source REF`, `Term Accession Number`, `Unit`, `Comment[…]`
  only (the ISA-Tab 1.0 sample-table vocabulary); `Source Name` and `Sample Name` mandatory.
- **Assay-sheet validator (`a_`).** Tabular; header from the ISA-Tab 1.0 assay-table vocabulary **[A-4]**: `Sample Name`,
  `Protocol REF`, `Parameter Value[…]`, `Extract Name`, `Labeled Extract Name`, `Label`, `Assay Name`, `MS Assay Name`,
  `NMR Assay Name`, `Hybridization Assay Name`, `Scan Name`, `Normalization Name`, `Data Transformation Name`,
  `Array Design REF`, `Raw Spectral Data File`, `Derived Spectral Data File`, `Raw Data File`, `Derived Data File`,
  `Metabolite Assignment File`, `Array Data File`, `Array Data Matrix File`, `Derived Array Data File`, `Derived Array
  Data Matrix File`, `Image File`, `Spot Picking File`, `Protein Assignment File`, `Peptide Assignment File`, `Post
  Translational Modification Assignment File`, `Free Induction Decay Data File`, `Acquisition Parameter Data File`,
  `Characteristics[…]`, `Factor Value[…]`, `Material Type`, `Term Source REF`, `Term Accession Number`, `Unit`,
  `Comment[…]`, `Date`, `Performer` only. Any other header, including a near-miss spelling, refuses the file.
- **Numeric guard (secondary to structure, never a substitute):** in `s_`/`a_`, > 50 % numeric cells outside
  `Characteristics[Age…]`, `Unit`, `Parameter Value[…]`, `Factor Value[…]`, `Comment[…]`, `Term Accession Number`
  columns → refuse **[A26, approved 2026-10-26 after pilot 25; was 5 % without `Factor Value`]**; the fraction is reported for
  every accepted sheet. Every accepted file: name, size,
  SHA-256, accepted columns logged. Any refusal → study `FORMAT-REFUSED`, counted, not screened.

**3.3a FTPDIR diagnostic classes (fixed before run 2; [A13]).** Order per run: phase 1 — for each exemplar, FILES then
FTPDIR; phase 2 — metadata FETCH. Classes: `INDEX_ACCESSIBLE` — index 200 with ≥ 1 `href`; report lists listed
metadata names present in / absent from the mirror, mirror `i_/s_/a_` names not in the listing, and forbidden names seen
(counted, never fetched); interpretation: a present name whose FETCH fails → the file URL behaviour needs investigation;
an absent name → a REST/mirror filename discrepancy is supported. `DIRECTORY_404` — directory template unresolved; the
mirror's contents are **not** established and no absence is claimed. `DIRECTORY_403` — index inaccessible; file existence
unresolved. `DIRECTORY_OTHER_<code>` — unresolved. Redirect → halt, destination never requested. 200 without a
recognisable `href` structure → halt, no names inferred. No class triggers any alternative path.

**3.4 Analysis choice.** Among `a_*.txt` whose `Study Assay Technology Type` is mass spectrometry or NMR, the
one with most sample rows; tie → lowest file name. Stage C is bound to it.

**3.5 Qualification (evaluation order; first failure recorded).**
- M1 human: `Characteristics[Organism]` = Homo sapiens / NCBITaxon:9606 on every biological row; an explicit other
  organism → EXCLUDED. **[A30]** An organism cell naming an analytical reference (`reference compound`, `reference
  material`, `reference standard`, `blank`, `qc`, `quality control`, `standard mixture`, `calibrant`) is a QC row on the
  strength of that explicit designation, provided no other metadata conflicts with it, and is removed under M2 — a
  missing taxonomy accession is not corroboration, only the absence of a conflict; an explicit analytical label that
  carries a taxonomy accession (a conflict) is QC only if another designated column affirms the analytical role, else
  the study is UNRESOLVED; a missing organism with no QC evidence → study UNRESOLVED; the cell's text takes precedence
  over its accession; raw organism values and the removed/missing/conflict counts are recorded. **[A30 as clarified by
  the owner, 2026-10-27]**
- M2 QC/blank/pool removal by a frozen pattern on `Sample Name`, `Characteristics[Material type]`, `Factor
  Value[…]`, `Comment[…]`; any omitted row must match it; unidentified rows → UNRESOLVED.
- M3 identity: `Sample Name` present, non-sentinel, unique among biological rows; `Source Name` present,
  non-sentinel; else UNRESOLVED.
- M4 matrix: a single `Characteristics[Organism part]` (or `[Sample type]`) value among biological rows; missing
  on any row → UNRESOLVED; > 1 → EXCLUDED.
- M5 subject independence (D1): R1 → CONFIRMED; else R2 → CONFIRMED; else `SOURCE_UNIQUE` → UNRESOLVED; a
  source or identifier with > 1 sample, or a pooling/splitting `Protocol REF` on the path → EXCLUDED; **[A27, approved
  2026-10-26 after pilot 25]** except the degenerate case of one `Source Name` for every sample, which is `SOURCE_DEGENERATE` →
  UNRESOLVED (the source node is uninformative, not evidence of repeated sampling).
- M6 temporal structure, **two levels** [A4]; A9: `CONFIRMED-with-note` = D1 established by R1/R2 **and**
  temporal dependence not demonstrated; a temporal indicator never compensates for missing R1/R2 evidence. *Strong* indicators (`timepoint`, `time point`, `visit`,
  `follow-up`, `collection time`, `sampling day/week/month`, `baseline/endpoint` as factor names, `day N`/`week N`
  level patterns; matched in word-boundary, underscore and camel-case forms) with > 1 level → the prespecified
  repeated-measures assessment: EXCLUDED if any source has > 1 sample across levels; otherwise (independent groups
  sampled at different times) CONFIRMED-with-note if and only if M5 already confirmed; else UNRESOLVED. *Ambiguous* indicators (`stage`, `phase`, `pre`, `post`,
  `collection`, `cycle`, `duration`) → UNRESOLVED unless the design descriptor (R2 set) or a one-sample-per-source
  structure under R1 establishes non-temporal use; then CONFIRMED-with-note. A temporal factor never by itself
  establishes repeated sampling.
- M7 n ≥ 50 L1-confirmed subjects [A5]; else EXCLUDED.
- M8 feature count [A3]: a repository-level named-feature count is used only if the schema check (§7.1)
  demonstrates the field exists and defines what it counts; otherwise `FEATURE_COUNT_UNKNOWN`. No
  metabolite-assignment file is opened in Stage A. ≥ 60 → L2 confirmed (preliminary); < 60 → L2 excluded;
  UNKNOWN → L2-preliminary unresolved. `L2_confirmed_status = NOT_ASSESSED`.
- M9 touched (§6) → EXCLUDED.

**3.6 Outputs.** `stageA_M_cohorts.csv` (rule reached per level; separate columns `design_evidence` (Layer 1
result and matched terms), `subject_evidence` (R1 / R2 with matched template and sentence offset / SOURCE_UNIQUE),
L1 status; counts; assay; matrix; platform class; feature count or UNKNOWN; L2-preliminary status), `stageA_M_units.csv`, `stageA_M_exclusions.csv`, `stageA_M_summary.json` (L1/L2/L3
confirmed/unresolved/excluded by stratum), SHA256SUMS, run log with every URL.

**3.7 Documentary-recovery route for degenerate-source studies (A31-D; DEVELOPMENT ONLY; no retrieval authorised).**
Scope: studies with status UNRESOLVED and grade `SOURCE_DEGENERATE` (one `Source Name` for all samples) or
`SOURCE_UNIQUE`, human (M1 passed), one matrix, MS/NMR assay, n_samples ≥ 50. Principle (owner): when the deposited
metadata cannot identify independence, obtain independent evidence; never weaken D1. Evidence that can recover a study:
(a) an author-confirmed sample-to-participant mapping for the chosen assay; (b) an explicitly documented
one-person-one-sample design with a verified participant count equal to the number of included biological samples and
assay membership established; (c) another independently auditable source (a linked publication's methods stating the
design and count, matched by exact integer to the assay). Procedure, frozen before any evidence is examined: the
question set and acceptance criteria are written and hashed first; candidates are selected by the census outputs
alone; each piece of evidence is recorded with source, date, verbatim text and the rule it satisfies; evaluation is
blind to any value. Outcome grades: `DOC_RECOVERED_STRICT` (a or b with exact count), `DOC_RECOVERED_COUNT` (c with
exact count only — sensitivity), `DOC_UNRESOLVED`, `DOC_CONTRADICTED` (evidence contradicts the claimed independent-sampling structure of the chosen assay; whether a
separately defined subset remains usable requires its own evidence). Reported
as `L1_documentary_recovered` columns, **never added to** `L1_strict` or `L1_sensitivity`. Author correspondence is a
legitimate source only with the questions and criteria frozen before the first message is sent and every reply filed
verbatim. Nothing in §3.7 is executable until its own protocol section, test suite and authorisation exist.

## 4. Route W — Workbench mwTab (secondary; [A6] tests only, no run before Route M is assessed)

**4.1** The `#SUBJECT_SAMPLE_FACTORS` block (subject ID, sample ID, factors, additional data) carries a `SUBJECT`
field that REST `factors` omits; the specification allows it to be absent; prevalence unknown. Values share the
file, so value-blindness is a parser property.

**4.2 Boundary reader.** Streams line by line; parses only the metadata sections (`#METABOLOMICS WORKBENCH` …
`#MS`/`#NMR`); **terminates and closes the connection** at the first line matching `#MS_METABOLITE_DATA`,
`#NMR_METABOLITE_DATA`, `#NMR_BINNED_DATA`, `*_DATA_START`, `#METABOLITES`, or a line whose fields from the third
onward are ≥ 90% numeric; records termination byte offset, lines read, SHA-256 of bytes read; writes nothing to
disk; a file whose factor block is incomplete at termination is `FORMAT-REFUSED`.

**4.3 Boundary tests (authorised; touched fixtures only).** `ST000568_AN000873/874/875`, `ST002733_AN004431`,
`ST002955_AN004852`, `ST003368_AN005520`, `ST004571_AN007678` mwTab files in `paper4_external_validation/`.
Assertions per fixture: termination offset < offset of the first `*_DATA_START` line found by an independent byte
search; no numeric matrix row in parsed output; parsed factor-row count = cached REST `factors` row count
(`b0_cache/http`); `SUBJECT` prevalence reported (a biased first look, labelled so). Synthetic files: marker at
line 1, mid-block, absent → refuse. Tests hashed before any run is proposed.

**4.4 Qualification.** W1–W9 mirror M1–M9 with `SUBJECT` in R1's role and the mwTab `#STUDY` summary text in R2's
role under the same frozen recognition rule; species, counts, analysis types and feature counts come from the B0
cache without new requests. Candidates for a future run: the 104 + 20 Workbench studies — not fetched here.

## 5. Units, strata and compatibility (L3)

After qualification, D3 merging uses normalised institute/submitter, title similarity, exact overlap of sample or
subject identifiers, and declared cross-references. Report units by (repository × modality × platform class ×
matrix). L3 is reported as counts per stratum with no stratum selected; Stage C's frozen protocol chooses the
task stratum before any value is opened and never after.

## 6. Touched set

The 217-ID Workbench union (`touched_union_strict_full.json` `8a42d061…`), the §4.3 fixtures, the two schema
exemplars (§7.1), and any MetaboLights study found by D3 to be a deposit of a touched Workbench cohort. Computed
and hashed before any route runs.

## 7. Locks and order (none of 7.2–7.5 authorised)

**7.1 Schema lock (first).** Order, with **no byte of any locked file changing after step 3**:
1. The owner writes `stageA_live_authorisation.txt` containing exactly one line of the form
   `LIVE SCHEMA INSPECTION AUTHORISED <YYYY-MM-DD> exemplars=2 protocol=EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0.6.7.md owner=<name>`
   and **nothing else** (surrounding whitespace only; a second line of any kind invalidates it). This file *is* the
   dated amendment; the protocol text is not edited for it. Any last protocol edit happens before step 2.
2. `python3 -I stagea_isa.py selftest && python3 -I stageA_schema_probe.py selftest` on exactly the bytes to be locked.
3. `python3 -I stageA_schema_probe.py plan --protocol <this file> --auth stageA_live_authorisation.txt` →
   `stageA_schema_lock.json` with four hashes (probe, isa, protocol, authorisation file), the endpoint templates
   (UNVERIFIED), exemplar rule, fields, outputs, failure handling, authorisation conditions, `owner_signoff: null`.
4. The owner adds `owner_signoff` to the lock (the lock is the only file edited after step 3, and it is not hashed into
   itself).
5. `run`: refused unless the four hashes equal the bytes present, the sign-off is present, and the authorisation line
   is valid and names this protocol file; otherwise runs the tested workflow and exits nonzero on any halt.
The lock binds the executable bytes that passed step 2, never an earlier development version. Before any request: the listing endpoint; the two exemplar accessions [A7] = the two
lowest-numbered public MTBLS accessions in the listing, chosen mechanically and thereafter excluded from the
census; the exact fields to be inspected (validator acceptance of `i_`/`s_`/`a_`; presence and definition of a
study-level metabolite count; `Study Design Type` vocabulary; `Characteristics` column names on the source node);
expected outputs (`stageA_schema_report.json`); failure handling (any refusal halts). Hashed with a dated
sign-off → **7.2 schema inspection** (exemplars only) → **7.3 test suites** (validators on synthetic sheets:
pooled sources, split samples, strong/ambiguous time factors, missing organism part, numeric injection; mwTab
boundary tests §4.3) → **7.4 Stage A freeze** (protocol + code hashed, dated sign-off, A1–A7 recorded) → **7.5
census** (Route M; Route W only on separate authorisation). **Precondition for 7.5 (owner, 2026-10-25):** listing
completeness must be ESTABLISHED, or a pagination procedure must be specified and tested so that the universe is
enumerated in full; a census over an UNRESOLVED listing is not authorised, because A7 and every count would be
relative to an unknown subset of the repository. No rule change after 7.4 without a dated amendment;
any change after an output is seen is a recorded deviation.

**7.1a Run 1 record (v0.6.2 lock; owner-signed 2026-10-25; executed by the owner from github.pathway).** Lock hashes:
probe v0.3.2 `4b88eba7…c060`, isa v0.2.0 `3c2feabd…b1be`, protocol v0.6.2 `65ebae7e…791c`, authorisation
`d3071cb4…4677` (line `LIVE SCHEMA INSPECTION AUTHORISED 2026-10-25 exemplars=2 protocol=…v0.6.2.md owner=Anas`).
Requests, all within the v0.6.2 allowlist, in order: (1) LIST `https://www.ebi.ac.uk/metabolights/ws/studies` → 200,
44 057 B, 3584 accessions, no total field → completeness UNRESOLVED; exemplars MTBLS1, MTBLS2 (A7, relative to the
returned set). (2) FILES `…/ws/studies/MTBLS1/files?include_raw_data=false` → 200, 779 B, shape recognised; one
forbidden name counted, never fetched. (3) FETCH `https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS1/a_MTBLS1_NMR_metabolite_profiling_NMR_spectroscopy.txt`
→ **404**; halt; report and fetch log written; exit 1. Outcome: endpoint_reachable = yes (LIST, FILES);
response_shape_recognised = yes (LIST, FILES); metadata_validated = **none** (no metadata file retrieved);
`a3_metabolite_count_field` absent in LIST/FILES → `FEATURE_COUNT_UNKNOWN` stands; MTBLS2 not reached. No URL was edited
and no request was retried; the correction is this dated amendment (A-1), which requires a new four-file lock and a
new owner sign-off before run 2. What run 1 does **not** show: whether the FTP mirror lacks the *path* or the *file
name* (the FILES listing is the working-directory name list; the mirror may name files differently) — the FTPDIR
index of A-1 is the names-only instrument for that question.

**7.1b A13 record (owner, 2026-10-25).** Owner verified the v0.6.3 bundle hashes (bundle `9f4fcf8f…990d`, protocol
`3fbec18a…ad44`, probe v0.3.3 `9aa8f022…c95c`, isa `3c2feabd…b1be`; both selftests passed) and **did not approve it for
live execution**. A13-a approved: diagnostic-only run 2 (LIST, FILES, FTPDIR, original FTP FETCH template; redirect
refusal, binary/bundle guard, size cap and two-exemplar limit retained; census and any MAF/raw access not authorised).
A13-b rejected: the unverified web FETCH route. Requirement: an unavailable directory index (403/404) must be
distinguished from a verified absence of the requested name, and the report must preserve that distinction without
guessing alternative paths. Staging: run 2 diagnoses the public file interface; run 3, only if necessary, tests a
retrieval mechanism supported by observed interface evidence; then Stage A cohort qualification proper. Owner's
constraint: this is not an invitation to indefinite engineering.

**7.1c Run 2 record (v0.6.4 lock; owner-signed `Anas 2026-10-25T20:37:01Z`).** Six requests, all within the allowlist, in
the two-phase order; LIST 200 (3584; completeness UNRESOLVED); FILES and FTPDIR 200 for both exemplars
(`INDEX_ACCESSIBLE`); FETCH of the REST-named MTBLS1 assay sheet 404 → halt. Diagnosis and two disclosed reporting
artefacts (transient `_meta` key in the halted report; exit status not captured under zsh — use `${pipestatus[1]}`) in
`STAGEA_RUN2_RESULT_v1.0.md`. No metadata file and no value retrieved.

**7.1d Run 3 attempt 1 record (v0.6.5 lock; owner-signed `Anas 2026-10-25T20:50:02Z`; lock `created_utc`
2026-10-25T20:48:03Z; authorisation `a657c977…dc81`).** First request LIST → `RemoteDisconnected` → uncaught → traceback;
`EXIT STATUS: 1` captured with `${pipestatus[1]}`; no report, no fetch log, no bytes received, nothing fetched. Not a
result about the repository; a defect of the probe's error handling, corrected prospectively by A-3.
**Census precondition added (owner, A14 qualification):** the fetch set may contain several files of one ISA class; the
two-exemplar inspection records each; the §3.4 one-assay choice is a census rule and the census must not assume one
`i_`/`s_`/`a_` per study. Provenance per fetched file: URL, status, size, SHA-256 (fetch log) and run timestamp (record).

**7.1e Run 3 attempt 2 record (v0.6.6 lock; owner-signed `Anas 2026-10-25T20:56:02Z`).** Six requests within the
allowlist: LIST 200; FILES + FTPDIR 200 for both exemplars (`INDEX_ACCESSIBLE`; MTBLS2 also lists
`MTBLS2.announcement.json`, `MTBLS2.mhd.json` — named, never fetched); FETCH of the mirror's
`a_MTBLS1_metabolite_profiling_NMR_spectroscopy.txt` → 200, 120 888 B, `cd615abd…de7a` → FORMAT-REFUSED on two ISA-Tab
process-node columns → halt, exit 1 (captured). A-2 verified live; no `i_`/`s_` retrieved; no value opened. Full record:
`STAGEA_RUN3_RESULT_v1.0.md`.

**7.1f Run 3 attempt 3 and close-out (owner, 2026-10-25).** Under the v0.6.7 lock (isa v0.2.1 `161036c5…9f48`, probe
v0.3.7 `bd676561…fea6`, protocol v0.6.7 `ca15b8e0…c7d`) the inspection completed: `halted: null`, exit 0 (captured),
`n_accessions_listed` 3584, 2 forbidden names counted and never fetched, `count_keys_seen` empty (A3 →
`FEATURE_COUNT_UNKNOWN`). Both exemplars: investigation, sample and assay sheets retrieved from the mirror under their own
names and validated. Findings as reported by the owner: MTBLS2 is a plant study (*Arabidopsis thaliana*; 16 sample rows,
16 assay rows; genotype and replicate factors) → not a human cohort; MTBLS1: human status and a person identifier not
demonstrated by the inspected fields; L2-confirmed NOT_ASSESSED for both; both excluded from every census count (A7).
These findings validate the inspection mechanism, not the availability of eligible cohorts. The owner closed the
schema-inspection phase: no further two-exemplar inspection; the run-3 executable and conclusions are not to be modified.
Formal record: `STAGEA_RUN3_RESULT_v1.1.md` (from `stageA_out3c/`: 11 requests, per-file SHA-256; R1 person key
absent on both exemplars; identical 86-label investigation template and identical source-node characteristics on both —
an observation that sets expectations for the strict/sensitivity L1 split, not a rule change).

**7.5a A17 — Public-universe verification (precondition of 7.5; separately locked; metadata-only).**
*Question (one):* does `{REST}/studies` enumerate the complete public accession universe, or a partial / paginated set?
*Requests (the only ones):* LIST `{REST}/studies` × 3 (≥ 2 s apart) with response headers recorded; FTPROOT `{FTP}/`
× 2 — the names-only HTML index of the public mirror root, whose `MTBLS\d+/` folders are the mirror's own public set (the
source of record since A14). No study folder, no study file, no value; redirects never followed; a non-200, oversize,
malformed or transport failure is a recorded halt.
*Verdict rules (fixed here, before execution):*
- `ESTABLISHED_BY_CONCORDANCE` — every LIST repeat byte-identical, every FTPROOT repeat identical, REST set = mirror set,
  no pagination-/total-like JSON key (`total*`, `page*`, `size`, `limit`, `offset`, `next`, `last`, `count`, …) and no
  pagination-like header (`Link`, `Content-Range`, `X-Total-Count`, `X-Pagination-*`), no duplicate, no non-conforming token.
- `BOUNDED` — repeats stable, sets differ: census universe **U = REST ∩ mirror** (both interfaces call the accession
  public); band [|∩|, |∪|]; every accession in the symmetric difference is listed by name in the manifest with status
  `UNRESOLVED`, counted, and never screened unless a later dated amendment resolves it (the B0 K_lo/K_hi convention).
- `UNRESOLVED` — anything else (unstable repeats, a pagination-like key or header, duplicates, non-conforming tokens,
  FTPROOT unavailable, any halt). The census remains blocked.
*Outputs:* `stageA_universe_manifest.csv` (accession, in_rest, in_mirror, status; sorted numerically),
`stageA_universe_report.json` (repeat hashes, headers, JSON keys, comparison, verdict, reasons, manifest SHA-256, run
timestamp), `stageA_universe_fetch_log.txt`, `SHA256SUMS.txt`. Exit 0 only on `ESTABLISHED_BY_CONCORDANCE` or `BOUNDED`.
*Lock:* three hashes (`stageA_universe_audit.py`, this protocol, `stageA_universe_authorisation.txt` holding exactly
`LIVE UNIVERSE AUDIT AUTHORISED <YYYY-MM-DD> protocol=EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0.6.9.md owner=<name>`) and an
owner sign-off; `run` refuses otherwise. **[A-5]** Root-index parsing: relative `MTBLS\d+/` hrefs are folders; absolute
or parent hrefs are navigation (listed, ignored); any other relative folder href is non-conforming. FTPROOT per-read
timeout 600 s; a halt after one LIST and one root success writes a `PROVISIONAL_` comparison and manifest, which no
census may read. *Run 1 (v0.6.8 lock, 2026-10-25):* UNRESOLVED by transport (second root read timed out); provisional
evidence in `STAGEA_A17_RUN1_RESULT_v1.0.md`. *What A17 does not do:* it does not read MetaboLights API documentation, probe
pagination parameters, or request any endpoint other than the two above; if the verdict is `UNRESOLVED` because a
pagination signal exists, the pagination contract is specified by a later amendment from documentation, not discovered
by trial. *Interpretation limit:* concordance of two interfaces operated by the same repository is evidence that both
enumerate the same public set, not proof that the repository has no further public studies; the protocol accepts that
limit explicitly and the census universe is defined as "public according to both interfaces on the audit date".

**7.5b A17 run 2 result (v0.6.9 lock; owner sign-off 2026-10-25T21:40:24Z; formal record
`STAGEA_A17_RUN2_RESULT_v1.0.md` `8cb7e258…eedf`, written from `stageA_universe_out2/` and verified).** LIST ×3 stable; mirror root ×2 stable;
no transport failure; `halted: null`; verdict **BOUNDED**; REST 3584, mirror 3596, both 3584, REST-only 0, mirror-only 12
(MTBLS9, MTBLS30, MTBLS202, MTBLS224, MTBLS351, MTBLS539, MTBLS557, MTBLS643, MTBLS728, MTBLS926, MTBLS4938, MTBLS14881),
union 3596; exit 0. Census universe U = 3584 (intersection), band [3584, 3596]; the twelve are UNRESOLVED and are not
screened under this protocol. Interpretation limit of §7.5a stands: U is "public according to both interfaces on
2026-10-25", not a proof of repository completeness. Screening pool = U minus the two A7 exemplars = **3582 accessions**
(an accession count, not a cohort count). The manifest's SHA-256 is the value the census authorisation must carry.

**7.4a A18 — Stage A census executable (`stageA_census.py` v0.1.0; offline only until 7.4 freeze and 7.5 authorisation).**
*Input:* `stageA_universe_manifest.csv` from A17, SHA-256 verified against the authorisation line and the lock; rows with
status `UNIVERSE` only; a `PROVISIONAL_` manifest is refused; exemplars excluded (A7). *Requests per study (mirror only;
REST never used):* the folder index (names only) and every `i_/s_/a_` it lists, under the §3.3 guards. *Rules:* §2 D1–D4,
§3.3 validators, §3.4 assay choice, §3.5 M1–M9 in order with the first failure recorded, §5 D3 merging, §6 touched set;
the D1 recogniser is `stagea_isa` v0.2.1 unchanged. *Implementation choices that the protocol text left open — each is a
decision the owner takes before the lock (A20–A23), not a rule discovered from data:*
- **A20 failure policy.** A transport failure, redirect, binary response, index failure or validator refusal on one
  study marks that study (`TRANSPORT_FAILED`, `FORMAT_REFUSED`, `INDEX_404/403/REDIRECT/UNRECOGNISED/OTHER`) and the
  census continues; nothing is retried within a run. A single, separately invoked and disclosed `--retry-pass` on the
  same lock re-requests only `TRANSPORT_FAILED` studies. Progress is checkpointed per study (JSONL); a restart never
  re-requests a recorded study. Every status is counted in the summary; none is silently dropped.
- **A21 pacing and duration.** Fixed 0.5 s pause before every request; per-read timeout 300 s. At ~4 requests per
  study and the mirror's observed latency, [Guessing] 3582 studies take several hours to a day; the run is resumable.
- **A22 M2 interpretation.** Biological rows = sample-sheet rows whose `Sample Name` is in the chosen assay and that do
  not match the frozen QC pattern (`qc|quality control|blank|pool(ed)|standard|std|reference|calibration|calibrant|
  solvent|wash|spike(d)`, word-bounded, on `Sample Name`, `Characteristics[Material/Sample type]`, every `Factor Value[…]`
  and `Comment[…]`). A sample-sheet row that is in no assay and not QC-matched → UNRESOLVED (M2 "unidentified row"); a
  row in another assay only is neither counted nor fatal; a chosen-assay `Sample Name` absent from the sample sheet →
  UNRESOLVED. M1 accepts `Homo sapiens`/`human` (normalised) or an accession containing `9606`, on every biological row.
  M4 uses `Characteristics[Organism part]`, else `[Sample type]`. M5 pooling: a `Protocol REF` value on a biological row
  matching `pool|pooled|pooling|split|splitting|aliquot` → EXCLUDED. M6 as §3.5 (a strong temporal factor with > 1
  level and one sample per source → CONFIRMED-with-note; an ambiguous factor without R1 or Layer 1 → UNRESOLVED).
- **A23 structural preconditions.** Exactly one `i_` and exactly one `s_` per study (else `FORMAT_REFUSED`, counted);
  any number of `a_`; §3.4 chooses among MS/NMR assays by row count, tie → lowest name; an `a_` on the mirror absent
  from the investigation's assay list is recorded (`assay_files_on_mirror_not_in_investigation`) and has no technology,
  so it is never chosen. D3 merging (§5) uses exact overlap of sample/source identifiers (not yet populated in v0.1.0 —
  identifiers are not retained across studies in this version; to be decided: retain hashed identifiers for D3) or same
  normalised first affiliation with title-token Jaccard ≥ 0.6; units reported per (repository × modality × platform
  class × matrix).
*Outputs:* `stageA_M_cohorts.csv` (one row per accession, all evidence columns, `design_evidence` and
`subject_evidence` separate), `stageA_M_exclusions.csv`, `stageA_M_units.csv`, `stageA_M_summary.json` (strict L1,
sensitivity L1, UNRESOLVED, status and rule counts, L2-preliminary all `FEATURE_COUNT_UNKNOWN`, `L2_confirmed_status =
NOT_ASSESSED`, L3 units by stratum, transport failures, request count), `stageA_M_files.csv` (every fetched file with
SHA-256), fetch log, checkpoint, `SHA256SUMS.txt`. *Lock:* five hashes (census, `stagea_isa`, protocol, authorisation
`LIVE STAGE A CENSUS AUTHORISED <date> universe=<manifest sha256> protocol=<this file> owner=<name>`, manifest) + sign-off;
`run` refuses if the authorisation's universe hash differs from the manifest present. *Offline tests (27):* R1, R2_T12,
R2_COUNT_SUPPORTED (sensitivity only), SOURCE_UNIQUE, non-human, two samples per source, two matrices, n = 49,
longitudinal descriptor, strong temporal factor, pooling protocol, no MS/NMR assay, index 404, transport failure, bad
sheet, redirect to a MAF (never followed), MTBLS1-layout fixture (→ UNRESOLVED SOURCE_UNIQUE, n = 132), summary counts,
D3 merge, resume, retry pass, manifest hash and provisional refusals, guards, authorisation line. **Expectation from
run 3 (not a rule):** MetaboLights' fixed template carries no person key; strict L1 will come from R2_T12 sentences,
and the UNRESOLVED band will be large. The census reports the band; it does not narrow it.

**7.4b A18 decisions (owner, 2026-10-25) and their implementation in census v0.1.1.**
- A20 APPROVED with the global-halt restriction → implemented: per-study classes `TRANSPORT_FAILED` (operational),
  `FORMAT_REFUSED` / `UNRESOLVED_STRUCTURE` (structural), `ACCESS_REFUSED` (server-side redirect/binary/oversize; never
  followed or parsed); a `Refused` from the URL guard is a protocol violation → `halted_protocol_violation` in the summary,
  outputs written, census stops. Tested.
- A21 APPROVED for a bounded pilot only → implemented as `pilot=N` in the authorisation line (first N accessions of the
  manifest in numeric order; a numeric-prefix pilot is operationally informative, not representative). The pilot's
  pacing/throughput findings are reported, not extrapolated. No pilot is authorised yet; N and its output directory are
  fixed by the owner before the lock.
- A22 CONDITIONALLY APPROVED → implemented as the two-tier QC rule above; tested on `reference group` in a Factor Value
  (→ UNRESOLVED) and `pool_1` / `standard_2` as Sample Names (→ removed). M1/M4/M5 unchanged; missing information
  yields UNRESOLVED, never eligibility.
- A23 APPROVED for deterministic handling → implemented: two sample sheets → `UNRESOLVED_STRUCTURE` (counted in the
  summary as `structure_unresolved`); assay reconciliation as in the header; both names preserved.
- D3 APPROVED with safeguards → implemented as above; `stageA_M_overlap_flags.csv` distinguishes flags from merges;
  `D3_suspected_overlap_flags` counted in the summary; no raw identifier in any artefact (tested).
- Universe: the census verifies the local manifest bytes against the SHA-256 in the authorisation line and the lock
  before any request; the result record is not a substitute (implemented since v0.1.0; unchanged).
**Open before the 7.4 freeze:** owner audit of the v0.1.1 bytes; pilot N; then the freeze record (protocol + `stagea_isa`
+ census hashed, dated sign-off, A1–A23 recorded) → pilot authorisation → pilot record → full-census authorisation.

**7.4c A24 pilot (owner, 2026-10-26): N = 25.** Selection: the first 25 non-exemplar `UNIVERSE` accessions of the A17
manifest in numeric order. Authorisation line: `LIVE STAGE A CENSUS AUTHORISED <date>
universe=3e483c7a2f5a4dd9968d404f78b45ec4d5583c15d75915550b315ee5720c9dc5 protocol=EXTERNAL_DATA_FEASIBILITY_PROTOCOL_v0.7.2.md
owner=<name> pilot=25`; fresh output directory; `run --pilot 25`. Evaluation criteria (all reported, none a pass/fail on
eligibility): exactly 25 attempted; manifest hash verified before any request; zero local access violations (a
violation → nonzero halt); zero forbidden retrievals; 25 terminal checkpoint records; SHA-256 for every retrieved ISA
file; resume verified (a second invocation requests nothing); every study classified; D3 confirmed vs suspected
distinguished; no measured value. Operational measurements: retrieval success rate, request durations (fetch log),
files per study, validator refusal rate, status distribution. A pilot that leaves studies UNRESOLVED, malformed or
server-failed has still succeeded operationally; it says nothing about L1 prevalence in the repository.

**7.4d Pilot 25 (v0.7.2 lock; owner sign-off `Anas 2026-10-26T17:08:09Z`; 17:08:15Z–17:11:52Z; exit 0; resume 0
requests).** 25 attempted; 112 requests; 87 files hashed; 0 access violations; 0 transport failures; statuses:
FORMAT_REFUSED 17 (numeric guard), EXCLUDED 5 (M5 ×4 single/multiple source, M1 ×1), UNRESOLVED 3 (M2); ELIGIBLE 0.
Throughput ≈ 1.9 s/request → [Likely] ≈ 9 h for 3582. Defect P1 (fetch log overwritten by the resume check) and
proposals A26/A27 in `STAGEA_PILOT25_RESULT_v1.0.md`. The v0.7.2 freeze record is withdrawn unsigned.

**7.4e Pilot 25b (v0.7.4 lock; owner sign-off `Anas 2026-10-27T02:38:59Z`).** 25 attempted; 111 requests journaled with
timing (median 1.0 s, p90 1.9 s); 0 access violations; 1 transport failure (MTBLS8, TLS handshake timeout 300 s);
FORMAT_REFUSED 0 (was 17); EXCLUDED 17 (M1 ×16, M5 ×1); UNRESOLVED 7 (`SOURCE_DEGENERATE` ×4, M2 ×3); ELIGIBLE 0.
MTBLS27 was EXCLUDED at M1 for blank / reference-compound organism rows → A30. Owner decisions 2026-10-27: A30
conditionally approved and implemented as specified; D1 relaxation rejected; §3.7 approved for development; MTBLS8
retry approved under the same lock (copy outputs first). Record `STAGEA_PILOT25B_RESULT_v1.0.md`.

**7.5c Census run 1 (v0.7.6 lock; owner sign-off `Anas 2026-10-27T03:07:32Z`; 03:07:32Z–06:16:40Z; exit 0).** 1810
studies screened (MTBLS3 … MTBLS10446): EXCLUDED 1253 (M1 1163, M4 42, M5 47, 3.4 1), UNRESOLVED 521 (`SOURCE_UNIQUE`
320, `SOURCE_DEGENERATE` 92, M2 97, M3 11, M1 1), FORMAT_REFUSED 24, ACCESS_REFUSED 9 (defect D-full-2),
UNRESOLVED_STRUCTURE 2, INDEX_UNRECOGNISED 1; ELIGIBLE 0. 459 reached M5; R1 never reached; R2 Layer 1 passed 8, Layer 2
qualified 0. From 06:01 the Mac's DNS failed: 1772 studies recorded `TRANSPORT_FAILED` without contacting the repository
(defect D-full-1). Interim record `STAGEA_CENSUS_RUN1_INTERIM_v1.0.md`.

**7.5d Census retry pass (v0.7.7 lock; owner sign-off `Anas 2026-10-27T17:26:59Z`; 17:26:59Z → 2026-10-28T03:51:13Z; exit 0).**
1781 studies re-requested (1772 outage casualties + 9 double-encoding); 8642 requests (8617 × 200, 23 ERR, 2 × 403).
Census totals over 3582: EXCLUDED 2436 (M1 2308, M4 67, M5 60, 3.4 1); UNRESOLVED 1052 (`SOURCE_UNIQUE` 612,
`SOURCE_DEGENERATE` 121, M2 302 of which 230 "no biological rows after QC removal", M3 13, M1 3, M4 1); FORMAT_REFUSED 58;
UNRESOLVED_STRUCTURE 9; INDEX_UNRECOGNISED 2; TRANSPORT_FAILED 23; FETCH_HTTP_403 2; ELIGIBLE 0. 793 reached M5; R1 0;
R2 Layer 1 passed 18, Layer 2 qualified 0. §3.7 candidates (n ≥ 50): 369. Caveats C1–C4 in
`STAGEA_CENSUS_RESULT_v0.9_PROVISIONAL.md`.

**7.5e A34 — diagnostic rescreen of the 230 "all rows QC" studies (proposed).** Same rules, census v0.1.8, new lock
(authorisation line naming this protocol, no `pilot=`), same output directory, `run --rescreen-allqc` (optionally together
with `--retry-pass` for the 25 operational failures if those have not been retried under v0.7.7). Outputs per study: status
(expected unchanged), `qc_triggers`, `organism_values_raw`, `person_key_columns`. Decision rule fixed now: if the triggers
show a column/token that marks biological samples in human studies, a rule amendment (A35) is drafted and the affected
studies are re-screened under it; if not, the 230 stay UNRESOLVED with the trigger recorded.

## 8. Transition to Stage B/C (separately authorised)

**Resolving `FEATURE_COUNT_UNKNOWN` prospectively.** Stage B may define a frozen, tested
*metabolite-assignment header reader* that reads the MAF's identification columns (`metabolite_identification`,
`database_identifier`, `chemical_formula`) and **no abundance column** — stopping at the column boundary, with the
abundance columns detected by header name and never parsed — producing a named-feature count and the RefMet/ChEBI
alignment for the moderator. This is a value-adjacent access and requires its own authorisation, lock and test
suite; it is not Stage A. Stage C opens values only under P4V v0.12 §6's split utility and ledger, with N, the
evaluation set, donor matching, stratum and decision rule frozen before the first design half. The internal
93-study archive stays unopened for P4V purposes until a Stage C design is frozen.

## 9. Decisions recorded (owner, 2026-10-25) and remaining items

A1 RefMet for moderator only · A2 MAF filename readable, file forbidden · A3 repository count if validated, else
UNKNOWN; no MAF in Stage A · A4 two-level temporal detector · A5 ≥ 50 independent subjects · A6 Route W tests
only · A7 two lowest-numbered public accessions, excluded thereafter · A8 two-layer R2 with the frozen descriptor
sets and templates of D1 · A9 `CONFIRMED-with-note` counts toward L1 only when D1 is established.
A10 approved: templates T1–T4 (+T4b) and guard lists are fixed verbatim in `stagea_isa.py` v0.2.0 with 29 offline selftests
(positive T1–T4, negation, subset, frequency, integer mismatch, verb-first count, contradictory sentences,
word-boundary `no`/`normal`, Layer 1 negatives, `cohort` alone, n < 50, duplicate source, sentinels, validators,
numeric guard, temporal camel/underscore). A11 approved: exact integer equality, no tolerance.
A12 given (owner sign-off 2026-10-25T20:06:48Z on the v0.6.2 lock; run 1 executed, §7.1a).
**A13 decided (owner, 2026-10-25): A13-a APPROVED, A13-b REJECTED** (§7.1b); run 2 executed (§7.1c).
**A14 decided (owner, 2026-10-25): APPROVED** (mirror index names; public release as source of record); run 3 attempt 1
aborted (§7.1d). **A15 decided (owner, 2026-10-25): A-3 accepted; run 3 attempt 2 executed** (§7.1e). **A16 decided (owner): A-4 accepted; run 3 attempt 3 completed, exit 0; schema inspection CLOSED** (§7.1f).
**A17 (owner, 2026-10-25): development APPROVED; run 1 executed under the v0.6.8 lock → UNRESOLVED by transport.**
**A19 decided (owner): A-5 accepted; A17 run 2 → BOUNDED (§7.5b); A17 COMPLETE.** **A18 (owner): offline development APPROVED; A20–A23 and D3 decided (§7.4b); live census NOT authorised.**
**A24 decided (owner, 2026-10-26): v0.1.1 audited; C1–C3 required and implemented in v0.1.2; N = 25.**
**A25 decided (owner): pilot 25 authorised and executed under the v0.7.2 lock** (§7.4d). **A26 / A27 decided (owner, 2026-10-26): APPROVED as prospective pilot amendments** (A26 an engineering threshold to be
observed, not a validated boundary; `SOURCE_DEGENERATE` outside both L1 counts). **A28 decided: same-25 repeat pilot
APPROVED in design; live execution on HOLD until C5–C7 and request-level provenance (this version) are audited.**
**A29 decided: v0.7.4 frozen (record v1.1 signed for the run); A28 executed → pilot 25b (§7.4e).** **A30 decided
(conditional) and implemented; D1 relaxation REJECTED; A31-M (census with A30) and A31-D (§3.7) APPROVED FOR
DEVELOPMENT.** **A32 decided (owner, 2026-10-27): CONDITIONALLY APPROVED — no further pilot; full census after the A30
clarification (this version), a signed freeze record v1.3 and a separate authorisation line without `pilot=`.** The
census must report separately: strict L1; sensitivity L1; `SOURCE_DEGENERATE`; `SOURCE_UNIQUE`; other UNRESOLVED
(structural / metadata); EXCLUDED by rule; operational failures; D3 metadata units; L2/L3 NOT_ASSESSED. **A32 executed: census run 1 (§7.5c).** **A33 decided and executed: retry pass under the v0.7.7 lock (§7.5d).** **A34 (open, owner):** (i) the final retry
invocation for the 25 operational failures under the existing v0.7.7 lock; (ii) accept census v0.1.8 and the diagnostic
rescreen §7.5e; sign freeze record v1.5. **§3.7 (open, owner):** whether to specify and authorise documentary recovery for
the 369 candidates. `run`. **A18 (owner): census executable development APPROVED offline** (synthetic studies and the two touched exemplar
fixtures only; no untouched cohort metadata during development); **execution NOT authorised** until A17 yields
`ESTABLISHED_BY_CONCORDANCE` or `BOUNDED` and the Stage A freeze (7.4) is signed. Interpretation constraint (owner): a 200 response establishes
reachability only; completeness of the listing is reported separately and, if UNRESOLVED, the A7 exemplars are the
minimum of what was returned, not of the repository.

**Next implementation priority (owner):** the separately locked schema-inspection utility (§7.1–7.2), not the
census.

Programme status unchanged: Paper 2 CLOSED; Track B P4V EXECUTION NO-GO on the Workbench REST route; v0.12 on
hold; six frozen programmes untouched; DEP-B working. Run 1 opened no metadata file and no value; this amendment
creates no access until the owner locks and signs it.
