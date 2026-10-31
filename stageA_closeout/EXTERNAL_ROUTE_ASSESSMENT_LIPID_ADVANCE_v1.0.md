# External-validation route assessment: LIPID / ADVANCE controlled-access dossier, v1.0

Date: 2026-10-28.

- **Scope:** an assessment of the project file "LIPID ↔ ADVANCE Dataset and Access-Route Dossier for Stage-A Compatibility and Stage-B Planning" (PDF, created 2026-08-13) against P4V v0.12.
- **Status:** assessment only. No request has been sent, no access has been granted, no value has been opened, and no protocol has been changed.

## Bottom line

1. **[Certain] The route cannot give a P4V confirmatory result.**
   - P4V v0.12's unit of analysis is the independent study, with one evaluation block per study.
   - LIPID and ADVANCE give N ≤ 2.
   - k_crit(n) := min{k : P(Bin(n, ½) ≥ k) < 0.05, k ≥ ⌈0.6n⌉} is **empty for n ≤ 4**, because P(Bin(4, ½) ≥ 4) = 0.0625. The primary therefore cannot reject at N ≤ 4 by construction.
   - Every tabulated planning scenario needs N_plan ≥ 30 (w = 0.80) to 60 (w = 0.70) independent studies.
   - One or two cohorts can support only a separately specified single-cohort analysis. Its result would be *supportive at most*, which is how P4U's two-study result was reported.
2. **[Certain] The dossier was written for Paper 2, which is CLOSED.**
   - It audits a *reference/target cohort pair* (`p2g_metadata_compatibility_audit.py`, "Stage B" cross-cohort transfer) in Paper 2's terms. Paper 2 was closed by `B0_CLOSURE_RECORD_v1.0.md`.
   - Its two draft emails describe that design: "reference/target cohort pair … candidate-generation procedure". They must not be sent as written for any P4V purpose.
3. **[Certain] LIPID is not one sample per participant as distributed.**
   - Verified today in the Nature Communications 2024 article, which describes LIPID lipidomics as "5991 participants for whom baseline and/or one-year fasting plasma samples were available", with "n = 5991 at baseline, n = 5782 at one-year followup".
   - The dossier reports only the 5,991. One sample per participant would hold only if the release carries participant ID and time point and the analysis is restricted to baseline, which must be prespecified before access. This is the same structure that left every MetaboLights candidate unresolved; the difference is that a custodian can supply the identifiers.
4. **[Certain] No record shows that either request was sent.** The dossier's manifest state is "route confirmed, application pending". The project holds no record of correspondence.

## What the dossier establishes, by evidence status

| Item | Dossier claim | Status here |
|---|---|---|
| LIPID access route | apply to Prof. John Simes / LIPID Study Investigators; analyses must not conflict with ongoing studies; "response … within two months" | **Verified** (Nature Commun. 2024 PDF, data-availability statement). Two months is a response time, not time-to-data. |
| LIPID samples | 5,991 participants, 342 species | **Verified with correction**: baseline 5,991 *and* one-year 5,782 (see 3 above). Panel size not checked here. |
| ADVANCE governance | George Institute policy: completed studies covered; proposal ≤ 3 pages; Data Sharing Agreement; requester ethics, possibly original-study ethics; cost recovery quoted after review; committee meets monthly; no turnaround; Institute representative expected as collaborator/co-author; manuscript to custodian ≥ 28 days before submission | **Verified** (georgeinstitute.org data-sharing policy page) |
| ADVANCE design | case-cohort 3,779 = random sub-cohort 3,154 + 625 outcome-enriched cases; 310 species | Relayed, not checked here |
| 3,154 indicator in an external release | unverified | Agree: unverifiable until the custodian answers |
| Shared profiling strategy (same lab and platform) | from JCI Insight 2018 | Relayed, not checked here |
| Feature overlap 265 / 342 / 310 | recomputed by the dossier's author | Relayed, not checked here. The artifacts are not in the project. |
| Artifact hashes | table of SHA-256 values | **Unusable**: every value in the PDF is 56 hex characters, not 64 (truncated in rendering). The files must be re-hashed from the originals before citation. The dossier itself flags its preregistration digest `1ba9d44e…62fa` as incomplete. |

[Likely] The dossier was produced by another assistant: its links carry `utm_source=chatgpt.com`. Its unverified rows should be treated as leads, not as evidence.

## Gate check against P4V v0.12

| Gate | LIPID | ADVANCE | Comment |
|---|---|---|---|
| Source admissible under §6 | no | no | §6 sources are (1) census studies, (2) post-census Metabolomics Workbench studies and (3) MetaboLights. A controlled-access cohort is route (c) of `TRACKB_P4V_STATUS.md` and needs v0.13. |
| B3 subject resolution | resolvable *if* the release has participant ID + time point; baseline restriction to be prespecified | baseline case-cohort, one sample per participant [Likely] | this is the only gate where controlled access is clearly better than repositories |
| Untouched provenance | no value ever opened by any branch [Likely: access never granted] | same | the feature lists (metadata) were inspected for Paper 2, which does not trigger provenance exclusion |
| Measurements | targeted LC–MS/MS, 342 species → five 60-feature hashed blocks | 310 → five blocks | [Likely] dimensionally compatible. Lipid-species RefMet naming must be frozen on the published list before access, and the targeted lipid panel is a different platform class from the Paper IV archive (B4). |
| Independence between the two | — | — | distinct participants, but the same laboratory strategy, so §6's shared-batch check applies if both are used |
| **Statistical sufficiency** | N = 1 | N ≤ 2 with LIPID | **fails**: k_crit empty for n ≤ 4; N_plan ≥ 30 |
| Outcome data | not needed | not needed | the P4V estimand is held-out prediction error among measured features, so clinical outcomes are not required. ADVANCE's 3,154 vs 3,779 question therefore concerns how the population is defined, not independence, and must be decided before the request. |
| Access burden | application; terms unspecified publicly | proposal, DSA, ethics, costs, co-authorship, 28-day review | time-to-data unspecified for both |

## Consequence

The owner's milestone was: *one externally independent cohort with verified participant structure, suitable measurements, and an authorised frozen predictive evaluation*. LIPID could plausibly meet the first two after access. The third would then be a **single-cohort protocol, not P4V**. Its claim is limited to one population, one laboratory and about five blocks, and P4V's own interpretation grid would grade it supportive at most.

## Decision A36 (owner): which claim is the external evaluation meant to support?

- **A36-1, confirmatory (P4V estimand kept):**
  - The cohort supply has to come from a route that can scale to N_plan ≥ 30–60 eligible, informative, independent studies.
  - The only candidate in hand is documentary recovery: 390 MetaboLights candidates, plus optionally the 104 Metabolomics Workbench PROXY_ONLY candidates from the B3 audit, subject to a v0.13 decision.
  - It needs v0.13 to admit documentary evidence as B3 route 2.
  - The first tranche then works as a futility test (`DOC_RECOVERY_PROTOCOL_DRAFT_v0.1.md` §7). Controlled-access cohorts cannot reach that N at the verified per-cohort burden [Likely].
- **A36-2, supportive single cohort:**
  - A new protocol (working name P4W) uses LIPID baseline as the cohort. It names its own unit and inference and is labelled supportive at most.
  - It must be frozen before the application text is written, because the application describes the analysis.
- **A36-3, both, in sequence (recommended):**
  - Send LIPID a *design-neutral feasibility inquiry* now. It would ask only about availability, participant identifiers, time-point metadata, terms and timelines, and make no analytical claims.
  - In parallel, freeze and run the documentary-recovery tranche.
  - Choose between A36-1 and A36-2 when both answers are in. The inquiry costs little and only starts the clock. It must be the owner's own email, drafted fresh, not the Paper 2 text.

**Do not:**
- send the dossier's emails as written;
- change any `access_status`;
- treat LIPID's 5,991 as one sample per participant without the time-point field;
- choose 3,154 vs 3,779, or remap lipid names, after any data are seen;
- describe a single-cohort result as P4V confirmation.

## Sources

- [Imputation of plasma lipid species to facilitate integration of lipidomic datasets (Nature Communications 2024), PDF](https://www.nature.com/articles/s41467-024-45838-3.pdf): data-availability statement and LIPID sample description, read 2026-10-28.
- [The George Institute for Global Health: data sharing policy](https://www.georgeinstitute.org/data-sharing-policy): read 2026-10-28.
- Project records: `TRACKB_P4V_PROTOCOL_DRAFT_v0.12.md` (§3 k_crit, §5c table, §6 sources); `TRACKB_P4V_STATUS.md`; `B0_CLOSURE_RECORD_v1.0.md` (Paper 2 closed); the dossier PDF.
