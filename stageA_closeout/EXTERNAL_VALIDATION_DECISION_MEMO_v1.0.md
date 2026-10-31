# External validation after T1: decision memorandum v1.0

Date: 2026-10-31. Status: advice to the owner. It changes no protocol, authorises no access, and opens no data.
Labels: **CERTAIN** = verified from records or reproduced; **SUPPORTED** = from records or sources not re-verified here; **HYPOTHESIS**; **UNKNOWN**.

## 1. The uncomfortable answer

Under the evidence standard the owner set (affirmative documentation that each participant contributes one sample), **every repository route to P4V's N ≥ 30 independent studies has now been tested and closed.** The only remaining route that could reach that N is a consortium of prospective cohorts, run federated. Its feasibility is unknown and will take months to establish.

Until then, the manuscript's third requirement (independent external validation) can be met only at the **supportive** level (P4U; possibly P4W), not the confirmatory level.

## 2. What is now established

| Route | Record | Result | Status |
|---|---|---|---|
| Metabolomics Workbench metadata (B3) | `P4V_B3_EVIDENCE_AUDIT_RESULT_v1.0.md` | 0 demonstrated / 104 proxy-only | CERTAIN (record); EXECUTION NO-GO 2026-10-25 |
| MetaboLights census (3,582 accessions) | `STAGEA_CENSUS_RESULT_v1.2.md` | L1 strict = L1 sensitivity = 0 | CERTAIN (record) |
| MetaboLights documentary recovery, stratum M | `DOCREC_T1_RESULT_v1.0.md` | x = 0/40; 95% upper bound 25/371; DR-7 reached | CERTAIN |
| Controlled-access cohorts (LIPID, ADVANCE) | `EXTERNAL_ROUTE_ASSESSMENT_LIPID_ADVANCE_v1.0.md` | N ≤ 2; k_crit empty for n ≤ 4 | CERTAIN |

**T1 also gives the reason** [CERTAIN]. Publications do not state one sample per participant, and 30% of sampled candidates document the opposite. More documentary tranches, or applying the same documentary route to the 104 Workbench proxies, would meet the same reporting convention [HYPOTHESIS, strongly suggested by T1].

## 3. What P4V tests, and what the development evidence says

- **P4V's hypothesis** [CERTAIN, P4V v0.10–v0.12 §4]. In an independent cohort, the own-design theorem pair beats a borrowed theorem pair on the same held-out operator, in aggregate across studies. The primary is a one-sided sign test on one win indicator per study.
- **P4V does not test** average superiority of a new acquisition selector over G1.
- **Selector development evidence** [SUPPORTED, project decision log]:
  - The predictive-superiority programme was closed after Gate 3A, where no predictor family qualified against G1.
  - NOVA-S reached +0.11% over G1 against a 0.29% qualification floor.
- **Consequence.** A new external selector-versus-G1 test would be confirming something that development has not shown. The externally testable claim with prior support is P4V's directional-geometry prediction.
- **Prior support for that claim** [SUPPORTED, Paper IV]. P4U was supportive: two studies, donor-randomisation p = 0.014.

## 4. Remaining routes, ranked

Criteria: (a) can reach N ≥ 30 independent study units; (b) participant independence establishable; (c) measurement compatibility; (d) access and time; (e) evidentiary class; (f) provenance risk.

| Rank | Route | (a) | (b) | (c) | (d) | (e) | (f) |
|---|---|---|---|---|---|---|---|
| 1 | **Federated consortium (COMETS)** | plausible: > 80 prospective cohorts, each ≥ 100 participants with metabolomics [SUPPORTED] | custodians hold participant IDs; baseline restriction can be prespecified [HYPOTHESIS] | platform heterogeneity (MS vendors, NMR) may confine the study to one platform family and cut N [HYPOTHESIS] | nomination by a Steering Committee member, proposal vote, then per-cohort uptake; months [SUPPORTED process; time UNKNOWN] | confirmatory, if v0.13 is frozen before any cohort result | untouched [HYPOTHESIS] |
| 2 | **Publish the documentation audit** (Stage A + B3 + T1) | n/a | n/a | n/a | ready now | a separate, complete contribution; does not test the hypothesis | none |
| 3 | **LIPID single cohort (P4W)** | N = 1 | baseline restriction with custodian IDs | targeted lipid panel | inquiry approved (A36-3), unsent; two-month response time | supportive only | untouched |
| 4 | **A35 verification** | ≤ 19 | documentary, as T1 | — | already fetched | completes a prespecified arm; cannot reach 30 | none |
| 5 | Author contact, stratum M | needs a large recovered pool; 30% already contradicted; recovered studies still face compatibility, provenance and informativeness gates | author assertion | — | slow; low and unknown response | UNKNOWN | — |
| 6 | Further documentary tranches; documentary route on Workbench proxies | DR-7 reached; same reporting convention | — | — | — | — | — |

### Why the consortium route fits P4V better than it first appears [HYPOTHESIS, to verify against v0.12 §3]

- **Little crosses cohort boundaries.** The primary needs from each study one win indicator Z_s. The only cross-study inputs are design-derived feature pairs (feature identities, not values).
- **A two-round federated execution** would work like this:
  1. Each cohort computes its own-design pair locally on its design partition.
  2. A coordinator applies the frozen donor-matching rule and returns each cohort its borrowed pair.
  3. Each cohort evaluates both pairs on its untouched evaluation partition and returns Z_s plus diagnostics.
- **No participant-level data leaves any cohort.** That is the governance model consortia already use for summary-level meta-analysis.
- **Main technical risk.** A borrowed pair must exist in the recipient's feature set, which requires harmonised metabolite identities across platforms.

## 5. Recommendation and GO/NO-GO criteria

1. **Accept the T1 record** (§9 of the record).
2. **Send the LIPID feasibility inquiry this week.** It is already approved and design-neutral, and it only starts a clock. Any result can be supportive at most.
3. **Make one design-neutral COMETS feasibility contact** (owner's own email; no protocol attached). Ask four things:
   - whether a federated, summary-only predictive evaluation run by member cohorts with a frozen container is admissible;
   - how a non-cohort investigator obtains nomination;
   - how many member cohorts profile ≥ 100 baseline participants on a common platform family;
   - typical time from approval to cohort results.
   - **GO** to drafting P4V v0.13 (consortium cohorts as study units, custodian-attested independence, federated execution) only if the answer makes ≥ 30 compatible, participating cohorts realistic after v0.13's own power and attrition calculation.
   - **NO-GO** otherwise: confirmatory P4V is closed for lack of data supply, and the programme reports supportive external evidence only. This would be a prospective, pre-outcome amendment; no external P4V outcome has been opened, so it does not replace a test after seeing results.
4. **Write the documentation-audit paper now.** Working claim: *public metabolomics repositories do not document participant independence well enough for study-level external validation.* The evidence:
   - the census: 3,582 accessions, L1 = 0;
   - a randomised documentary audit: 0/40 recoverable; 30% documented dependence; duplicate deposits in the pool;
   - the Workbench audit: 0/104.
   - A prior-art check against published repository-compliance assessments is needed before any novelty claim [UNKNOWN].
5. **A35.** Complete V1–V4, which is cheap and prespecified, if the audit paper goes ahead; otherwise record it as not executed after T1, without reinterpretation.
6. **Stop:** further documentary tranches (DR-7); author contact for stratum M (low, unknowable yield against a large required pool); any documentary route on the Workbench proxies.

## 6. The manuscript decision this forces

There are two defensible options for the Springer *Metabolomics* manuscript:

- **(i) Submit now.** Report external validation as supportive (P4U; P4W if LIPID data arrive), and document the confirmatory infeasibility as a finding.
- **(ii) Hold for the consortium route.** Its feasibility is unknown and the timeline is at least months.

What is not defensible is describing any single-cohort or two-study result as confirmation, or relaxing the independence standard to manufacture N.

## Sources

- [NCI DCEG: Consortium of Metabolomics Studies (COMETS)](https://dceg.cancer.gov/research/who-we-study/cohorts/comets): > 80 cohorts, > 2 million participants; membership, nomination and proposal process. Read 2026-10-31.
- [The Consortium of Metabolomics Studies (COMETS): Metabolomics in 47 Prospective Cohort Studies (Am J Epidemiol 2019)](https://pubmed.ncbi.nlm.nih.gov/31155658/): title only; abstract not retrieved.
- Project records cited in §2–§3.
