# External relevance factor: owner confirmation record v1.0 (5 Nov 2026)

**Purpose.** `EXTERNAL_VALIDATION_PREREG_STUB_v0.1.md` (sha256 `d1cb1f82…`) left one element of its success criterion
open: "[Owner to confirm the 0.10 factor before hashing; it must not be changed after any outcome is seen.]" This
record confirms that factor prospectively and separately, so that the stub is preserved unedited with its history
intact. This record does not modify the stub; it supplies the confirmation the stub requested.

## 1. The confirmed criterion

The external success criterion stands as written in stub §4:

> one-sided 97.5 % lower bound of Δ_ext > 0 **and** Δ_ext ≥ 0.10 × Δ_G1→∅ (ten per cent of G1's own measured gain
> over no measurement in the same cohort and blocks)

with the baseline condition that Δ_G1→∅ has a one-sided 97.5 % lower bound > 0.

**Owner decision, 5 November 2026 (Anas Enoch):** retain the factor 0.10. It is not recalibrated, and in particular
it is not replaced by the 0.29 % archive detectability floor. Owner's stated reason: a detectable-effect floor is not
automatically a scientifically meaningful effect threshold, and the relevance threshold must not be lowered because
the current method struggles to reach it.

## 2. Evidence the factor was confirmed against

The confirmation was made with the following development figures in view, all from frozen records:

| Quantity | Value | Source |
|---|---|---|
| NOVA-S development improvement Δ_acq | +0.112 % [+0.029, +0.221] | `SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md` §1 |
| Hindsight maximum over the threshold family τ ∈ {0,…,4} | +0.158 % at τ = 1.5 | matrix v0.4, replay #10 |
| Archive detectability floor (Gate 2 MDE, K = 512) | 0.29 % | matrix v0.4 |
| Δ_G1→∅ measured on the one external cohort examined | 15.8 % [LB 13.1 %] | `NHANES_FA_RESULT_v1.0.md` |
| Implied bar on a cohort with that baseline gain | 1.58 % | 0.10 × 15.8 % |

The mismatch between the implied bar and the development figures is recorded as substantial and acknowledged. The
confirmation is therefore made in the knowledge that the development evidence does not demonstrate an advantage of
the required magnitude. It is not a prediction that the criterion will be met.

**Correction of record.** An earlier assessment by Claude stated that the criterion "cannot be passed at any sample
size." That claim was too strong and is withdrawn here: the criterion is a ratio and is therefore cohort-dependent,
and both Δ_G1→∅ and the attainable nonlinear advantage vary with the data-generating distribution. The defensible
statement is the row-by-row comparison above, on a cohort with a NHANES-like baseline gain.

## 3. Timing attestation

- Confirmed on 5 November 2026, before any Phase II / public-matrix candidate was screened and before any new
  metabolite measurements were accessed. The public-cohort screening cycle closed at commit `bc6aba2` with
  **no eligible cohort and no measurements accessed**; NHANES FAS_H/DEMO_H were never downloaded.
- No external or candidate outcome has been observed under this criterion. Nothing in this record may be changed
  once any such outcome is seen.

## 4. Open item requiring a separate owner decision (NOT settled by this record)

Because the criterion is a ratio, Δ_G1→∅ is part of the success definition, and it is computable from a cohort's
**design partition alone**, with no NOVA-S involvement and no access to evaluation rows. Cohort eligibility could
therefore drift toward cohorts with small Δ_G1→∅, which lowers the bar without touching NOVA-S, any outcome, or any
evaluation data. Such drift would pass every leakage control the programme currently has.

Given the figures in §2, the two conditions are jointly satisfiable only in a narrow band — roughly
Δ_G1→∅ ≲ 1.0–1.6 % if the attainable advantage is of the order of the development estimate — and that band is
reachable by selection.

**Proposed mitigation, for owner decision, to be fixed before any cohort is examined:** a prespecified eligibility
window for Δ_G1→∅, justified on scientific grounds about what constitutes a meaningful acquisition problem, rather
than leaving the denominator unconstrained. This record does not set that window. Until it is set, no cohort
eligibility decision that depends on Δ_G1→∅ may be made.

Status of this item: **OPEN**. Raised by Claude; not an owner decision; no value proposed here.

## 5. Scope

This record confirms one numerical factor in one frozen criterion. It makes no claim about whether NOVA-S meets it,
authorises no data access, and changes no implementation. The stub's hypothesis, estimands, frozen implementations,
partitioning and inference are unchanged.

Author of record: Anas Enoch. Compiled by Claude (Opus 5) from the frozen records cited above.
