# NOVA external validation: preregistration stub v0.1

Programme: decision-optimal metabolite acquisition (NOVA-S vs G1). Date: 3 November 2026. Author: Anas Enoch.
Status: **hypothesis and access-decision record, frozen by hash; not yet an executable protocol.** Development
archive (186 records): closed after 11 documented uses. External routes: LIPID and COMETS, inquiries sent 31 Oct 2026.
**Data-access decision deadline: 15 December 2026.**

## 1. Evidentiary status
Under P4V the independent study is the replication unit. One external cohort gives supportive evidence only (P4W),
never confirmatory evidence of generalisable selector superiority. This stub concerns NOVA-S vs G1 acquisition
performance; it does not substitute for, and does not inherit the confirmation rule of, P4V's directional-geometry
hypothesis. The 0.29 % floor of Gate 2 was the development programme's detectable effect at 512 clusters; it is
reported alongside the external result as context and is not the external success criterion.

## 2. Primary hypothesis
At budget B = 1, NOVA-S yields lower held-out predictive loss than G1 on an independently sampled external cohort:
Δ_ext = 1 − L_NOVA-S / L_G1 > 0, both losses pooled over the same evaluation cells under the frozen procedures.

## 3. Frozen implementations (to be copied by hash, not re-described)
- G1: Paper IV plug-in rule (linear shrinkage α = min(1, max(0.05, d/(d+n))), Schur complements, ridge 1e-6), as
  re-estimated per training split in `p6_n2.py`.
- NOVA-S: `p6_n2.py` v0.1.0 SHA-256 `681ce852609e6ab4fcaf0f0ff603477c64e2351e10bd119811d45f008d4bbe8e` under
  `N2_PROTOCOL_v1.0.md` `6bbafb3c…`: RBF kernel-ridge nuisance on log-rule inputs, N2 grids, 2-SE paired departure
  rule (τ = 2), floored squared loss. No variant.
- Block construction (20 O / 20 T / 20 P), target standardisation and loss: the P4D rule and salt as applied to the
  186 archive records (`design_model_manifest.json`, `master_salt` PAPER3_OPEN_BENCHMARK_V1_20260928), applied once
  to the external cohort's design partition. The same blocks are used by both arms.

## 4. Primary evaluation
- Participant-level sampling structure must be established from the cohort's own documentation before any
  partition is made; all samples of one participant fall in one partition (LIPID: baseline and 1-year samples of
  a participant are one cluster, never two participants).
- Design/evaluation partition by the frozen P4D rule; evaluation rows never touched before the freeze of §7.
- Both arms receive identical design rows, blocks, preprocessing and tuning budget; no evaluation information in
  selection or tuning.
- Uncertainty: participant-cluster bootstrap (10,000 resamples), seed fixed in the executable protocol.
- **Success criterion (supportive), fixed now:** one-sided 97.5 % lower bound of Δ_ext > 0 **and** Δ_ext ≥ 0.10 ×
  Δ_G1→∅ (ten per cent of G1's own measured gain over no measurement in the same cohort and blocks). Failure:
  otherwise. [Owner to confirm the 0.10 factor before hashing; it must not be changed after any outcome is seen.]
- Reported with the result, whatever it is: Δ_G1→∅, the agreement rate of the two selectors, the departure rate,
  and the archive floor 0.29 % as context.

## 5. Secondary hypothesis (unvalidated; exploratory on any cohort)
G1 selection instability identifies decisions where an alternative pick improves prediction. The N2 diagnostic
(archive use 11, `N2_REGIME_DIAGNOSTIC_RESULT_v1.0.md`) found no record-level stability association; its A4′
G1-plausible-set stratification is invalid for inference because the reference set used folds whose training
included the scored fold. Any external evaluation needs a fold-nested stability measure built without the scored
fold. No stability-gated policy is qualified; if one is to be tested, its algorithm and rule are frozen before §7,
else the analysis is labelled exploratory.

## 6. Data-access decision, 15 December 2026, each route separately
GO: usable data access secured, with documentation sufficient to establish participant independence and to judge
feasibility of the frozen task (panel overlap, targets, n). NO-GO: no usable access or documentation by the
deadline. A positive reply without usable access is NO-GO. No further repository-search cycle is authorised.
If both are NO-GO: external validation is unperformed; the programme's outputs are the NOVA negative-result
manuscript and the documentation-audit manuscript; this stub stays in the repository as the prediction a future
cohort would test.

## 7. Freeze requirements before any external analysis
SHA-256 of the executable protocol, selector code, preprocessing spec, block rule, evaluation code and this stub,
recorded before evaluation rows are read. Any change after outcomes are accessible makes the analysis exploratory
and is labelled so. Negative results are reported without adjustment of selector, subgroups or criteria.

## 8. Interpretation
A positive result in one cohort is supportive; a negative one is reported as such. Neither alters any frozen
archive number. The empirical claim of superior acquisition remains unproven until an independent test under
these rules is performed.
