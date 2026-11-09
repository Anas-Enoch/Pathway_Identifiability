# From structural fidelity to decision value: acquisition programme synthesis

9 November 2026. Retrospective synthesis based on versioned study records
and manuscript summaries. It adds no experiment or inferential result, makes no
claim of publication or novelty, and changes no frozen decision. The selection
and evaluation archives were not reopened.

The [scope addendum](ACQUISITION_SUCCESSOR_SCOPE_v1.1.md) connects this synthesis
and the translation audit to the original successor scope record.

The programme asked which extra metabolite to measure when only part of a panel
can be observed. It established why the answer depends on the prediction task,
how some scoring rules can agree or disagree, and where evidence about the value
of a replacement breaks down. The completed work leaves G1 (A0) as the comparator
and the selector-improvement slot empty. That position reflects the tests and
decisions below; it is not a claim that G1 is unbeatable.

## From a structural score to a prediction decision

**Paper I separated fidelity from utility.** The commutator score admits an
exact algebraic reduction. For the rank-one form in the revised manuscript,
its vanishing has an invariant-direction interpretation. The score's ability
to recover a representation-derived ordering did not establish an advantage
on independent classification or reconstruction endpoints. The disposition
record archives Paper I as a negative audit and records the label and oracle
provenance limitations. The lesson is that reproducing an oracle defined by
the same representation cannot stand in for measuring downstream value.
[Paper I disposition](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/PAPER_I_DISPOSITION.md); [algebra and downstream manuscript](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/FINAL_MANUSCRIPT.md).

**Paper II moved the endpoint downstream.** The exploratory Stage A work
separated the observed, candidate and target panels. Candidate self-prediction
and overlap between selection and scoring were identified as concrete sources
of optimism. In ST001849 lipid partitions, the realised gain landscape was
substantially less reliable than its Gaussian counterpart. It was possible to
narrow the candidate pool without reliably naming its best member. The
subsequent Stage B branch was dormant, not a completed negative experiment.
The Stage A planning results should remain exploratory, including the overlap
qualifications recorded at handover.
[Paper II closeout](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/README_paper2_stageA.md); [Stage B status](</Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/paper3_feasibility 2/stageB_status_dormant.md>).

**Paper III tested policies on held-out reconstruction.** Its budget-one primary
compared the Gaussian selector against marginal variance in reference-standardized
RMSE. PAPER3_REML_CORRECTION.md reports the converged REML estimate across 186
benchmark records as +.00867 (95% interval [.00659, .01075]). The full interval
falls below the .014 practical threshold in p3_protocol_preregistration.md §3.
The numerical correction is part of
the provenance. The estimate is not a confirmatory population statement beyond
the benchmark; its sampling assumptions need to be read alongside the later
documentation audit. Paper III did not execute the dormant Stage B hypothesis.
[Paper III estimand](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/p3_protocol_preregistration.md); [corrected result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/PAPER3_REML_CORRECTION.md).

**Paper IV explained objective dependence.** For a common Gaussian rank-one
update, G1 values trace reduction and RIG ranks determinant reduction, which
is equivalent to conditional Gaussian information gain. Strict disagreement
certifies Loewner incomparability: the preferred candidate can change with
the quadratic task weighting. This is a population result for the specified
model and objective. It does not say an estimated selector will win held out.
RIG superiority was not supported on the reused 186-record split. The separate
stability analysis had an association with reproducibility, with provenance
exceptions and a Gaussian calibration that limit biological interpretation.
The two-study directional external result was supportive, not confirmatory
and not evidence of selector superiority.
[Paper IV manuscript and qualifications](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/paper4_v1.1/source/paper4_manuscript_v1.1.md).

## What the attempted improvements added

**NOVA distinguished predictive improvement from acquisition improvement.**
The programme separated the Gaussian linear score, population Bayes acquisition
value and the gain of the fitted prediction procedure. Better prediction from
the observed panel does not imply a better choice of extra measurement. The
selective NOVA-S rule's +.112% development gain fell below the .29% development
qualification scale; threshold replay did not repair that comparison. The
certification work also exposed the gap between covering fitted-procedure gain
and covering population Bayes value. The NHANES 2011–12 development exercise
gave agreement in all 25 block-fold decisions and therefore no discriminatory
contrast. It did not test superiority on the untouched 2013–14 cycle.
[NOVA synthesis](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_PAPER_v0.4.md); [selector ledger](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md); [NHANES result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NHANES_FA_RESULT_v1.0.md).

**NOVA-D tested recovery where an advantage was planted.** The frozen rule
made no departures in either designated condition (0/500 fold decisions in
each; folds within replicate are dependent), so both recovery claims failed.
Near-tie screening did not yield a qualifying improvement on that grid, and
the validity controls passed. The versioned corrections distinguish a gain
over the estimated comparator from recovery against a population comparator.
The later operative-shrinkage diagnostics showed that joint regularization
can reorder population Gaussian scores. That is a statement about the score,
not a causal attribution of every finite-sample error. This diagnosis motivated
G1-R as a different experiment.
[NOVA-D corrected record](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_RESULT_v1.2.md); [operative-shrinkage note](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_EXACT_ORDERING_NOTE_v1.2.md).

**G1-R tested decision loss directly.** In its Gaussian family, the selector's
regret is computed from the population covariance, without an evaluation
predictor. The matched conditional control reproduces joint shrinkage once the
conditioning ridge, variance floor, scale and clipping are accounted for. That
equivalence verifies the mechanism's representation, not a performance gain.
The decoupled GCV arm supported H1 at D1; H2 at D2 was not established, and its
interval did not establish harm beyond the margin. Both were required. Thus
the joint criterion is unmet. Grid averages and sparsity descriptions do not
replace the designated comparisons.
[G1-R result v1.1](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/G1R_RESULT_v1.1.md).

## External evidence supply and the successor decision

The **documentation audit** addresses whether public records establish the
sampling units needed for an independent comparison. The census subsection of
audit_paper/audit_results.tex reports that 797 of the 3,582 MetaboLights
accessions reached the prespecified independence test and zero qualified.
DOCREC_T1_RESULT_v1.0.md §1 reports 0/40 strict recoveries in a random sample of
the 371 unresolved stratum-M candidates. DOCREC_A35_RESULT_v1.0.md §1 reports
3/19 passing V1–V4 pre-verification in the separately flagged A35 pool. The A35
record leaves strict documentary recovery unassessed. The Workbench
factor-metadata check covers another source and sampling frame. Each finding
concerns evidence available in the documentation under its stated rules.
[Census results](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/audit_paper/audit_results.tex:3),
[T1 results](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/DOCREC_T1_RESULT_v1.0.md:6),
[A35 results](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/DOCREC_A35_RESULT_v1.0.md:6).

## Claim ledger

The evidence column distinguishes internal participant-grouped evaluation where
the grouping is documented, tests whose independent study is the unit of
replication, and external generalisation to a new population. Theory, synthetic
experiments and documentary audits have their own evidence types. A study-indexed
analysis can carry unresolved participant structure and reused outcomes. The
selector evidence matrix supplies the history of benchmark reuse and the
qualification decisions. [Selector evidence matrix v0.4](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md).

| Claim | Status | Evidence level and scope | Record |
|---|---|---|---|
| The Paper I structural surrogate improves held-out classification or reconstruction. | Not established | Internal held-out audit in the studied cohorts. | [Paper I disposition](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/PAPER_I_DISPOSITION.md) |
| A single winning candidate was reliably identifiable in Paper II Stage A. | Not established | Exploratory within-cohort split-half assessment in ST001849. The landscape narrowed the candidate pool. The external Stage B test remains open. | [Paper II closeout](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/README_paper2_stageA.md) |
| The Gaussian selector clears the .014 practical RMSE margin over variance at budget one. | Not established | Study-level benchmark analysis of 186 records. The corrected interval is [.00659, .01075]. The later documentary audit qualifies claims about independent sampling units. | [Paper III correction](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/PAPER3_REML_CORRECTION.md:50) |
| RIG improves held-out RMSE over G1. | Not established | Internal held-out comparison on a reused evaluation split, analysed across records. | [Paper IV abstract](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/paper4_v1.1/source/paper4_manuscript_v1.1.md:9) |
| The external directional comparison yielded donor-randomization p = .014. | Established within that comparison | Supportive two-study evidence, with one study added after preregistration. The endpoint was directional agreement. | [Paper IV abstract](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/paper4_v1.1/source/paper4_manuscript_v1.1.md:9) |
| The frozen selectors agreed in all 25 NHANES 2011–12 development block-folds. | Established | Internal participant-disjoint folds in a documented survey cohort. The 25 decisions share participants and analytes. | [NHANES development result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NHANES_FA_RESULT_v1.0.md:18) |
| NOVA-S clears the .29% development qualification scale. | Not established | Reused development records. The observed gain was +.112% on the stated endpoint. | [Selector evidence matrix v0.4](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/SELECTOR_CLAIM_EVIDENCE_MATRIX_v0.4.md:113) |
| NOVA-D recovers the planted advantage in its designated tests. | Not established | Frozen synthetic evaluation with nonlinear alternatives and Gaussian controls. The gate made zero departures in both designated conditions. | [NOVA-D result v1.2](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_RESULT_v1.2.md) |
| Frozen G1 joint shrinkage changes the population Gaussian argmax in 27 of 72 conditions at operative training sizes. | Established | Population-level calculation in the NOVA-D synthetic family using exact population covariances and the frozen ridge. The argmax agrees across both actual fold training sizes in each condition. | [Exact ordering note v1.2](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_EXACT_ORDERING_NOTE_v1.2.md:39) |
| G1-R's A2 meets both the D1 improvement and D2 no-harm criteria. | Not established | Frozen Gaussian synthetic study. H1 supported and H2 not established. | [G1-R result v1.1](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/G1R_RESULT_v1.1.md) |
| Zero of 797 screened MetaboLights records met the census independence criterion. | Established in the census | Documentary evidence under the prespecified eligibility rules. | [Census results](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/audit_paper/audit_results.tex:3) |
| Zero T1 records passed strict documentary recovery. | Established in the sample | Random documentary assessment, 0/40. A35 pre-verification has a different endpoint. | [T1 result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/DOCREC_T1_RESULT_v1.0.md:6) |
| A successor yields a scientifically relevant gain over A0 in a new population with acceptable harm. | Open | External generalisation. Requires a declared measurement decision, a loss bridge and independent evaluation with verifiable sampling units. | [Successor scope record](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/ACQUISITION_SUCCESSOR_SCOPE_v1.0.md) |

## Mathematical results and assumptions

| Result | Assumptions and interpretation | Record |
|---|---|---|
| **Signed approximation gap**, $V_j-G1_j=A_j-A_0$. | Classical squared-loss projection and subtraction. $A(S)$ is the excess risk of the best affine predictor over the Bayes predictor with inputs $S$. The targets are square-integrable and $O\subseteq O\cup\{X_j\}$. Here $A_0=A(O)$ denotes an approximation gap. Each gap is nonnegative and their difference can have either sign. The proposed $H_j=V_j-G1_j\ge0$ assertion is withdrawn. Joint Gaussianity gives equality of the acquisition values. | [NOVA draft, Theorem 2](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_PAPER_v0.4.md:107) |
| **Loewner incomparability under strict selector disagreement.** | Classical rank-one update and matrix-order argument applied to one-step Gaussian acquisition. With a positive definite residual target covariance $K_{TT}$ and admissible updates, trace and determinant rankings use the $I$ and $K_{TT}^{-1}$ weightings on the same update. Strict reversal supplies opposite trace signs under those positive weights and forces the update difference to be indefinite. | [Paper IV derivation](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/paper4_v1.1/source/paper4_manuscript_v1.1.md:64) |
| **Joint shrinkage equals a matched ridge-conditioned score.** | Schur-complement algebra with unit-diagonal standardization and $0\le\alpha<1$. The observed-block ridge becomes $(\alpha+\epsilon)/(1-\alpha)$, the residual variance floor becomes $\alpha/(1-\alpha)$, and the score has a common factor $1-\alpha$. A denominator clip is rescaled by the reciprocal factor. The population checks found reorderings at the operative training sizes. The identity uses the Schur form. | [Equivalence derivation](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/g1_shrinkage_equivalence.py:2), [population ordering checks](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_EXACT_ORDERING_NOTE_v1.2.md) |
| **OLS residual-Wishart and conditional pair law.** | The residual scatter $Q\sim W_{\nu}(K)$, with $\nu=n-p-1$, is a classical consequence of projecting Gaussian errors off an intercept and $p$ regressors with full column rank. With positive definite residual covariance, $\nu\ge2$ and fixed target metric, conditioning on the two-candidate scatter gives a matrix-normal cross block and an exact known-parameter law for the paired score difference. The confidence problem requires treatment of unknown parameters and the selected reference. | [Residual and pair-law derivations](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/EXACT_WISHART_PAIR_NOTE_v0.1.md) |
| **Variance-null fixed-pair power ceiling.** | Classical Neyman–Pearson bound applied to a programme-specific equal-gain null. In the fixed-metric Gaussian model with free candidate residual variances, a change of two variances can equalize the gains while preserving total residual risk. With $\nu=389$ and level .05, the D1-like fixed-pair ceiling is about 7.42%. Different variance constraints or a permitted-harm null change the comparison. | [Boundary-null derivation and arithmetic](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/CERTIFICATE_POWER_SCOPE_v0.1.md) |

## Evaluation coverage still open

- **NHANES 2013–14** is the only untouched, publicly available evaluation reserve with participant identifiers documented in the programme records reviewed here. CDC documents linkage of FAS_H and DEMO_H through the survey participant key SEQN. The recorded fatty-acid exercise left the cycle unopened, preserving it for a future independently justified and prospectively declared test if its feasibility is established. The 0/25 agreement figure refers to 2011–12 development. [NHANES result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NHANES_FA_RESULT_v1.0.md:3), [CDC fatty-acid codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/FAS_H.htm), [CDC demographic codebook](https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2013/DataFiles/DEMO_H.htm)
- **Cross-cohort selector transport** still needs a valid held-out risk comparison under verified assay alignment and participant linkage. Paper IV records a transfer screen with identical picks and exploratory risk computations under unverified unit conversions. The directional external comparison addressed a separate endpoint. [Paper IV transfer record](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/paper4_v1.1/source/S3_Text.md)
- **Non-Gaussian regularization** remains open for the G1-R decoupled-versus-frozen joint criterion. G1-R tested a Gaussian family. NOVA-D tested synthetic nonlinear acquisition alternatives under a separate protocol. [G1-R result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/G1R_RESULT_v1.1.md), [NOVA-D protocol](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_PROTOCOL_v0.3.md)
- **External nonlinear selector superiority** under the relevance criterion remains unevaluated. Completed work outside the 186-record archive includes synthetic tests and the NHANES development feasibility exercise. [NOVA-D result](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NOVA_D_RESULT_v1.2.md), [NHANES feasibility decision](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/NHANES_FA_RESULT_v1.0.md:35)

## Why the successor slot remains empty

The **scope record** consolidates why a standalone positive-gain certificate
line was deprioritized. The covariance-region reversal witness and elementary
F-pivot width concern particular constructions. The variance-null power bound
concerns a fixed-pair comparison with unknown candidate variances. The
approximately .00034 D1 aggregate screen plugs published sample summaries into
an optimistic bound; it is not a guaranteed population ceiling or dominance
over A2. Together they justify prioritization without proving impossibility.
[Scope record and supporting-note ledger](/Users/nour/Desktop/Springer/Pathway_Identifiability/github.pathway/ACQUISITION_SUCCESSOR_SCOPE_v1.0.md).

The [translation audit](ACQUISITION_EXTERNAL_TARGET_CHARTER_v0.1.md)
adds a units check to this decision. The recorded external relevance factor
uses G1's remaining loss for the successor contrast and no-acquisition loss
for G1's own gain. At the recorded 15.8% baseline gain the bar is 1.58% of
G1's loss, or 1.33036% of no-acquisition loss. The synthetic averages make a
case for caution but do not identify an eligible or ineligible population.
No external action has been identified in the reviewed records whose value
justifies opening a successor study on the strength of those averages.

The defensible closeout is an evidence chain, not a ranking of universal
winners: structural fidelity is weaker than independent utility; the desired
utility determines the population optimum; estimating that optimum introduces
decision risk; and independent validation requires verifiable sampling units.
Any future candidate enters the same prospective average-regret comparison
with an absolute expected-harm cap. A named external decision, loss bridge,
condition distribution, protected class and externally justified budgets must
precede that comparison. None is declared by this retrospective report.
