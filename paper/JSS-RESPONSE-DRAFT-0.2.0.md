# Response to the JSS editorial assessment — working draft

**Not ready for submission.** This draft maps the editor's points to the current
implementation and manuscript. Final page references, full-scale analyses, candidate distribution checks
and remaining validation must
be completed before removing this status. It does not claim that the editor's
concerns have been accepted or closed.

## 1. Software contribution beyond the original method's supplementary code

We have reorganized the presentation around a reusable analysis workflow for
partially overlapping platforms. The package provides validated data objects,
list and formula interfaces, three outcome families, inspection and diagnostic
methods, uncertainty summaries, model comparisons, prediction and configurable
cross-validation. The revised account distinguishes the original mathematical
method, conventions in the released implementation, and package extensions.
Options share the same fitting and prediction engines rather than duplicating
an implementation for each convention.

Locations: manuscript Introduction, Analysis strategies, Data workflow, and
Choosing analysis settings (`sec:intro`, `sec:strategies`, `sec:workflow`,
`sec:settings`); package `inst/METHOD-COVERAGE.md`. These additions expand the
software's scope; they do not themselves establish reproduction of the published
historical tables.

## 2. Dedicated validated data type and newdata representation

The `imr_data` class stores outcomes, platform data, covariates and subject
availability. Construction/validation checks IDs and dimensions, while print,
summary and extraction methods expose the structure. Fitting and prediction
reuse this representation, reducing repeated manual alignment of data frames.

Locations: manuscript Data workflow (`sec:workflow`), package data-interface
methods and tests. The complete installed candidate suite includes the
`data-interface` checks; see the local validation evidence below.

## 3. Inspection, comparison, covariate choice and uncertainty

The fitted-object section documents validation, summaries, selection plots,
interaction/log-posterior traces and fit comparisons. The covariate example
uses paired and nested refit CV to compare formula specifications. Clinical
covariates remain forced predictors within a given specification; this example
is not a claim that they undergo within-model selection.

Conditional coefficient and predictive intervals are available through
`posterior_draws()`. The manuscript explains their conditional-chain and Laplace
selection-weight limitations. Explicit `initial` selection/interaction matrices
now support different chain starts. Short interface tests and component sampler
checks are not convergence evidence for the original-scale fits.

Locations: `sec:methods`, Coefficient and predictive intervals,
`sec:covariate-comparison`, and `sec:diagnostics`.

## 4. All manuscript code and printed output in the replication material

`execute-manuscript.R` extracts and executes the manuscript's CodeInput blocks,
including `print(fit, threshold = 0.5)`, and saves individual code/output files,
a command map and combined HTML. The older duplicate `make_analysis.R` is now
a forwarding entry point to the maintained replication script; its original
implementation is archived.

Location: `sec:reproducibility`, `paper/execute-manuscript.R` and
`paper/REPLICATION-README.md`. All nine displayed blocks have now executed
successfully against the shared-history candidate. The eight nonempty outputs
match the manuscript; missing outputs for the covariate comparison and reduced
KIRC summary were inserted, and explicit console width keeps them readable.
The audit saves source/library hashes and runs relative-output commands in an
isolated workspace. See [current displayed-code and PDF audit](../output/manuscript-accepted-20260921/PROGRESS.md).
The latest 34-page working draft has been visually reviewed after the structural
and paragraph-cohesion edits; its exact source/PDF hashes are recorded in the
[latest source and visual review](../output/citation-audit-20260921/PROGRESS.md).
All displayed code and outputs were preserved through those edits. The draft
reports R 4.4.2 for the computational environment. Further edits require renewed
checks; executing displayed code does not by itself reproduce all scientific
figures and tables.

## 5. Binary AUC discrepancy between supplied and reproduced tables

We acknowledge the reported discrepancy. The historical scripts and reference
CSV files have been retained rather than overwritten. The exact causal chain
for the old reported values has not been uniquely identified because the
original run's complete build/configuration manifest is unavailable. We do not
attribute it conclusively to a newly identified sampler difference.

The revised workflow records settings, source/data hashes, folds, full-precision
outputs and versions. Two complete manuscript runs of the same validated
shared-history candidate pass independent score and nested-fold audits. All 23 CSV/RDS artifacts agree at
tolerance 1e-12 (21 are R-identical), and three figures are pixel-identical at
150 dpi. See [current replication status](../output/manuscript-full-20260921/PROGRESS.md). The accepted values now have a separate reference in
`paper/expected-results-0.2.0-paper/`, with package and evidence checksums.
They are not a retroactive explanation of the historical mismatch.

Locations: `sec:reproducibility`, `paper/historical-audit/`, historical
`expected-results-0.1.4/`, and handoff comparison records.

## 6. Related functionality outside R

The revised overview includes MOFA2/mofapy2 (R/Python), mvlearn (Python),
Similarity Network Fusion (MATLAB/R), and MultivariateStats.jl (Julia), with
references to project documentation or the methods publication. It compares
modelling scope rather than asserting absence of competing capabilities or
claiming an unperformed cross-software performance benchmark.

Locations: Introduction and `tab:cross-language`. These four cross-language
entries have been checked against primary documentation and the SNF paper;
see the [source audit](../output/citation-audit-20260921/PROGRESS.md).
The wider bibliography still requires its final consistency review.

## 7. Full-precision prediction return values

Predictions are returned and written to numerical artifacts at full precision;
rounding belongs to display formatting. Prediction and CV tests exercise this
separation. No historical CSV was rounded or refreshed to make a mismatch pass.

Locations: `sec:predict`, prediction methods and associated precision tests.

## 8. Standard formula/data regression interface

The formula interface uses the standard R terms/model-frame/model-matrix
machinery, retains factor levels and contrasts for prediction, and supports
transformations. List/data-class routes remain available for multi-platform
inputs. Training-fold preprocessing is rebuilt for refit CV, and paired/nested
examples demonstrate comparisons of covariate specifications.

Locations: `sec:workflow`, `sec:covariate-comparison`, and formula/prediction/CV
tests. Formula interface support does not imply that all CV modes are refit;
the manuscript separately describes post-fit and training-refit estimands.

## Current evidence and final release gates

Local ordinary and ASAN/UBSAN suites for the shared-history candidate each
passed 2812 assertions. Frozen default regression covers 156 fields across six
fits and eighteen CV runs, all passing at tolerance 1e-10 (115 bitwise identical).
A same-chain 50,000-draw code2017 comparison has bitwise-identical posterior,
model, preprocessing and complete CV output, while peak RSS decreased from
1,754,890,240 to 673,693,696 bytes. All retained states and their order remain
available; no thinning or model dropping supplies this reduction. Lossless
research checkpoint recovery and lightweight fold metadata support resumption.
Separately scoped component checks address sampler mathematics and controlled
original-C predictions. These scopes are not interchangeable. macOS
LeakSanitizer is disabled; the separate strict Linux Valgrind evidence is
described below.
Current-candidate GitHub Actions now pass all four cross-platform R checks
(Status: OK), all four Linux GCC/Clang sanitizer combinations and macOS ARM
ASAN/UBSAN, including its vignette rebuild. Each separate full-suite verification
passes 2812 assertions. Full logs are retained in the CI evidence directory.

The earlier local source-package check had 0 errors, 1 warning and 3 notes;
those local diagnostics remain recorded alongside the subsequent cross-platform
Status: OK checks. The 72 injected allocation-failure recovery cases and
sanitizer vignette rebuild were also refreshed successfully.

The fixed candidate subsequently passed strict Linux Valgrind validation
(run 35843055431), using R-devel r90579 compiled with level-2 instrumentation
and Valgrind 3.27.1. Its source-package check passed 2765 assertions and rebuilt
the vignette; the only check NOTE concerns example runtime under instrumentation.
A separate instrumented parent/PSOCK suite passed 2811 assertions. The example,
test, vignette and independent-parent logs have zero Memcheck errors under
the archived graphics suppressions, as do all 274 completed worker logs.
Fork/exec logs without final summaries are counted separately.

The graphics rules were checked against package-free controls, including nine
scoped replays. Negative controls detected the old package's 7200-byte leak
and an independent R-heap uninitialized-read probe. The historical old-package
uninitialized-read stack was not reproduced in this environment; the probe
does not establish that reproduction. These checks support native-memory
validation of the candidate, while the scientific convergence and full-study
requirements below remain separate.

The original-scale code2017 Table 1 run uses 350,000 retained states, 50,000
burn-in and 10 folds x 10 rounds. All five model comparisons have completed
and passed independent checks of settings, folds, subject coverage, scores and
summary uncertainty. This establishes consistency of the documented rerun,
not recovery of historical digits or chain convergence. The clinical Cox
benchmark issued 200 warnings in two subgroups whose Stage 2 subjects had no
events. We retain the unpenalized benchmark and flag the separation problem;
its extreme coefficients are not reliable hazard-ratio estimates.

The original-scale simulation resource run has completed one replicate in each
of the three main scenarios. All predictive results and Bayesian selection AUCs
passed validation. Separate baseline refits reproduced the L1-CPH selection
paths and both L1-CPH and Uni-CPH selection AUCs. A truncated inner glmnet path
still covered the selected lambda; the warning remains recorded. A univariate
Cox warning also recurred on refitting. Agreement with saved values does not
establish the validity of that fit's Wald inference. These checks support the
executed manifest, not the full simulation study, which requires 50 main and
30 correlated replicates (Web Appendix E, Figure 5).

All eight paper-sampler chains completed with distinct recorded starts and seeds
at 350,000 retained states after 50,000 burn-in. Checkpoint integrity and exact
diagnostic projections passed, but the joint diagnostics show slow mixing.
The methylation interaction has classical PSRF 1.114 (upper bound 1.243) and
rank-normalized split/folded Rhat 1.201, with bulk ESS 27 and tail ESS 85.
The miRNA interaction has rank-Rhat 1.028. Feature-mPIP correlations range from
0.557 to 0.945, with a maximum absolute difference of 0.531; some methylation
top-ten lists share only two features. These findings prevent a convergence or
stable-ranking claim from the original iteration budget. A separate, fixed-budget
eight-chain rerun at 1,000,000 retained states per chain is underway, preserving
original evidence. Its trajectories and joint diagnostics remain to be checked.

Appendix G selection rankings can be exported by platform and subgroup, with
original feature names and the published top-list sizes. Its coefficient columns
are explicitly posterior modes. The released code computes model-conditional
coefficient vectors but does not supply a coefficient-table export or establish
the table-specific conditioning and aggregation. We therefore do not substitute
mPIPs or posterior mean coefficients for the published mode estimates.
Proprietary IPA networks are not regenerated.

Current evidence links:

- [Latest manuscript structure and visual review](../output/editorial-cohesion-20260921/PROGRESS.md)
- [Accepted numerical reference and earlier PDF](../output/manuscript-accepted-20260921/PROGRESS.md)
- [Full-run acceptance, warning audits and strict CI follow-up](../output/ci-repair-20260922/PROGRESS.md)
- [Current candidate GitHub Actions](../output/github-actions-20260921/PROGRESS.md)
- [Shared-history validation and resource comparison](../output/history-sharing-20260921/PROGRESS.md)
- [Current package/fault/vignette check](../output/shared-candidate-check-20260921/PROGRESS.md)
- [Lossless checkpoint and resume checks](../output/checkpoint-sharing-20260921/PROGRESS.md)
- [Original-scale code2017 run status](../output/original-full-code2017-20260921/PROGRESS.md)
- [Completed eight-chain review](../output/appendix-eight-chain-review-20260923/PROGRESS.md)
- [Interim six-chain mixing and ranking audit](../output/appendix-six-chain-interim-20260922/PROGRESS.md)
- [Eight original-scale diagnostic chains](../output/appendix-chains-full-20260921/PROGRESS.md)
- [Original simulation scale check](../output/original-simulation-scale-20260921/PROGRESS.md)
- [Original appendix design audit](../output/appendix-audit-20260921/PROGRESS.md)
- [Appendix G selection and coefficient audit](../output/appendix-g-audit-20260921/PROGRESS.md)
- [Candidate and initial-state validation](../output/initial-state-20260921/PROGRESS.md)
- [Prior sanitizer/fault checks and limits](../output/sanitizer-20260921/PROGRESS.md)
- [Pilot and clinical Cox warning analysis](../output/pilot-20260921/PROGRESS.md)
- [Random-marker implementation and checks](../output/marker-design-20260921/PROGRESS.md)
- [L1-CPH warning audit](../output/glmnet-audit-20260921/PROGRESS.md)

A separate clean-unpack candidate material snapshot installs successfully and
runs the full manuscript workflow. Independent score and holdout checks pass;
all 23 numerical artifacts agree with the accepted reference run at tolerance
1e-12, and the three figures are pixel-identical at 150 dpi. Input/reference
hashes are unchanged. This validates the candidate distribution's manuscript
workflow, while final scientific and release gates remain open. See the
[clean-unpack evidence](../output/candidate-materials-check-20260922-v2/PROGRESS.md).

Before submission: complete original-scale experiments and multi-chain checks,
run final-distribution package checks and all manuscript code,
verify numerical references in clean directories, inspect the compiled PDF,
resolve or explicitly scope all remaining limitations, and assemble/test the
final source and replication archives. The candidate has been pushed to a
separate GitHub testing branch under the user's authorization. No main-branch
merge, package release or journal submission has been made.

### Unpenalized original-data limitation found in the larger pilot

A 20,000-draw paper-sampler IMR clinical+molecular pilot on the bundled KIRC
inputs fails unpenalized all-draws CV. In subgroup 011, 22.065% of retained states
exceed the 64 training rows in five folds, and 28.060% exceed the 62 rows in the
other five folds. An independently reconstructed failing design has 69 columns,
64 rows and rank 64. These fractions are dimensional lower bounds; remaining
designs were not exhaustively rank-tested. A ridge=.001 control completes but
changes the estimator. Do not silently add a penalty, drop states or replace the
inverse when claiming the unpenalized paper reference. The separate code2017
sampler/ranked-model configuration requires its own evidence. Details and saved
artifacts: `output/memory-scale-20260921/PROGRESS.md` in the research workspace.


### Generated-outcome confirmation of the unpenalized CV limitation

A separate paper-sampler/random-marker simulation pilot completed scenario 1
but failed all-state, zero-ridge CV in scenario 2 at 20,000 retained draws.
Independent design reconstruction found 25.2%–30.0% of subgroup 011 states
with more columns than training subjects, including a 64-by-65 matrix of rank
64. Thus this is not limited to the original KIRC outcome. Enlarging the
replicate count does not repair an undefined solve. Regularization, state
truncation or training-refit validation would be separately labelled alternatives,
not silent implementations of the original evaluation. The checkpoint, failed
attempt and reconstructed design are retained in the
[simulation feasibility evidence](../output/paper-simulation-feasibility-20260922/PROGRESS.md).
