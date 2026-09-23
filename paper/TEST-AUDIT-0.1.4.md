# Test and API audit — IntegMultiReg 0.1.4

Date: 2026-09-08. Scope: public API calls, meaningful boundary cases,
statistical regression tests, editor-comment requirements and native-code CI.

## Assessment

The revised suite exercises every exported function and registered S3 method,
and tests the main fit–predict–validate workflow across binary, continuous and
right-censored responses. It is substantially stronger than a collection of
smoke tests: several results are checked against independent references or
invariance properties. It is not exhaustive over every argument combination,
every possible input or every statistical regime.

The audit found and fixed observable defects, not just missing assertions:
integer conversion could warn/overflow; duplicate feature names could silently
lose columns; a custom identifier could collide with a feature named `id`;
and some invalid saved-fit iteration counts and subgroup indices were accepted.
These supplement the 0.1.4 CV leakage, probability-averaging, concordance and
theta-label corrections described in REVISION-0.1.4.md.

## Measured local results

- 467 passing expectations; 0 failures, warnings or skips.
- 14/14 exported functions and 21/21 registered S3 methods executed under covr.
- Instrumented line coverage: 89.42% overall (baseline 89.06%), 91.28% in R,
  87.63% in C. Coverage measures execution, not branch completeness or accuracy.
- Full source-package `R CMD check`: Status OK, including examples, rebuilt
  vignette and PDF reference manual; 0 errors, warnings and notes.
- Environment: macOS ARM, R 4.4.2, GSL 2.8.

The smallest change in overall coverage should not obscure the value of boundary
regressions: many invalid inputs share already-executed validation lines. The
previously uncalled `print.summary.imr_data` is now explicitly tested.

## Public API and S3 test mapping

Paths below are relative to the package's `tests/testthat/` directory. The
machine-readable coverage artifact `public-api.csv` lists each of the 35 entry
points separately. `.github/scripts/coverage.R` derives that list from NAMESPACE
and fails if a future exported function or registered S3 method is unexecuted.

| API / method family | Main evidence |
| --- | --- |
| `imr` default, formula and `imr_data` methods | `test-fit.R`, `test-outcome-types.R`, `test-data-interface.R`, `test-formula-prediction.R`, `test-api-boundaries.R` |
| `imr_data`, `validate_imr_data`, data print/summary/extraction | `test-data-interface.R`, `test-api-boundaries.R` |
| `validate_imr` | `test-diagnostics.R`, `test-api-boundaries.R` |
| `compare_imr` | `test-diagnostics.R` |
| fit print/summary and `summary_imr` | `test-methods.R` |
| fit coefficients and `coef_imr` | `test-methods.R` |
| `predict_imr`, `predict.imr` | `test-predict.R`, `test-formula-prediction.R`, `test-review-regressions.R` |
| `cv_imr` | `test-cv.R`, `test-review-regressions.R`, `test-api-boundaries.R` |
| `posterior_summary`, its print method, `confint.imr` | `test-diagnostics.R`, `test-review-regressions.R` |
| `posterior_draws`, posterior print/summary/coef/confint/predict | `test-posterior-draws.R` |
| `plot_imr`, `plot.imr` | `test-methods.R`, `test-plots.R` |
| `plot_subgroup_sizes`, `plot_top_features` | `test-plots.R` |

## Boundary and numerical evidence

| Risk | Tested property or reference |
| --- | --- |
| Invalid computational controls | Missing/nonfinite, empty, nonscalar, fractional and overflowing integers fail before sampling; combined retained + burn-in counts cannot overflow native integers. |
| Minimum supported chain | One retained draw with zero burn-in fits and predicts for all three response families. This verifies dimensions, not chain convergence. |
| Data alignment and naming | Duplicate/missing IDs, reordered IDs, nonnumeric/nonfinite features, duplicate/blank/missing column names, identifier collisions and corrupted availability metadata fail or align as specified. |
| Formula design and new data | Factors, contrasts, polynomial bases, custom IDs, unseen levels and unsupported offsets/no-intercept models have explicit behavior. |
| Corrupted saved fit | Invalid subgroup indices, missing draws, inconsistent matrix dimensions and iteration-trace lengths are rejected. |
| Prediction | No routed subjects, platform/covariate mismatches, row ordering, response scales and full numerical precision. |
| CV leakage | Perturb held-out outcomes: their predictions do not change. Perturb held-out features: the training fit does not change. Full-fit selection/latent draws are not reused. |
| CV metrics | MSE reconstructed from out-of-fold records; tied/censored concordance compared against `survival::concordance`; single-class AUC, empty scores and no comparable survival pairs return NA. |
| Binary averaging | Average model-specific probabilities checked against an explicit two-model numerical oracle, distinguished from applying the link after averaging. |
| Theta uncertainty | Six explicitly labelled pairs for a platform in four subgroups match the native draw order. |
| Posterior uncertainty | Response/mean prediction distinction, residual uncertainty, model-selection zero mass, response bounds/scales, conditional diagnostics and seed reproducibility. |
| R/native memory boundary | Frequent garbage collection exercises fit, prediction and CV; separate sanitizer and Valgrind CI instruments native execution. |

## Editor-comment coverage

All seven response sections are mapped in `response-to-editor.tex`: data class,
fit methods/uncertainty/covariates, executable manuscript code, the prediction
table discrepancy, software outside R, full-precision prediction, and formula
syntax. Clinical covariates remain forced effects; automatic clinical-variable
selection through a new inclusion prior is a model change and is not claimed.
The new installed example compares prespecified formulas on identical folds and
uses nested CV for predictive formula selection. Related-software
positioning is a manuscript review item, not something unit tests can prove.

The replication script includes all nine manuscript R input blocks, prints
outputs and checks versioned numerical references. Default full replications
are distinguished from quick smoke runs. The historical 0.1.3 archive is not
used as validation of the new CV algorithm.

## Remaining limits

- Execution coverage is not branch coverage; 100% correctness cannot be inferred
  from 35 exercised API entry points or 467 expectations.
- Rare numerical failure/cleanup branches and legacy utility paths remain only
  partially covered (`src/utils.c` has 62.80% local line coverage). Exhaustive
  allocation-failure and very-large-dimension stress tests were not performed.
- Short deterministic fixtures do not establish convergence or high-dimensional
  frequentist calibration. Conditional diagnostics do not assess selection-chain
  convergence, and model weights retain a Laplace approximation.
- The obsolete full-scale KIRC performance table has been removed from the
  manuscript and current replication entry point and archived separately.
  No new full-scale KIRC performance study is claimed.
- The full manuscript simulation tables were replicated in the stated local
  environment. Cross-platform CI exercises package checks and tests, not the
  entire full-default manuscript simulation; cross-platform table identity is
  not claimed.
- CV uses fixed settings. The formula-selection example supplies nested
  validation; other adaptive choices must also remain inside training folds.

GitHub run identifiers, final status and any platform-specific notes are recorded
in `output/revision-0.1.4/CI-RESULTS.md`; all validation claims must be read with
that source revision and scope.

Current addendum: `../output/editor-completion-0.1.4/COMPLETION.md`. Earlier CI
evidence remains in the revision-0.1.4 output directory.
