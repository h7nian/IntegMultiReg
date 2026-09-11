# IntegMultiReg 0.2.0

## Breaking API cleanup

* `imr()` is now an S3 generic with list, formula, and `imr_data` methods. Its
  first argument is `x`; fitting arguments now use `covariates`, `outcome_type`,
  `min_subgroup_size`, `draws`, `burnin`, `molecular_prior_scale`,
  `forced_prior_scale`, `residual_prior`, and `interaction_prior`.
* Method values are lowercase (`"imr"`, `"bms"`). Printed output keeps the
  familiar uppercase labels. Legacy argument aliases are intentionally not
  accepted.
* Predictions use the column name `prediction`. `cv_imr()` now returns the
  named fields `pooled`, `fold_mean`, `predictions`, `metric`, and `validation`.
* Removed the duplicate `predict_imr()`, `summary_imr()`, `coef_imr()`, and
  `plot_imr()` wrappers. Use the standard S3 generics.

## Fit schema and native safety

* New fits use schema version 2 with four named sections: `control`, `model`,
  `preprocessing`, and `posterior`. `upgrade_imr_fit()` explicitly upgrades a
  structurally complete 0.1.x fit; public methods do not upgrade automatically.
* Formula fits retain only the identifier, response, and variables used by the
  formula. Fit validation now checks native-bound types, dimensions, mappings,
  indices, and draw counts before compiled code is called.
* MRF normalization now enumerates states with an unsigned counter and uses
  log-sum-exp. A platform is limited to 16 modelled subgroups, with an R error
  before native code for larger models.
* Native entry points are named `imr_fit`, `imr_predict`, and
  `imr_concordance`; unreachable legacy numerical helpers were removed.

See `inst/MIGRATION.md` for a complete old-to-new API table.

# IntegMultiReg 0.1.4

* Reject platform names `id` and `subgroup` in data construction and validation
  to prevent collisions with availability metadata and misleading CV errors.

* Add a runnable paired-CV and nested formula-selection example, with saved
  subject/fold assignments and repeatability checks in CI.

* Reject overflowing integer controls, ambiguous data-frame column names, and
  corrupted saved-fit iteration counts or subgroup mappings before computation.
* Add public API execution auditing and boundary regressions; CI retains coverage
  and native sanitizer logs as downloadable artifacts.

- Replaced the post-fitting CV approximation with independent training-fold
  MCMC fits, preprocessing and model weights. Removed the native CV path that
  used held-out responses and assigned increasing weight to larger errors.
  Runtime now scales with folds times rounds; old prediction tables must be
  regenerated. CV preserves the caller's RNG and exposes out-of-fold records.
- Formula fits retain raw formula data so CV can rebuild transformations using
  only training subjects. Older formula fits require refitting for CV.
- Corrected survival concordance for tied predictions and incomparable pairs;
  undefined fold metrics return NA.
- Binary point prediction now averages model-specific probit probabilities.
- Corrected theta interval labels when a platform appears in four or more groups.
- Added deterministic leakage, probability-averaging and pair-order regressions,
  with concordance comparisons against survival when available.

# IntegMultiReg 0.1.3

* Explicitly initialize conditional outcome pointers for strict compiler diagnostics.

- Corrected survival preprocessing to log positive event/censoring times once,
  matching the original supplementary C code. Survival fits must be rerun;
  `survival_scale = "identity"` explicitly reproduces the historical raw-time
  implementation. The CV partitioning and selection sampler order are unchanged.
- Added optional `posterior_draws()` with subgroup coefficient intervals and
  posterior predictive intervals. All active coefficients, including intercepts
  and clinical effects, use the original pMOM prior. Clinical effects remain
  always included; automatic clinical variable selection is not implemented.
- Conditional sampling reaugments binary and censored responses and reports
  classical split R-hat. Model averaging retains the fitted selection sampler's
  Laplace approximation; it is not an exact model-weight calculation.
- Survival posterior prediction uses time-scale medians as point summaries
  because inverse-gamma variance mixtures need not have finite time-scale means.
- Added independent analytic/integration checks and survival-scale regressions.

# IntegMultiReg 0.1.2

## Native-code reliability

- Fixed premature unprotection of R return objects in the fitting and
  prediction entry points.
- Fixed an allocation leak in binary cross-validation prediction.
- Initialized cross-validation work arrays and moved input-sized buffers off
  the C stack.
- Added dedicated valgrind CI and macOS ARM vignette checks under ASAN/UBSAN.

## User interface

- Added the validated `imr_data` class for training and prediction inputs.
- Added an `imr()` formula/data method while preserving the original interface.
- Formula predictions reuse training transformations, factor levels and
  contrasts; subject IDs are excluded from dot-formula expansion.
- Unsupported no-intercept and offset formulas now fail explicitly instead
  of silently fitting a different model.
- Added fitted-object validation, fit comparison, posterior uncertainty
  summaries, credible intervals and parameter-specific trace plots.
- `predict.imr()` now returns values at full numeric precision.

## Reproducibility

- Expanded the standalone replication script to execute every R command shown
  in the manuscript and to fail when full-run reference results drift.
