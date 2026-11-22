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
