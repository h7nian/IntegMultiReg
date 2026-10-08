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
