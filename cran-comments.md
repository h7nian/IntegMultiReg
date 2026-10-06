## Submission

This is a feature release. The archived 0.1.1 memory-access problems were
corrected in 0.1.3; this release adds the interface and analysis features
described below.

## What changed

0.2.0 reorganizes the user-facing interface and adds the analysis options the
package previously lacked. The changes that reviewers are most likely to look
at are:

* `imr()` is now an S3 generic with list, formula and `imr_data` methods, and
  the duplicated `predict_imr()`, `summary_imr()`, `coef_imr()` and
  `plot_imr()` wrappers are removed in favour of the standard generics.
  Legacy argument aliases are deliberately not accepted. `inst/MIGRATION.md`
  gives the complete old-to-new table.
* `cv_imr()` gains `cv_method`, and independent `ridge`, `model_set`,
  `df_method`, `score_method`, `folds` and `fold_rng` arguments. The previous
  behaviour is the default.
* `sampler_method = "paper"` adds the symmetric MRF conditional, the boundary
  Hastings correction and the corrected Gamma log-density rate sign. The
  default `"legacy"` retains the historical updates.
* Fits use a named four-section schema; `upgrade_imr_fit()` converts a
  structurally complete 0.1.x fit explicitly rather than automatically.
* `posterior_draws(latent = TRUE)` retains the augmented response draws for
  binary and right-censored fits, and `summary()` and `confint()` take
  `parm = "latent"`. The default is `FALSE` and the coefficient draws are
  bit-identical with and without it.

These are breaking changes. There are no reverse dependencies on CRAN.

## Current candidate and validation

The package remains version 0.2.0. The local source tarball checked on
October 5, 2026 has SHA256
`fd0fdbe1a311090c70869759f9772cd2c16852540d467dccd98d53898dc7aa08`.
All results below refer to that candidate's R/C sources and generated help.

* Linux, R 4.4.1, GSL 2.4: the installed suite passes 2879 expectations,
  with zero failures, errors, warnings or skips. The test driver checks errors
  separately from failed expectations and exits unsuccessfully for either.
* Independent installations of the preceding source snapshot and the refactor
  give identical posterior objects, preprocessing, model metadata, predictions
  and CV metrics across six outcome/model combinations and eighteen CV runs.
* The native Gamma-density, prior-indexing and stationary-distribution checks
  pass. Conditional predictions agree with the independently compiled archived
  C reference to the recorded tolerance of 1e-12.
* Roxygen 7.1.1 regenerates NAMESPACE and all help files without changes.
  The new documentation guard also rejects a deliberately misplaced helper
  that takes over a public function's export annotation.

## R CMD check results

Linux, R-devel 2026-09-21 r90579, `R CMD check --as-cran`:

0 errors | 0 warnings | 1 note

The note reports unavailable HTML checking tools (`tidy` and R package `V8`).
Examples, including `--run-donttest`, package tests, vignette rebuilding and
PDF manual generation pass. Remote CRAN incoming checks were disabled for
this local run.

## Native memory check

Valgrind 3.27.1 with level-2 instrumented R-devel runs a focused check of the
refactored native paths: binary, continuous and censored fits with both sampler
conventions, prediction, and all three CV modes. It reports:

    definitely lost: 0 bytes in 0 blocks
    indirectly lost: 0 bytes in 0 blocks
    possibly lost: 0 bytes in 0 blocks
    ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

This focused check does not replace the complete sanitizer and Valgrind suites.
The previous cross-platform and memory results belong to earlier sources.
Before submission, run those suites and the four-platform checks on the final
candidate, resolve the HTML-check tooling note, and update this record with
that candidate's source identity and results.
