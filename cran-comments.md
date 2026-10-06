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

## Checked source

The package remains version 0.2.0. The implementation was validated at
`f90a7d350f9c68fa745d397594c08ee08a505ded` in
https://github.com/h7nian/IntegMultiReg/pull/6.
The subsequent website changes do not alter R/C code, the namespace, data,
help topics, tests or vignettes. This submission record is excluded from the
source package by .Rbuildignore.

## Test environments and R CMD check

All thirteen candidate checks completed successfully on October 6, 2026:

* Ubuntu R-release and R-devel, Windows R-release, and macOS ARM R-release:
  each R CMD check log reports `Status: OK` (0 errors, 0 warnings, 0 notes).
  Linux and macOS also check the PDF manual; Windows checks without it.
* The installed suite reports 2880 passing expectations on each platform,
  with zero failures, errors, warnings or skips.
* Linux GCC UBSAN, GCC ASAN+UBSAN, Clang UBSAN and Clang ASAN+UBSAN pass.
* The macOS ARM ASAN+UBSAN job passes.
* Coverage, documentation-export verification and the complete covariate
  comparison workflow pass.

The earlier local HTML-tooling note is absent from these CI checks, which
provide tidy and V8. The local R 4.4.1 run used an older testthat and reported
2879 expectations; that count is not substituted for the CI result.

## Memory checks

The full installed suite and instrumented PSOCK workers pass Valgrind.
The parent log reports:

    definitely lost: 0 bytes in 0 blocks
    indirectly lost: 0 bytes in 0 blocks
    possibly lost: 0 bytes in 0 blocks
    ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

The retained full-suite artifact contains 177 logs with zero error summaries;
the dedicated worker artifact contains another 36. These are the current
candidate's results, not results copied from a previous release.

## Numerical and documentation checks

Independent installations of the pre-refactor and candidate sources give
identical results across six fits and eighteen CV combinations on the same
runtime. Native Gamma-density, prior-indexing and stationary-distribution
checks pass, as does comparison with independently compiled archived-C
conditional predictions at tolerance 1e-12.

Two independent manuscript-scale predictive runs on Linux give identical
saved fit/CV objects for all six outcome/model configurations. The eight
covariate-comparison CSV tables also agree with an independent execution of
the displayed manuscript commands. The recorded ARM Mac numerical reference
is preserved: its binary MCMC results are platform-sensitive, and the archived
source reproduces the current Linux values. No reference table was overwritten
to claim cross-platform bitwise agreement.

The website is built from maintained documentation, with generated HTML
excluded from the source package. Local-link, desktop search, math rendering
and mobile-navigation checks pass at https://h7nian.github.io/IntegMultiReg/.
