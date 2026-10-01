## Submission

This is a feature release. The archived 0.1.1 memory-access problems were
corrected in 0.1.3, which is the version currently on CRAN and passing on all
flavors; nothing in this submission relates to that archival.

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

## Test environments

* GitHub Actions, `--as-cran` with `error-on = "warning"`:
  Ubuntu R-release, Ubuntu R-devel, Windows R-release, macOS R-release
* GitHub Actions sanitizers: Linux GCC UBSAN, Clang UBSAN, Clang ASAN+UBSAN,
  and macOS ARM ASAN+UBSAN, each running the installed tests and rebuilding
  the vignette
* GitHub Actions Valgrind on R-devel with level-2 instrumentation
* Linux, R 4.4.1, GSL 2.4: installed test suite, frozen default-regression
  comparison against the previous release, and independent sampler and
  original-C prediction checks

## R CMD check results

0 errors | 0 warnings | 0 notes

`Status: OK` on all four `--as-cran` platforms: Ubuntu R-release, Ubuntu
R-devel, Windows R-release and macOS R-release. Those runs use
`error-on = "warning"`, so a warning on any platform would have failed them.
The test suite reports `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 2864 ]`.

## Memory checks

Installed tests run under Valgrind with `--leak-check=full` and
`--track-origins=yes` and no package suppressions, and under ASAN+UBSAN on
Linux and macOS ARM.

Valgrind on R-devel reports:

    definitely lost: 0 bytes in 0 blocks
    ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

GCC UBSAN, GCC ASAN+UBSAN, Clang UBSAN and Clang ASAN+UBSAN on Linux, and
ASAN+UBSAN on macOS ARM, all pass with no sanitizer diagnostic.

As a control, the archived 0.1.1 tarball still reproduces its exact
7,200-byte / 60-block leak under these same settings, so the configuration is
known to detect the defect that led to that archival.

The test suite includes a regression test that exercises fitting, prediction
and cross-validation under frequent garbage collection, and post-fit
cross-validation now releases native resources and reports round, fold and
subgroup context when a Cholesky decomposition or solve fails.
