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
* `sampler_method = "corrected"` adds the symmetric MRF conditional, the boundary
  Hastings correction and the corrected Gamma log-density rate sign. The
  default `"original"` retains the historical updates.
* Fits use a named four-section schema; `upgrade_imr_fit()` converts a
  structurally complete 0.1.x fit explicitly rather than automatically.
* `posterior_draws(latent = TRUE)` retains the augmented response draws for
  binary and right-censored fits, and `summary()` and `confint()` take
  `parm = "latent"`. The default is `FALSE` and the coefficient draws are
  bit-identical with and without it.

These are breaking changes. There are no reverse dependencies on CRAN.

## Checked source

The package remains version 0.2.0. The convention-name and tutorial snapshot is
`641ca7a5bdb5b4004442732c5359ddeffa0e972e`. The README pins that source for the
examples. Later documentation records the replay procedure and excludes a
local Git worktree pointer from source builds.

This snapshot changes public labels and saved-object conversion, and retains
the previous algorithm defaults. Its native sources, bundled data, namespace
and dependency declarations are unchanged from `9b552fa7a766ef01507414dc18f351629b2bddbe`.

## Current local checks

* R-devel (2026-09-21 r90579), GCC 11.2.0 and GSL 2.4 on x86_64 Linux:
  2904 installed expectations pass, with zero failures, warnings or skips.
* After mapping renamed metadata, 24 fit/prediction configurations and 18 CV
  results are bit-identical to an independent pre-rename installation.
* Two complete corrected-sampler covariate comparisons produce eight identical
  CSVs. Independent checks reconstruct outer MSE and verify that no outer test
  subject enters an inner fold. Both executed tutorials render successfully.
* The full local source check passes the tests, examples, rebuilt vignettes and
  PDF manual. It initially reports missing checkbashisms/tidy/V8 tooling and a
  bundled Git worktree pointer. Excluding that pointer leaves every other
  archived file byte-identical. A packaging recheck with checkbashisms passes;
  the existing website checks validate TeX with KaTeX independently of V8.
  The CI source checks supply tidy and V8.

The clean distribution archive has SHA-256
`e6d214298cac9ae8ea2dbbd8fbd324b47fc1635a1d565849b6b7400dd3c5487e`.
The validation record distinguishes this packaging-only change from the source
archive used by the full replication runs; all 147 retained files match.

## Earlier computational baseline

The computational baseline at `fe9ce3b7a286f2cee8f69ac0b8f25456212a79e1`
passed the cross-platform source checks, GCC/Clang and macOS ARM sanitizer
checks, and full Valgrind checks. The full-suite artifact contains 177 process
logs and the separate worker artifact 36, with zero reported errors or
unfreed definite allocations. These are identified as baseline checks rather
than relabelled as runs of the later naming revision. Native C sources are
unchanged by the naming revision.

## Numerical reference scope

The frozen ARM Mac manuscript reference remains unchanged and has its own
source-package and script checksums. Earlier independent Linux runs agree
with each other, while binary AUC differs from that Mac reference by up to
0.002924. The archived reference source itself reproduces this Linux behavior.
New source revisions and reference environments are recorded separately;
matching a version label or replacing an expected CSV is not evidence that an
older discrepancy has been explained.

This submission record is excluded from the source package by .Rbuildignore.
