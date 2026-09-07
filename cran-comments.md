## Resubmission

IntegMultiReg was archived on CRAN on 2026-09-05 because the memory-access
problems reported by the supplementary Valgrind checks had not been corrected
before the deadline. This version corrects them.

### Fixes for the reported problems

* `main_function_prediction()`: the cross-validation status buffer is now
  initialized before it is copied in the binary and continuous outcome
  branches, removing the reported uninitialized read.
* Binary cross-validation prediction now frees the otherwise unused `yhat`
  allocation before returning, removing the reported 7,200-byte leak.
* The fitting and prediction entry points keep returned R objects protected
  across `PutRNGstate()` and diagnostic output.
* Candidate prediction arrays are initialized, and large input-dependent
  buffers were moved off the C stack to heap allocations that are freed after
  use.

## Test environments

* local: macOS 26.5.2, aarch64, R 4.5.2, gcc-15
* GitHub Actions: Ubuntu R-devel, Windows R-release, macOS R-release
* Ubuntu, R-devel (2026-09-04 r90492), Valgrind 3.22.0
* Ubuntu, R-devel, GCC UBSAN / Clang UBSAN / Clang ASAN+UBSAN
* macOS ARM, R-devel, ASAN+UBSAN (tests and vignette rebuild)

## R CMD check results

0 errors | 0 warnings | 2 notes

* New submission / Package was archived on CRAN. Expected for this
  resubmission.
* HTML version of manual: HTML Tidy is not recent enough and V8 is
  unavailable on the local machine. This note does not arise on Windows,
  which reported Status: OK.

Every platform ran the 338 test expectations with no failures, warnings or
skips.

## Memory checks

Installed tests under Valgrind (`--leak-check=full`, `--track-origins=yes`,
no package suppressions):

    [ FAIL 0 | WARN 0 | SKIP 0 | PASS 338 ]
    definitely lost: 0 bytes in 0 blocks
    indirectly lost: 0 bytes in 0 blocks
      possibly lost: 0 bytes in 0 blocks
    ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)

Examples and tests were additionally checked against a level-2 instrumented
R-devel build (r90498) with Valgrind 3.27.1, again reporting zero errors and
zero definitely, indirectly or possibly lost bytes.

As a control, the archived 0.1.1 tarball still reproduces its exact
7,200-byte / 60-block leak under these same settings, confirming that this
configuration does detect the defect that was reported.

Sanitizer checks passed on Linux and macOS ARM. The test suite includes a
regression test that exercises fitting, prediction and cross-validation under
frequent garbage collection.

## Other changes in this version

This version also adds a validated multi-platform data class, a formula/data
interface that preserves training transformations and factor coding at
prediction time, fitted-model diagnostics, and optional conditional
coefficient and posterior predictive intervals. Survival preprocessing now
logs positive event and censoring times once, matching the original
supplementary code of the reference paper; `survival_scale = "identity"`
reproduces the previous behaviour. See NEWS.md for the full list.
