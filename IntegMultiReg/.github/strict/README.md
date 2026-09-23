# Instrumented CRAN-style release audit

Official basis: [Writing R Extensions, 4.3.2](https://cran.r-project.org/doc/manuals/r-devel/R-exts.html#Using-valgrind).
The manual documents `--with-valgrind-instrumentation=2`, `--track-origins=yes`,
and `R CMD check --use-valgrind` (including vignette rebuilds since R 4.2).

This manually dispatched release audit supplements the ordinary per-push
Valgrind and cross-platform jobs. It builds Valgrind 3.27.1 from its official
source and R-devel r90498 with level-2 instrumentation and reference BLAS.
It checks the exact 0.1.3 submission tarball, including `--run-donttest` examples,
tests, vignette code/rebuilds and the PDF manual. No package suppression is used.

The archived CRAN 0.1.1 tarball must reproduce its exact 7,200-byte/60-block
allocation leak. Its historical uninitialized-read report did not reproduce
with the Ubuntu compiler; this limitation is recorded explicitly, not reported
as reproduced. A separate deliberately invalid R-heap read verifies level-2
detection capability, but does not reproduce that package call stack. A random
failure or missing dependency does not count as a successful control.

CLI_NO_THREAD=1 disables the test dependency cli's optional progress timer,
using its upstream switch, to avoid a 336-byte loader/TLS allocation at exit.
This does not disable package computation or filter Valgrind diagnostics.
The initial controls and their unsuccessful acceptance are retained in the
audit report. Candidate acceptance remains zero errors across all stages.

Fixtures are excluded from CRAN source packages by the existing `.Rbuildignore`
rule for `.github`. SHA-256 pins both fixtures and the Valgrind download. The
candidate verifier checks packaged tracked files against the checkout, allowing
only R build's line-ending/final-newline normalization of Makevars.win; generated
DESCRIPTION metadata is treated separately. Refresh the candidate and checksum
before auditing a future package version. Do not reuse a prior green run for
new source. The R download hash, configuration proof, versions, old failure
log and full candidate check tree are uploaded even if a check fails.

This follows official techniques but does not claim identical CRAN hardware,
compiler flags, system libraries or all private CRAN settings.
