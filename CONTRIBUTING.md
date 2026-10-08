# Contributing to IntegMultiReg

## Local checking tools

The CRAN checking tools `checkbashisms` and `qpdf` are not runtime
dependencies. Package users do not need them. Maintainers running
`R CMD check --as-cran` locally can install them with
`brew install checkbashisms qpdf` on macOS or
`sudo apt-get install devscripts qpdf` on Debian/Ubuntu.

## Source layout and development checks


The public API uses `snake_case` names and S3 methods. Internal R functions use
`.imr_`; `%||%` is the internal null-default operator. Statistical priors and simulation controls use
`imr_priors()` and `imr_mcmc()` respectively. Their resolved values are
stored in the fitted object's `control` component.

| Responsibility | Source |
| --- | --- |
| Validated data and subject identifiers | `R/imr-data.R` |
| Fit dispatch and sampler inputs | `R/fit.R` |
| Formula response and design matrices | `R/formula.R` |
| Subject alignment and matrix scaling | `R/preprocessing.R` |
| Prediction and new-data routing | `R/predict.R` |
| CV folds, training fits and PSIS reweighting | `R/cv*.R` |
| Fit printing, summaries and feature rankings | `R/methods.R` |
| Fit validation and MCMC diagnostics | `R/diagnostics.R` |
| Plotting and shared appearance settings | `R/plots.R`, `R/plot-theme.R` |
| Joint chain preparation and stored draw access | `R/joint-sampling.R`, `R/posterior.R` |
| Shared chain/fold worker orchestration | `R/workers.R` |
| Joint sampler and pairwise concordance | `src/joint_sampler.c`, `src/concordance.c` |

Native calls list their arguments explicitly so reviewers can compare each call
with its registered C signature. The independent reference calculations in the
tests deliberately retain their own implementations; sharing the production
calculation would remove the comparison they are intended to make.

After editing roxygen comments, regenerate the documentation with
`roxygen2::roxygenise()` and inspect the `NAMESPACE` and `man/` changes. Run
`Rscript .github/scripts/verify-documentation.R` to check that the source
annotations reproduce the recorded exports. This catches an internal helper
accidentally inserted between a public function and its documentation.

Install the candidate with `R CMD INSTALL --install-tests .`, then run
`Rscript .github/scripts/verify-tests.R`. This checks the installed package and
fails on test errors, failed expectations, warnings or skips. Build the source
package and run `R CMD check --as-cran` on that tarball. Record the source
revision or snapshot hash with each result; a previous candidate's CI run does
not validate later changes. The short help examples demonstrate the interface;
their MCMC budgets are not evidence of convergence.


The joint sampler is checked against independent small-model integrals with
`Rscript .github/scripts/verify-joint-reference.R`. Maintained integration and
quadrature sources are in `.github/reference/`; their frozen fixtures are
independent of production transition code. Regenerate fixtures only with the
reference calculations, never from MCMC output. Compare the new marginal
samplers with those same independent targets. Compare retained Laplace outputs
with their pinned original source, preserving its approximation and precision;
those outputs are not an exact-posterior oracle. Current native tests cover GC protection, numerical-error
cleanup, RNG restoration, serial/parallel equivalence and all outcomes.

The manually dispatched workflows named **Historical** audit pinned older
release artifacts. Current per-push Valgrind and sanitizer workflows build and
test the checked-out source. Neither historical evidence nor a passing toy
reference certifies a new real-data posterior analysis.


For strict serial/parallel replay, start the parent R process with the same
one-thread BLAS/OpenMP environment used by workers. The CI workflow-level
settings do this before R loads its math libraries. The unit tests compare
worker computations under matching conditions; `verify-marginal-workers.R`
also compares public serial and parallel fits in the controlled parent.
Do not replace exact assertions with a tolerance to hide runtime differences.
