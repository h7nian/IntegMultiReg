# Contributing to IntegMultiReg

## Local checking tools

The CRAN checking tools `checkbashisms` and `qpdf` are not runtime
dependencies. Package users do not need them. Maintainers running
`R CMD check --as-cran` locally can install them with
`brew install checkbashisms qpdf` on macOS or
`sudo apt-get install devscripts qpdf` on Debian/Ubuntu.

## Source layout and development checks


The public API uses `snake_case` names and S3 methods. Internal R functions use
`.imr_`; `%||%` is the internal null-default operator. Fitting options remain
explicit arguments, grouped by purpose in `?imr`. Their resolved values are
stored in the fitted object's `control` component.

| Responsibility | Source |
| --- | --- |
| Validated data and subject identifiers | `R/imr-data.R` |
| Fit dispatch and sampler inputs | `R/fit.R` |
| Formula response and design matrices | `R/formula.R` |
| Subject alignment and matrix scaling | `R/preprocessing.R` |
| Prediction and new-data routing | `R/predict.R` |
| CV settings, refit/post-fit algorithms and workers | `R/cv*.R` |
| Fit printing, summaries and feature rankings | `R/methods.R` |
| Fit validation, migration and posterior summaries | `R/diagnostics.R` |
| Plotting and shared appearance settings | `R/plots.R`, `R/plot-theme.R` |
| Conditional coefficient and latent-response sampling | `R/posterior-draws.R`, `R/posterior-kernel.R` |
| Explicit R-to-C argument conversion | `R/native-adapters.R` |
| Native sampler, prediction and CV | `src/` |

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
