# IntegMultiReg

[![R-CMD-check](https://github.com/h7nian/IntegMultiReg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/h7nian/IntegMultiReg/actions/workflows/R-CMD-check.yaml)
[![Package website](https://github.com/h7nian/IntegMultiReg/actions/workflows/pkgdown.yaml/badge.svg)](https://h7nian.github.io/IntegMultiReg/)

Integrative Bayesian regression for molecular platforms measured on overlapping
sets of subjects. Each availability subgroup has its own regression; a Markov
random field shares feature-selection information across subgroups. Continuous,
binary and right-censored outcomes use non-local pMOM coefficient priors.

The development version, **0.3.0**, offers four sampling choices through
`imr(..., marginalize = ...)`. The default samples coefficients and residual
variance together. Two alternatives integrate out either parameter block, and
a fourth retains the original Laplace-based selection calculation. The first
three store regression draws for intervals and prediction. The [tutorial](https://h7nian.github.io/IntegMultiReg/articles/IntegMultiReg.html)
introduces the workflow, and the [reference](https://h7nian.github.io/IntegMultiReg/reference/index.html)
describes each function.

## Installation

```r
install.packages("remotes")
remotes::install_github("h7nian/IntegMultiReg")
```

Source installation needs the usual R compilation tools (Rtools on Windows)
and GSL for the retained Laplace engine. Install `libgsl-dev` on Debian/Ubuntu,
`gsl-devel` on Fedora/RHEL, or `gsl` with Homebrew on macOS. Windows source
builds use the GSL libraries supplied with current Rtools. The CRAN release, 0.1.3, has an older interface; these examples
require the development package. See the
[migration guide](https://h7nian.github.io/IntegMultiReg/migration.html).
For numerical replay, install the exact source snapshot recorded with the
analysis, as explained in the
[reproducibility guide](https://h7nian.github.io/IntegMultiReg/articles/reproducibility.html).

## Fit and inspect

```r
library(IntegMultiReg)
data("simIMR")

analysis <- imr_data(
  simIMR$platforms,
  outcome = simIMR$outcome.continuous,
  covariates = simIMR$covariates,
  outcome_type = "continuous"
)
validate_imr_data(analysis)

fit <- imr(
  analysis,
  priors = imr_priors(nu = c(-4, -3, -4)),
  mcmc = imr_mcmc(draws = 4000, burnin = 2000, chains = 4, seed = 1)
)
fit
coef(fit)                         # model-averaged regression effects
confint(fit)                      # coefficient credible intervals
inclusion_probabilities(fit)      # molecular selection probabilities
summary(fit, parm = "variance")
mcmc_diagnostics(fit)
```

`imr_priors()` contains statistical assumptions; `imr_mcmc()` contains sampling
and storage choices; `imr_control()` contains numerical settings. `draws` is the retained count **per chain**. Inspect the
actual diagnostics before interpreting an analysis. Constant indicators have
undefined R-hat and are reported separately. An illustration budget is not a
convergence guarantee.

## Choose which regression parameters to integrate out

| `marginalize` | Main-chain regression parameters | Computation |
|---|---|---|
| `"none"` (default) | Coefficients and residual variance | Joint pMOM sampling |
| `"variance"` | Coefficients | Analytic variance integral and weighted-t coefficient updates |
| `"coefficients"` | Residual variance | Polynomial-exact Gaussian integration; small models only |
| `"coefficients_and_variance"` | Neither | Original Laplace-based selection sampler |

The first three target the same model without Laplace model scores. The last
retains the original approximation, numerical tolerances and proposal rules.
Identical seeds need not give identical draws across different algorithms.
The [sampler guide](https://h7nian.github.io/IntegMultiReg/articles/marginalization.html)
explains parameter recovery, computational limits and comparisons.

For a Laplace fit, `predict()` retains the original point calculation and
`sample_regression_posterior()` runs optional conditional regression chains.
Their uncertainty retains the fitted selection weights and their approximation.

## Predict and validate

```r
new_data <- imr_data(simIMR$platforms[1:2], covariates = simIMR$covariates)
predict(fit, new_data, interval = TRUE)
predict(fit, new_data, quantity = "new_observation", interval = TRUE)

cv <- cv_imr(fit, k = 3, rounds = 1)
cv
cv$diagnostics
```

Prediction reuses the training transformations and all retained joint samples.
`conditional_mean` describes expected responses; `new_observation` includes
outcome noise. For log-time survival, response-scale summaries are medians;
`type = "link"` gives the working log-time scale.

CV refits the model and preprocessing within every training fold by default.
Reuse `cv$control$folds` to compare prespecified models on the same subjects and
partitions. The alternative `cv_method = "reweight"` applies PSIS to joint
observed-data likelihood ratios for fits that store regression draws and reports unstable weights. It retains
full-fit preprocessing and does not replace training-fold validation.

The [covariate comparison guide](https://h7nian.github.io/IntegMultiReg/articles/covariate-comparison.html)
shows formula/data fitting, paired CV and nested selection.

## Posterior samples and plots

```r
beta <- posterior_draws(fit)      # iteration x chain x coefficient arrays
plot(fit, type = "coefficient_trace", subgroup = 1, parameter = 1)
plot(fit, type = "selection")
plot_top_features(fit)
```

The fit stores exclusion zeros in molecular coefficient draws. Selection
histories and probabilities are derived from those zeros, avoiding duplicate
histories. Optional `imr_mcmc(keep_latent = TRUE)` retains augmented responses
for binary and censored models. `validate_imr_object()` checks object structure;
`compare_fit_summaries()` compares descriptive fit summaries.

## Data and methodology

`simIMR` is a synthetic example with known feature truth. `kircIMR` is a reduced
public UCSC Xena TCGA-KIRC illustration with mRNA, miRNA, methylation and clinical
covariates. Its globally screened feature panel is suitable for learning the
interface; it is not an unbiased predictive benchmark. The package uses
internal subject IDs and does not distribute a TCGA barcode mapping.

Chekouo T, Stingo FC, Doecke JD, Do K-A (2017). “A Bayesian Integrative Approach
for Multi-Platform Genomic Data: A Kidney Cancer Case Study.” *Biometrics*,
**73**(2), 615–624. [DOI](https://doi.org/10.1111/biom.12587) ·
[Publisher page](https://academic.oup.com/biometrics/article/73/2/615/7537638).

The [methods guide](https://h7nian.github.io/IntegMultiReg/method-coverage.html)
distinguishes the statistical model, the four sampling choices and archived
computational conventions. Historical results remain tied to their recorded source; the new
sampler is a substantive inference change and can change numerical results.
