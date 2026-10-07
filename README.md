# IntegMultiReg

[![R-CMD-check](https://github.com/h7nian/IntegMultiReg/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/h7nian/IntegMultiReg/actions/workflows/R-CMD-check.yaml)
[![Package website](https://github.com/h7nian/IntegMultiReg/actions/workflows/pkgdown.yaml/badge.svg)](https://h7nian.github.io/IntegMultiReg/)

**Integrative Bayesian Multiple Regression for Multi-Platform Biomarkers.**

`IntegMultiReg` implements the integrative multi-regression (IMR) model of
Chekouo, Stingo, Doecke and Do (2017, *Biometrics*) and extends it from
time-to-event outcomes to continuous (Gaussian) and binary (probit) outcomes.

Given several molecular platforms measured on overlapping but partially missing
sets of subjects, IMR partitions subjects into the availability subgroups of a
Venn diagram, fits one regression per subgroup, and shares information across
availability subgroups through a **Markov random field (MRF) prior** on the
variable-selection indicators. Each regression uses **non-local product-moment
priors** on its active coefficients. Subjects can therefore contribute when
some platforms are missing.

Read the [getting started guide](https://h7nian.github.io/IntegMultiReg/articles/IntegMultiReg.html)
or browse the [function reference](https://h7nian.github.io/IntegMultiReg/reference/index.html).
These pages describe the GitHub development version, 0.2.0. The
[reproducibility guide](https://h7nian.github.io/IntegMultiReg/articles/reproducibility.html)
explains how to retain the source, settings, folds and results for an analysis.

## Installation

The package contains C code that links against the
[GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/), which must be
installed first:

* macOS: `brew install gsl`
* Debian/Ubuntu: `sudo apt-get install libgsl-dev`
* Fedora/RHEL: `sudo dnf install gsl-devel`
* Windows: GSL is provided by Rtools.

To install the development interface from GitHub:

```r
install.packages("remotes")
remotes::install_github("h7nian/IntegMultiReg")
```

For numerical replication, use the source snapshot and checksum recorded with
the analysis. The [replay guide](https://h7nian.github.io/IntegMultiReg/articles/reproducibility.html)
explains which inputs and settings to retain.

The CRAN release, 0.1.3, uses the earlier API. The 0.2.0 examples on this
site require the GitHub version. See the
[migration guide](https://h7nian.github.io/IntegMultiReg/migration.html)
when updating existing code. To install the CRAN release:

```r
install.packages("IntegMultiReg")
```

Alternatively, install a local source tarball:

```r
install.packages("IntegMultiReg_0.2.0.tar.gz", repos = NULL, type = "source")
```

## Quick start

Fitting uses the symmetric MRF conditional and Hastings-adjusted selection
moves. `model_variant = "imr"` shares selection information across subgroups;
`"bms"` fits them independently. CV refits each training fold by default.

```r
library(IntegMultiReg)
data("simIMR")

analysis <- imr_data(
  platforms = simIMR$platforms,
  outcome = simIMR$outcome.binary,
  covariates = simIMR$covariates,
  outcome_type = "binary"
)
analysis
stopifnot(validate_imr_data(analysis))

fit <- imr(
  analysis, nu = c(-4, -3, -4),
  draws = 2000, burnin = 1000,
  min_subgroup_size = 30,
  seed = 1
)

fit                       # short summary
summary(fit)              # selected biomarkers per platform
inclusion_probabilities(fit)                 # per-platform mPIP matrices
plot(fit, type = "selection")
plot_top_features(fit)    # ranked biomarker bar chart
new_data <- imr_data(simIMR$platforms[1:2], covariates = simIMR$covariates)
predict(fit, newdata = new_data)
cv_imr(fit, k = 3, rounds = 1, cv_method = "refit")
```

## Real-data example

`kircIMR` is a reduced public UCSC Xena TCGA-KIRC survival example aligned with
the Biometrics kidney cancer case study: mRNA expression, miRNA expression, DNA
methylation, clinical covariates and right-censored survival.  It is derived
from public UCSC Xena TCGA-KIRC sampleMap files, not from controlled-access
TCGA/GDC files, and contains only a reduced Cox-screened feature panel.

The package replaces TCGA barcodes with package-internal IDs such as `KIRC001`
and does not distribute a barcode mapping.  Users should not attempt
participant re-identification or linkage to external resources.

```r
data("kircIMR")
sapply(kircIMR$platforms, dim)
kircIMR$model_subgroup_sizes

kirc_fit <- imr(
  kircIMR$platforms,
  kircIMR$outcome.survival,
  covariates = kircIMR$covariates,
  outcome_type = "right.censored",
  nu = c(-4, -3, -4),
  draws = 4000, burnin = 1000,
  min_subgroup_size = 30,
  seed = 1
)
```

See the package vignette `vignette("IntegMultiReg")` for a complete walk-through.

## Reference

Chekouo T, Stingo FC, Doecke JD, Do K-A (2017). "A Bayesian Integrative
Approach for Multi-Platform Genomic Data: A Kidney Cancer Case Study."
*Biometrics*, **73**(2), 615–624. <https://doi.org/10.1111/biom.12587>

[Read paper](https://academic.oup.com/biometrics/article/73/2/615/7537638) ·
[Publisher PDF](https://academic.oup.com/biometrics/article-pdf/73/2/615/55973435/biometrics_73_2_615.pdf)

When using `kircIMR`, please also acknowledge TCGA, the National Cancer
Institute Genomic Data Commons, and UCSC Xena as the public data sources.

## Coefficient and predictive uncertainty

`coef(fit)` does not return inclusion probabilities or silently run another
sampler. Use `inclusion_probabilities(fit)` for variable-selection probabilities.
After fitting, `sample_regression_posterior(fit)` explicitly adds conditional pMOM coefficient and
variance draws to the retained selection models. `summary(draws)` and
`confint(draws)` report coefficient intervals including point mass at zero for
inactive molecular features. Clinical covariates remain always included.

```r
# Use an adequately explored fit and inspect both stages of diagnostics.
draws <- sample_regression_posterior(fit, seed = 2)
draws$diagnostics
confint(draws)
predict(draws, simIMR$platforms, covariates = simIMR$covariates,
        quantity = "conditional_mean")       # uncertainty in the conditional response mean
predict(draws, simIMR$platforms, covariates = simIMR$covariates,
        quantity = "new_observation")   # uncertainty in a future outcome
```

These are approximate model-averaged intervals: selection weights retain the
original Laplace approximation. Conditional split R-hat does not assess the
original selection chain; increase simulation effort when diagnostics are poor.
Coefficients use subgroup-standardized predictor scales. Binary probability
intervals use `quantity = "conditional_mean"`; binary new-observation draws are zero/one, with interval endpoints obtained by interpolated empirical quantiles.

## Survival migration from 0.1.2

Version 0.1.3 logs positive event and censoring times once, matching the original
Biometrics supplementary code. Supply raw times, including positive times below
one; refit previous survival models. The optional `survival_scale = "identity"`
reproduces the historical raw-time implementation. CV partitioning is unchanged.

For log-time fits, `predict(fit, ...)` returns a log-time point prediction.
`predict(draws, ...)` returns time-scale intervals and median point summaries.
The time-scale posterior mean need not exist under the variance mixture; it is
not estimated by averaging exponentiated draws. With `quantity = "conditional_mean"` the interval
summarizes conditional mean time, whereas `quantity = "new_observation"` includes future
outcome variability. The censoring process for future observations is not modeled.

## Inspect numerical computation

`fit$control$laplace_diagnostics` records fitting-stage calls, iteration-limit
hits and numerical failures by subgroup. Review these with the selection-chain
diagnostics. The counters neither establish MCMC convergence nor cover later
prediction-stage optimization.

## Cross-validation algorithms

`cv_imr(fit, cv_method = "refit")` is the default. It rebuilds preprocessing,
selection MCMC and prediction inside every training fold, at approximately
`k * rounds` full fits.

`cv_method = "reweight"` reuses the full-data fit with inverse-density
importance weights. The separate `model_set` argument chooses the state collection:

```r
cv_imr(fit, cv_method = "reweight", model_set = "all_draws")
cv_imr(fit, cv_method = "reweight", model_set = "top_unique", max_models = 100)
```

`all_draws` preserves every retained draw and its multiplicity. `top_unique`
uses at most `max_models` ranked distinct selection models. Reweighting conditions
on full-fit preprocessing and augmented-response means. Held-out responses enter
its importance correction; their presence alone does not show an algorithm error.
This approximation does not refit the training folds. Its default ridge penalty
is 0.001; `ridge = 0` requires every training design to have full column rank.
`df_method = "fractional"` keeps the numeric predictive degrees of freedom;
`"integer"` truncates them. Both CV algorithms use the same standard scores.

The [methods guide](https://h7nian.github.io/IntegMultiReg/method-coverage.html)
describes these approximations and separates them from archived computations.
Prior scales, `laplace_max_iter` and `laplace_tolerance` are recorded in the fit.

Use `initial=list(selection=..., interaction=...)` for different chain starts.
Selection matrices follow `inclusion_probabilities(fit)` and contain zero/one
entries. Interaction matrices are symmetric, with positive off-diagonals and
zero diagonal. Refit CV reuses specified starts. Different starts and longer
chains support convergence assessment; they do not establish convergence by
themselves.

`cv$control` records effective settings and actual folds. Refit CV also records
`refit_seeds`. Use `folds = cv$control$folds` to replay partitions and fitting
seeds within the same runtime and RNG kind. Matching seed numbers across the
R and GSL generators does not imply identical partitions. Data-driven tuning
requires an outer validation layer.

Both methods accept `workers = 2L` (or another positive integer) for
PSOCK process parallelism. The default `workers = 1L` remains serial. Each
algorithm preserves its own partitions, seeds, prediction order and scoring;
changing the worker count does not select a different validation algorithm.
Process startup may outweigh the benefit for short runs, and each worker needs
its own fit/workspace memory. Avoid nesting CV workers inside parallel experiment
runs. For custom formula functions, use a serializable local formula environment
or a package-qualified function name; the global workspace is not exported.


### Compare prespecified covariate formulas

A runnable example compares two formulas on identical subject/fold assignments,
then uses nested cross-validation to evaluate formula selection using only
inner training data. Clinical terms remain forced within each candidate model.
Required adjustment terms must appear in every candidate; the example is not
causal confounder selection.

```r
source(system.file("examples", "compare-covariates.R", package = "IntegMultiReg"))
comparison <- run_covariate_comparison()
comparison$paired_summary
comparison$nested_summary
```

Follow the [worked guide](https://h7nian.github.io/IntegMultiReg/articles/covariate-comparison.html)
for paired folds, inner selection and outer evaluation. The full example uses
synthetic data, writes fold and seed records, and is repeated in CI. Use
`quick = TRUE` only for a smoke run. Earlier covariate-example results generated
with historical updates remain in their versioned replication archive.

## Development

See [CONTRIBUTING](https://github.com/h7nian/IntegMultiReg/blob/main/CONTRIBUTING.md)
for source organization, documentation generation and local validation.
