# Migrating to version 0.3.0

Version 0.3.0 groups fitting controls and offers four sampling choices. The
original Laplace selection algorithm remains available. New joint and marginal
samplers use the same pMOM regression and MRF selection model without Laplace
model scores. Preserve earlier source archives, scripts and references when
replaying historical results; fitted-object renaming cannot change a posterior
target or create missing parameter draws.

## Fitting controls

```r
fit <- imr(
  analysis,
  marginalize = "none",
  priors = imr_priors(nu = c(-4, -3, -4)),
  mcmc = imr_mcmc(draws = 4000, burnin = 2000, chains = 4, seed = 1),
  control = imr_control()
)
coef(fit)
confint(fit)
```

`draws` is retained iterations per chain. Total transitions per chain are
`burnin + draws * thin`. Initialization and sampling use separate recorded
seeds. `workers` inside `imr_mcmc()` distributes chains; `cv_imr(workers = ...)`
distributes folds and uses serial chains within each fold.

| Earlier use | Current use |
|---|---|
| Separate prior arguments in `imr()` | `priors = imr_priors(nu, molecular_scale, forced_scale, residual, interaction)` |
| Separate sampling arguments in `imr()` | `mcmc = imr_mcmc(draws, burnin, chains, seed, initial)` |
| Laplace selection fit | `marginalize = "coefficients_and_variance"` |
| `laplace_max_iter`, `laplace_tolerance` in `imr()` | `control = imr_control(...)`; original numerical defaults retained |
| `max_models` in point prediction | Still supported by `predict()` for a Laplace selection fit |
| `sample_regression_posterior(fit, ...)` | Still performs optional conditional sampling for a Laplace selection fit |
| `coef(fit)` unavailable on selection-only fits | Still unavailable there; returns coefficient posterior means for the other three choices |
| `selection_summary(fit)` | `summary(fit, parm = "selection")` or `parm = "interaction"` |
| `plot(fit, type = "theta")` / `"theta_trace"` | `type = "interaction"` / `"interaction_trace"` |
| `upgrade_imr_object(old_fit)` | Refit saved inputs using the required sampler and settings |
| CV with `ridge`, `model_set`, `df_method` or `fold_rng` | Full refits, or PSIS reweighting when regression draws are stored |

Obsolete arguments fail rather than act as aliases. Binary latent variance
retains the anchored inverse-gamma prior; a conflicting explicit residual prior
is rejected. The four choices are described in the
[sampler guide](https://h7nian.github.io/IntegMultiReg/articles/marginalization.html).

## Retaining the original calculation

For a one-chain replay of the previously corrected Laplace calculation, specify
`marginalize = "coefficients_and_variance"` and
`mcmc = imr_mcmc(chains = 1, initial = NULL, thin = 1, ...)`, with the original
seed, burn-in, retained count, data and priors. `initial = NULL` requests the
original native initialization. The default dispersed starts and four chains
will produce different draws even though the transition algorithm is retained.
Use the recorded toolchain and source for cross-version numerical comparisons.

A Laplace fit has class `c("imr_selection", "imr")`. Its `predict()` method
retains rescored-model point prediction: binary probabilities and continuous
or survival working-scale predictions. For a log-time model, this means log
predicted time. `sample_regression_posterior()` adds conditional regression
sampling with empirical selection weights. Prediction from that returned
`imr_posterior` object uses the response/link convention of regression draws;
log-time response-scale summaries are medians of transformed draws.

The three other choices store coefficient and variance draws directly,
recovering a marginalized parameter from its conditional distribution.
`posterior_draws()` extracts these arrays without additional MCMC.

## Summaries, diagnostics and storage

All summary and interval quantiles now use the empirical inverse CDF (type 1),
with a floating-point boundary tolerance. Earlier conditional summaries used
type 7 and selection summaries used type 8, so interval endpoints may differ
even for the same draws. Type 1 preserves binary support and exclusion zeros.

Rank-normalized split/folded R-hat and bulk/tail ESS assess fitted chains.
Conditional-regression diagnostics use their separate fixed-model chains;
returned mixture rows are not extra chains. Inspect the originating selection
fit separately. Constant parameters remain unassessed.

Fits with regression draws derive molecular selection as `beta != 0`.
Laplace fits instead store integer selection histories without coefficient
arrays. `max_draw_memory_mb` limits estimated retained-array storage, not total
process memory. Fitting, parallel workers and summaries need additional space.
Thinning reduces storage but does not improve mixing or eliminate intervening
updates.

## Earlier names

The CRAN 0.1.3 release and earlier 0.2.0 snapshots used different names. Current
names include `model_variant = "imr"/"bms"`, `inclusion_probabilities()`,
`compare_fit_summaries()` and `validate_imr_object()`. Historical `paper`,
`legacy`, `original`, `corrected` and precision/scoring selectors belong to
their archived versions. Current fitting uses the symmetric MRF, the boundary
Hastings correction where required, and coefficient-role prior indexing.


## Consistent extraction, prediction and conditional controls

`posterior_draws(fit)` now returns the same named family list for every sampler:
coefficients, variance, selection, interaction and latent. Unavailable families
are `NULL`. Previously the default selected coefficients for some fits and
selection indicators for others. Specify `parm = "coefficients"` or another
family when the caller needs its arrays directly.

Prediction methods share the same formals and default to `type = "response"`
and `interval = FALSE`. The `prediction` attribute records the calculation,
quantity, scale and summary. Laplace selection fits provide `model_average`
point predictions. For log-time survival the default exponentiates the original
working-scale point; `type = "link"` preserves its old numeric value. Binary
link predictions average the working predictors, while response predictions
retain the original weighted probabilities. Obtain conditional regression draws
before requesting posterior intervals or future observations from a selection
fit. No additional sampler runs implicitly during prediction.

Conditional sampling now uses the same settings constructor:

```r
posterior <- sample_regression_posterior(fit, output_draws = 1000,
  mcmc = imr_mcmc(draws = 200, burnin = 1000, chains = 2,
    seed = 1, keep_latent = TRUE))
summary(posterior)$parameters
predict(posterior, new_data, interval = TRUE)
```

Move the former `min_draws_per_model_chain` to `mcmc$draws`, and `latent` to
`mcmc$keep_latent`; `output_draws` retains its mixture-size meaning. Here the
chain budget is a minimum per fixed model. `thin` changes storage, `workers`
applies to large diagnostic blocks, and the conditional updates retain their
serial random-number sequence. Old flat arguments are not aliases.
Conditional summaries now use the same `summary.imr` container as fitting
summaries, with `draw_type = "conditional mixture"` and no fictitious chain count.

Non-default settings that do not apply to the selected algorithm now fail at
fitting. The supplied control objects remain in the fit for replay; the active
controls are listed separately in `control$effective`. `mcmc_diagnostics(fit,
type = "sampler")` returns a common table of update methods, proposal counts,
accepted counts and rates. Inapplicable rates remain `NA`. Laplace counters are
exported from the original computation without changing its proposals.
