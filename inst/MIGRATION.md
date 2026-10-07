# Migrating to IntegMultiReg 0.2.0

The development interface names the model, computation and returned quantities
separately. Retired argument and function names are not accepted as aliases.
Use the source snapshot recorded with an analysis when replaying old results.

## Data and fitting

| Earlier interface | Current interface |
|---|---|
| `platform_data_list` | generic argument `x` |
| `cov` | `covariates` |
| `type_outcome` | `outcome_type` |
| `ssize` | `min_subgroup_size` |
| `sample_mcmc = c(draws, burnin)` | separate `draws` and `burnin` |
| `hh` / `h0` | `molecular_prior_scale` / `forced_prior_scale` |
| `sig_alpha_psi` | `residual_prior = c(shape = ..., rate = ...)` |
| `thet_alph_bet` | `interaction_prior = c(shape = ..., rate = ...)` |
| `method = "IMR"`, `"BMS"`, `"imr"` or `"bms"` | `model_variant = "imr"` or `"bms"` |
| `sampler_method = "paper"` or `"corrected"` | omit the argument; fitting always uses the symmetric MRF and Hastings correction |
| `sampler_method = "legacy"` or `"original"` | replay with the archived source; a current fit has a different target |
| `prior_indexing = "standard"` | omit the argument; prior precision follows coefficient roles |
| `prior_indexing = "code2017"` or `"original"` | replay with the archived source |

These changes include a default-algorithm correction, not just renamed labels.
Ordinary fitting no longer offers the known historical departures from the
stated model. The priors and outcome definitions remain unchanged, including
the concentrated inverse-Gamma binary variance prior rather than a fixed value
of one. Old draws cannot be corrected by changing their metadata.

## Returned quantities and additional computation

| Earlier interface | Current interface |
|---|---|
| `coef(fit)` or `coef_imr(fit)` for inclusion probabilities | `inclusion_probabilities(fit)` |
| `posterior_draws(fit)` | `sample_regression_posterior(fit)` |
| `draws` in that additional sampler | `output_draws` |
| `conditional_draws` in that sampler | `min_draws_per_model_chain` |
| `posterior_summary(fit)` | `selection_summary(fit)` |
| `compare_imr(...)` | `compare_fit_summaries(...)` |
| `validate_imr(fit)` | `validate_imr_object(fit)` |
| `upgrade_imr_fit(object)` | `upgrade_imr_object(object)` |
| `type = "mean"` for regression-posterior prediction | `quantity = "conditional_mean"` |
| `type = "response"` for regression-posterior prediction | `quantity = "new_observation"` |
| prediction-list key `model:011` | `subgroup:011` |
| prediction column `predict` | `prediction` |
| `predict_imr()`, `summary_imr()`, `plot_imr()` | `predict()`, `summary()`, `plot()` |

`coef(fit)` now stops with an explanation: the selection fit does not store
regression coefficients. It neither returns probabilities under a coefficient
name nor starts MCMC implicitly. Call `sample_regression_posterior()` explicitly;
`coef()` on its result returns coefficient posterior means.

`confint(fit)` requires `parm = "selection"`, `"theta"` or `"all"`. These are
intervals for selection indicators and interactions. Regression coefficient,
variance and latent-response intervals come from the `imr_posterior` object,
using `parm = "coefficients"`, `"variance"` or `"latent"`.

The extra sampler's control record uses `output_draws` and
`min_draws_per_model_chain`. Its diagnostics use `selection_model`,
`output_draws` and `draws_per_model_chain`; the last is the actual chain length,
which can exceed the minimum. `selection_draw_index` identifies source draws.
Prediction tables retain full numeric precision and now print the availability
subgroup's measured platforms.

## Cross-validation

`cv_method = "refit"` is the default. This rebuilds preprocessing, formula
encoding and MCMC inside each training fold, at roughly `k * rounds` fits.
Choose `"reweight"` explicitly to reuse a full-data fit. Its state collection is
an independent choice:

```r
cv_imr(fit, cv_method = "reweight", model_set = "all_draws")
cv_imr(fit, cv_method = "reweight", model_set = "top_unique", max_models = 100)
```

The old `"importance"` calculation corresponds to reweighting with all draws.
The `model_set` names `"draws"` and `"ranked_unique"` become `"all_draws"` and
`"top_unique"`. `df_method = "legacy_integer"` becomes `"integer"`.
Both CV algorithms now use the standard metric definitions. Historical
`legacy`/`postfit_original` presets and `score_method` are not ordinary API
options; use their archived source for exact replay. Substituting a new state
collection alone does not reproduce an old bundle of scoring and df rules.
The redundant fitted-`method` assertion was removed from `cv_imr()`; fit a
separate BMS object when that model variant is wanted.

Results retain `pooled`, `fold_mean`, `predictions`, `metric`, `validation` and
`control`. The control record states `model_variant`, `selection_update` and
`score_rule`, as well as the effective numerical settings and actual folds.

## Saved objects

New fits use schema version 3. Upgrade earlier fits or regression-draw objects
explicitly before using them:

```r
updated <- upgrade_imr_object(saved_object)
validate_imr_object(updated)
saveRDS(updated, "updated-object.rds")
```

The upgrade preserves posterior draws and the original call. It does not repair
incomplete or corrupted objects, rerun sampling or replace a historical target
with the stated posterior. Fits made with historical update or precision rules
remain inspectable, but new prediction, CV and conditional sampling require
refitting with the current interface. See the
[methods guide](https://h7nian.github.io/IntegMultiReg/method-coverage.html) for the statistical distinctions.
