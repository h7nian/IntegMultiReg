# Migrating to IntegMultiReg 0.2.0

Version 0.2.0 deliberately removes the old aliases. Update calls as follows:

| 0.1.x | 0.2.0 |
|---|---|
| `platform_data_list` | generic argument `x` |
| `cov` | `covariates` |
| `type_outcome` | `outcome_type` |
| `ssize` | `min_subgroup_size` |
| `sample_mcmc = c(draws, burnin)` | separate `draws` and `burnin` |
| `hh` | `molecular_prior_scale` |
| `h0` | `forced_prior_scale` |
| `sig_alpha_psi` | `residual_prior = c(shape = ..., rate = ...)` |
| `thet_alph_bet` | `interaction_prior = c(shape = ..., rate = ...)` |
| `method = "IMR"` / `"BMS"` | `method = "imr"` / `"bms"` |
| prediction column `predict` | `prediction` |
| `predict_imr()` | `predict()` |
| `summary_imr()` | `summary()` |
| `coef_imr()` | `coef()` |
| `plot_imr()` | `plot()` |

`cv_method = "postfit_original"` is the default. Use `cv_method = "refit"` to keep
the 0.1.4 training-fold refitting workflow, or `cv_method = "importance"`
for the paper-derived empirical importance average. See `?cv_imr` for scope,
scoring and partition differences. `validation` records the selected method.

Cross-validation results now expose ordinary named fields:

```r
cv <- cv_imr(fit)
cv$pooled
cv$fold_mean
cv$predictions
cv$metric
cv$validation
cv$control
```

Saved 0.1.x fits are not upgraded implicitly. Convert a complete legacy object
once and save the result:

```r
fit_v2 <- upgrade_imr_fit(fit_v1)
saveRDS(fit_v2, "fit-v2.rds")
```

If `upgrade_imr_fit()` reports missing or inconsistent state, refit the model.
It intentionally does not guess values or repair corrupted objects.

Optional paper/released-code choices are documented in `METHOD-COVERAGE.md`.
Old fits lack a saved native random state and must be refitted before using
`fold_rng="continue"`. Missing sampler/numerical metadata retain the historical
package defaults; upgrading an object never changes its posterior draws.

## Convention names in development snapshots

Earlier 0.2.0 snapshots used labels tied to the article and implementation
history. Update explicit calls with this one-to-one mapping:

| Earlier value | Current value | Behavior |
|---|---|---|
| `sampler_method = "paper"` | `"corrected"` | Symmetric MRF update and boundary Hastings correction |
| `sampler_method = "legacy"` | `"original"` | Original selection updates, for historical reproduction |
| `cv_method = "legacy"` | `"postfit_original"` | Original post-fit defaults for model states, df and scores |
| `prior_indexing = "code2017"` | `"original"` | Original precision-block boundary |
| `df_method = "legacy_integer"` | `"integer"` | Integer-truncated predictive degrees of freedom |
| `score_method = "legacy"` | `"original"` | Original AUC and concordance rules |

`refit`, `importance`, `standard` and `fractional` retain their meanings.
The retired argument values are not accepted as aliases. To read a saved fit
with the earlier labels, call `upgrade_imr_fit()` explicitly. It changes the
metadata labels, preserves the original call and all posterior draws, and
validates the result. It does not rerun the sampler or convert an original
fit into a corrected one. Saved fold tables can still be supplied through
`folds`; archive the original CV result with its source version.
