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

Cross-validation results now expose ordinary named fields:

```r
cv <- cv_imr(fit)
cv$pooled
cv$fold_mean
cv$predictions
cv$metric
cv$validation
```

Saved 0.1.x fits are not upgraded implicitly. Convert a complete legacy object
once and save the result:

```r
fit_v2 <- upgrade_imr_fit(fit_v1)
saveRDS(fit_v2, "fit-v2.rds")
```

If `upgrade_imr_fit()` reports missing or inconsistent state, refit the model.
It intentionally does not guess values or repair corrupted objects.
