# Migrating to joint inference

Version 0.3.0 fits selection and regression parameters in one joint MCMC
process. It changes the inference engine and fitted-object schema, not only
function names. Preserve earlier source archives, scripts and references when
replaying historical results. Existing selection-only fits lack joint
coefficient samples and must be refitted for the new methods.

## One fitting call

Priors and MCMC settings have their own validated constructors:

```r
fit <- imr(
  analysis,
  priors = imr_priors(nu = c(-4, -3, -4)),
  mcmc = imr_mcmc(draws = 4000, burnin = 2000, chains = 4, seed = 1)
)
coef(fit)
confint(fit)
```

`draws` is retained iterations per chain. Total transitions per chain are
`burnin + draws * thin`. Chain initialization and sampling use separate recorded
seeds. `workers` inside `imr_mcmc()` distributes chains; `cv_imr(workers = ...)`
distributes folds and uses serial chains within each fold.

| Earlier use | Current use |
|---|---|
| `imr(..., nu, molecular_prior_scale, forced_prior_scale, residual_prior, interaction_prior)` | `priors = imr_priors(nu, molecular_scale, forced_scale, residual, interaction)` |
| `imr(..., draws, burnin, seed, initial)` | `mcmc = imr_mcmc(draws, burnin, chains, seed, initial)` |
| `sample_regression_posterior(fit, ...)` | Fit once with the required MCMC budget; `posterior_draws(fit)` extracts samples |
| `coef(fit)` unavailable on selection-only fits | `coef(fit)` returns regression posterior means |
| `selection_summary(fit)` | `summary(fit, parm = "selection")` and `summary(fit, parm = "interaction")` |
| `confint(fit)` required a selection/interaction family | Defaults to regression coefficients; `parm` selects other families or coefficient names |
| `plot(fit, type = "theta")` / `"theta_trace"` | `type = "interaction"` / `"interaction_trace"` |
| `imr_posterior` methods | Standard methods on the `imr` fit itself |
| `upgrade_imr_object(old_fit)` | Refit saved inputs; an object conversion cannot create missing joint samples |
| `laplace_max_iter`, `laplace_tolerance`, ranked-model prediction limits | Removed with the Laplace and coefficient-mode engines |
| Reweighting with `ridge`, `model_set`, `df_method` or `fold_rng` | PSIS on joint observed-data likelihood ratios, with per-fold diagnostics |

Obsolete arguments are rejected rather than treated as aliases. Binary latent
variance retains its anchored inverse-gamma prior; a conflicting explicit
residual prior is now rejected rather than silently overwritten.

## Interpretation and storage

Coefficient intervals and predictions use the same joint draws as selection.
`posterior_draws()` returns iteration-by-chain-by-parameter arrays. Molecular
selection is recoverable as `beta != 0`, so the fit stores no duplicate selection
history or inclusion-probability cache. Use `inclusion_probabilities()` for
per-platform probability matrices.

Summary and interval quantiles use the empirical inverse CDF (type 1), with a
floating-point boundary tolerance. This preserves the support of binary
indicators and point masses at zero. Rank-normalized split/folded R-hat and
bulk/tail ESS replace the separate conditional-chain diagnostic. Constant
parameters remain unassessed; this is not a claim of convergence.

The joint arrays use more storage than selection-only draws. Inspect the
reported retained-array estimate; `max_draw_memory_mb` is an explicit retained
storage budget, not total process memory. Fitting, parallel workers and summary
operations need additional space. Thinning reduces storage but does not improve
mixing or reduce the required intervening updates.

## Earlier API names

The CRAN 0.1.3 release and earlier 0.2.0 snapshots used different argument and
function names. The maintained interface uses `model_variant = "imr"/"bms"`,
`inclusion_probabilities()`, `compare_fit_summaries()` and
`validate_imr_object()`. Historical `paper`, `legacy`, `original`, `corrected`
and precision/scoring selectors belong to their archived versions. They are
not choices in the current joint sampler.
