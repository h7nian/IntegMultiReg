# Methods and reproducibility

`IntegMultiReg` separates the model variant from its numerical calculation.
`model_variant = "imr"` couples subgroup selection indicators through an MRF;
`"bms"` fits the availability subgroups independently. Both use the same outcome
families and coefficient priors.

## Fitting and uncertainty

Ordinary fitting uses the stated symmetric-MRF conditional log-odds,
`nu + 2 * sum(theta * gamma)`, and the Hastings correction for flip proposals at
empty/full selection boundaries. Its Gamma log-prior score has the negative
rate term and includes the log term for every positive interaction. Prior
precision follows coefficient roles: intercept and clinical effects use the
forced scale; molecular effects use the molecular scale.

The fitting engine integrates coefficients and residual variance using its
Laplace-based model scores. `inclusion_probabilities()` extracts retained
selection frequencies. `selection_summary()` summarizes selection indicators
and MRF interactions. Neither returns regression effects.

`sample_regression_posterior()` performs additional conditional pMOM sampling.
Its model weights inherit the fitted selection chain and its Laplace
approximation. `output_draws` controls returned samples;
`min_draws_per_model_chain` controls the minimum conditional-chain length.
The recorded `draws_per_model_chain` is the actual length. Conditional split
R-hat does not diagnose the original selection chain or remove its uncertainty.

## Predictive validation

| Choice | Computation | Scope |
|---|---|---|
| `cv_method = "refit"` (default) | Refit preprocessing, formula encoding, selection MCMC and prediction in every training fold | Evaluate the fitting procedure with fixed hyperparameters |
| `cv_method = "reweight"` | Reuse the full-data fit and apply inverse-density weights to its selection states | Post-fit approximation conditional on full-fit preprocessing and latent-response means |
| `model_set = "all_draws"` | Retain every sampled state and its multiplicity | Default state collection for reweighting |
| `model_set = "top_unique"` | Keep at most `max_models` ranked distinct states | A truncated collection for the same reweighting engine |

Held-out outcomes enter reweighting as an importance correction; their use alone
does not identify an error. The full-fit augmented-response mean plug-in is
stated in Section 4.1 of Chekouo et al. (2017). It does not make this procedure
an independent training-fold refit.

For reweighting, `ridge = 0.001` is the default coefficient penalty, including
the intercept. The article's unpenalized estimate requires `ridge = 0` explicitly
and a full-column-rank training design for every state. Singular systems stop
with fold/subgroup context. A positive ridge can stabilize them but changes the
estimator. `df_method = "fractional"` retains numeric predictive degrees of
freedom; `"integer"` truncates them. Scores use the same standard pair and tie
rules in both public CV algorithms.

Actual folds, summation order and refit seeds are recorded for replay. Refit uses
R's RNG and reweighting uses GSL; identical seeds across those generators do not
imply identical partitions. Supplied folds match subject IDs. Tuning or formula
selection requires an outer validation layer, illustrated in the
[covariate comparison guide](https://h7nian.github.io/IntegMultiReg/articles/covariate-comparison.html).

## Historical implementations

The released selection update used a neighbour coefficient of one and omitted
boundary Hastings ratios. These changes affect its stationary distribution.
Its diagnostic log score also used a positive Gamma rate term and omitted log
terms below 0.001. The theta acceptance calculation used the negative rate sign;
the score error should not be described as a reversed Gamma sampling prior.
An additional precision-index boundary assigned molecular precision to the last
forced coefficient in the released calculation.

These historical conventions are not alternative choices in the ordinary fitting
interface. Retain the exact archived source and its scripts for replay.
`upgrade_imr_object()` converts stored structure and labels without recomputing
draws. It does not turn historical draws into samples from the current target.
The [migration guide](https://h7nian.github.io/IntegMultiReg/migration.html) identifies renamed quantities and controls.

## What validation establishes

Independent tests check the MRF increments and proposal ratios against enumerated
conditional targets, Gamma increments against log densities, and conditional CV
against independently compiled archived code. Other tests cover fold replay,
worker equality, numerical failures and object boundaries. Source-specific
manifests and logs identify which build and environment each result belongs to.

Those checks do not establish reproduction of an entire study or convergence of
its long chains. Historical Table 1 and Figure 3 additionally depend on data,
preprocessing, starting states, seeds, replicate design and aggregation. The
supplement contains 778 gene columns whereas the article reports 776; an analysis
must identify its data source. Historical starting values and replicate seeds
are not fully recorded. A new reference table does not retroactively explain an
old unexplained discrepancy.

The [replay guide](https://h7nian.github.io/IntegMultiReg/articles/reproducibility.html) describes a record connecting
source identity, data, parameters, actual seeds/folds, generating code and
unrounded results. Computational repeatability and scientific validation remain
different claims.
