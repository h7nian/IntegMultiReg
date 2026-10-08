# Methods and reproducibility

## Statistical model

The package retains the availability-subgroup regressions, symmetric MRF
selection prior and first-order pMOM coefficient priors. Clinical coefficients
and the intercept are always active. Molecular coefficients have an exclusion
point mass at zero. Residual variances use inverse-gamma shape/rate priors and
MRF interactions use Gamma shape/rate priors.

Continuous, binary and right-censored outcomes share this structure. Binary
utilities and censored working responses are updated under the chosen sampler's
joint or marginal target. Positive survival times are logged by default. Binary variance
uses the original concentrated prior near one; it is not a point mass.

## Sampling choices

`marginalize` names the regression parameters integrated out of the main-chain
target. `"none"` retains coefficients and residual variance; `"variance"`
retains coefficients; `"coefficients"` retains residual variance. Each supplies
regression draws, using conditional recovery for a marginalized parameter.
The coefficient integral uses positive Gaussian quadrature that is exact for
the pMOM polynomial, subject to floating-point error. Its node count grows as
`(d + 1)!`, so the current implementation is for small candidate models.

`"coefficients_and_variance"` retains the original Laplace calculation with
its corrected symmetric-MRF and boundary-Hastings updates and standard prior
indexing. The variance integral, approximation, proposal rules and numerical
defaults are preserved. It is not the newly derived exact full integral.
`sample_regression_posterior()` provides optional conditional regression draws
with the empirical model weights from this selection fit.

All four use the same model and priors. The retained Laplace path approximates
model weights, so its finite-run differences from the three exact-target paths
may include approximation error as well as Monte Carlo error. The
[sampler guide](https://h7nian.github.io/IntegMultiReg/articles/marginalization.html)
gives the targets and computational limits. All four sampling loops run in C,
with shared R preprocessing and result methods. Dense operations in the two
single-block marginal samplers retain R's BLAS/LAPACK interfaces. Frozen R
implementations are used only for exact-output regression tests.

## Joint posterior computation

With `marginalize = "none"`, version 0.3.0 samples selection, coefficients,
variances, interactions and latent responses together. The update of one feature integrates its coefficient in each
available subgroup, samples the full subgroup inclusion pattern, then samples
the active coefficients immediately. The scalar pMOM Bayes factor is analytic.
Whole-feature exchanges improve movement between correlated predictors.

Variance and augmented-response updates use full conditionals. Interaction
updates include the exact MRF normalizer and log-normal proposal correction.
The engine requires no Laplace model scores, coefficient modes or second-stage
regression MCMC. The
[joint posterior guide](https://h7nian.github.io/IntegMultiReg/articles/joint-posterior.html)
gives the formulas and validation scope.

`coef()`, `confint()`, `summary()` and `predict()` operate on the fitted joint
samples. `posterior_draws()` only extracts them. All chain identities are kept
for rank R-hat, ESS and MCSE from the posterior package. Undefined diagnostics
are retained; software completion and diagnostic cutoffs do not establish
scientific convergence or interval calibration.

## Predictive validation

| Choice | Computation | Scope |
|---|---|---|
| `cv_method = "refit"` | Refit preprocessing, formula encoding and the chosen sampler within each training fold | Evaluate the fitting procedure with fixed hyperparameters |
| `cv_method = "reweight"` | PSIS on inverse held-out observed-likelihood weights from the full joint posterior | Posterior approximation conditional on full-fit preprocessing |

Reweighting requires stored regression draws and is unavailable for a Laplace
selection fit. It uses the observed-data likelihood, integrating binary latent
utilities and using survival probabilities for censored observations. Held-out
outcomes belong in those importance ratios. Pareto k and effective sample size
are reported per fold; unstable weights require refitting. The procedure is
neither the old ranked-model calculation nor an OLS/ridge plug-in approximation.

Actual partitions and separate fold/chain seeds are saved. Matching folds pairs
candidate comparisons. The
[covariate guide](https://h7nian.github.io/IntegMultiReg/articles/covariate-comparison.html)
also separates inner formula selection from outer performance assessment.

## Archived calculations

The released 2017 selection update used a neighbour coefficient of one rather
than the symmetric MRF coefficient of two and omitted a boundary Hastings
correction. Its logged Gamma score also differed from a Gamma density, although
the interaction update itself used the negative rate sign. The released prior
precision indexing assigned the last clinical column to the molecular block.

Earlier package snapshots exposed those historical behaviors and later added
matching-model corrections. Version 0.2.0 still used Laplace-based selection and
optional conditional regression sampling. The new joint engine targets the
stated joint posterior directly and can produce different model probabilities
and predictions. Preserve exact archived source and scripts for numerical
comparisons; relabelling saved draws cannot change their sampling distribution.

Historical original-study replays, current-model examples and scientific
validation are separate evidence. Unknown original seeds/build settings, global
screening and finite-chain limitations remain relevant. A new expected CSV does
not explain a historical discrepancy. The
[replay guide](https://h7nian.github.io/IntegMultiReg/articles/reproducibility.html)
describes the source/data/settings/output record for each analysis.
