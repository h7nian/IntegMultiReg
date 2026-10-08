# IntegMultiReg 0.3.0

* All samplers use one draw-family layout, prediction interface and sampler
  diagnostic table. Prediction results record their calculation and scale;
  Laplace survival predictions default to the response scale, with the original
  values available through `type = "link"`.
* Conditional computation uses `imr_mcmc()` and shares summary and interval
  handling with fitted draws. Non-default inapplicable controls fail explicitly.
  See the migration guide for changed extraction defaults and argument grouping.
* Removed inactive unadjusted selection updates, historical prior indexing and
  local-prior likelihood branches. The original ranked-score arithmetic remains
  explicit and separate from the posterior score.

* The variance- and coefficient-marginal samplers now run their updates,
  conditional recovery and sample storage in C. They share native MRF and
  storage helpers and preserve the R reference's arithmetic and random stream.
  Dense algebra continues through R's BLAS/LAPACK interfaces. The reference
  implementations are retained only in regression tests.
* `verbose = TRUE` now reports progress for both single-block marginal samplers.
  Reporting does not consume random numbers or change the sampled states.

* `imr()` exposes four `marginalize` choices: `"none"`, `"variance"`,
  `"coefficients"` and `"coefficients_and_variance"`. The original Laplace
  selection calculation remains available with its numerical defaults and
  proposals. The other paths target the same model through joint or marginal
  updates. Exact coefficient integration currently supports small models.
* `imr_priors()`, `imr_mcmc()` and `imr_control()` separate statistical,
  sampling and numerical settings. Four separately seeded chains are the
  default. Recorded starts, chain seeds and thinning define replay.
* Fits from the first three paths store regression draws. `coef()`, `confint()`,
  `summary()` and `predict()` use these draws; `posterior_draws()` only extracts
  them. Laplace fits retain their original point prediction and optional
  `sample_regression_posterior()` computation.
* Quantiles use an empirical inverse CDF (type 1), preserving point masses.
  Rank-normalized split/folded R-hat, bulk/tail ESS and MCSE assess chains.
  Conditional-model diagnostics remain separate from source selection chains.
* `cv_imr()` defaults to full training-fold refits with the selected sampler.
  Reweighting uses PSIS on joint observed-data likelihoods and requires stored
  regression draws. Earlier post-fit ridge and ranked-model CV are retired.
* A native fitted-mean cache and parallel diagnostic tasks reduce computation
  without changing the joint sampler's arithmetic, random draws or precision.
  GSL remains required for the Laplace engine.
* Earlier fitted-object schemas require refitting. Preserve archived source
  and outputs for historical replay; renaming an object cannot generate missing
  draws or change their target. See the migration and sampler guides.

# IntegMultiReg 0.2.0

## Interface and statistical meaning

* `imr()` accepts validated `imr_data` objects, platform lists and formula/data
  calls. `model_variant = "imr"` enables MRF borrowing; `"bms"` fits subgroups
  independently. Normal fitting uses the symmetric MRF factor, boundary
  Hastings correction, negative Gamma rate term and coefficient-role precision
  assignment. Historical update and precision choices are no longer public
  fitting options.
* `inclusion_probabilities()` extracts selection probabilities. `coef(fit)`
  reports that a selection fit has no stored regression coefficients; it does
  not return probabilities under a coefficient name or launch MCMC implicitly.
* `sample_regression_posterior()` explicitly performs additional conditional
  coefficient/variance sampling. Its `output_draws` and
  `min_draws_per_model_chain` distinguish returned samples from chain effort.
  Diagnostics identify the `selection_model` and actual `draws_per_model_chain`.
  `coef()` on the returned object reports coefficient posterior means.
* `selection_summary()` summarizes selection indicators and MRF interactions.
  `confint(fit)` requires an explicit parameter category. On regression-draw
  objects, `summary()` and `confint()` support coefficients, residual variance
  and retained latent responses.
* `compare_fit_summaries()` describes fit structure, response scale and selected
  feature counts. `validate_imr_object()` checks stored object structure; neither
  operation claims predictive superiority or convergence.
* Prediction lists use `subgroup:...` keys and print measured platform names.
  Regression-posterior predictions choose `quantity = "conditional_mean"` or
  `"new_observation"` to distinguish parameter uncertainty from outcome noise.
  Log-time survival point summaries remain medians because a time-scale
  posterior mean need not exist.

## Cross-validation

* `cv_method = "refit"` is the default. It recomputes preprocessing, formula
  transformations, selection MCMC and prediction within each training fold.
* `cv_method = "reweight"` reuses a full-data fit with inverse-density weights.
  Its independent `model_set` choice is `"all_draws"` or `"top_unique"`.
  `ridge`, predictive df and fold replay remain explicit numerical controls.
  The public algorithms use fixed standard AUC, concordance and MSE definitions.
* CV results have class `imr_cv`. Printing describes the actual calculation,
  state collection, model cap and scoring rule. Controls record model variant,
  update rule, effective settings, actual folds and refit seeds.
* Both algorithms support PSOCK workers. Their own seeds, partitions and output
  order remain stable across worker counts. Bounded caches reuse repeated
  states without changing their empirical multiplicities.

## Objects, documentation and native validation

* Fit schema 3 uses the named `control`, `model`, `preprocessing` and `posterior`
  sections. `upgrade_imr_object()` explicitly converts complete saved fits and
  regression-draw objects without resampling. Historical draws remain
  inspectable but require their archived source for numerical replay.
* Short verbose chains report progress correctly. Censored working times that
  exceed the historical latent-proposal bound are rejected before native
  workspace allocation. Native checks cover numerical failures, allocation
  ownership and cleanup. An unused positive-coefficient branch was removed.
* The reference includes equations and plain-text alternatives. Long web
  equations use display layout. Tutorials lead with data objects, formula
  fitting and paired CV, and include a complete nested-selection guide.
  Website checks validate links and parse TeX with KaTeX.

The function and argument renames are intentional breaking changes within the
0.2.0 development series. See `inst/MIGRATION.md`; retired names are not aliases.

# IntegMultiReg 0.1.4

* Reject platform names `id` and `subgroup` in data construction and validation
  to prevent collisions with availability metadata and misleading CV errors.

* Add a runnable paired-CV and nested formula-selection example, with saved
  subject/fold assignments and repeatability checks in CI.

* Reject overflowing integer controls, ambiguous data-frame column names, and
  corrupted saved-fit iteration counts or subgroup mappings before computation.
* Add public API execution auditing and boundary regressions; CI retains coverage
  and native sanitizer logs as downloadable artifacts.

- Replaced the post-fitting CV approximation with independent training-fold
  MCMC fits, preprocessing and model weights. The previous inverse-density
  weights followed an importance-sampling correction; their use of held-out
  responses alone does not establish an implementation error.
  Runtime now scales with folds times rounds; old prediction tables must be
  regenerated. CV preserves the caller's RNG and exposes out-of-fold records.
- Formula fits retain raw formula data so CV can rebuild transformations using
  only training subjects. Older formula fits require refitting for CV.
- Corrected survival concordance for tied predictions and incomparable pairs;
  undefined fold metrics return NA.
- Binary point prediction now averages model-specific probit probabilities.
- Corrected theta interval labels when a platform appears in four or more groups.
- Added deterministic leakage, probability-averaging and pair-order regressions,
  with concordance comparisons against survival when available.

# IntegMultiReg 0.1.3

* Explicitly initialize conditional outcome pointers for strict compiler diagnostics.

- Corrected survival preprocessing to log positive event/censoring times once,
  matching the original supplementary C code. Survival fits must be rerun;
  `survival_scale = "identity"` explicitly reproduces the historical raw-time
  implementation. The CV partitioning and selection sampler order are unchanged.
- Added optional `sample_regression_posterior()` with subgroup coefficient intervals and
  posterior predictive intervals. All active coefficients, including intercepts
  and clinical effects, use the original pMOM prior. Clinical effects remain
  always included; automatic clinical variable selection is not implemented.
- Conditional sampling reaugments binary and censored responses and reports
  classical split R-hat. Model averaging retains the fitted selection sampler's
  Laplace approximation; it is not an exact model-weight calculation.
- Survival posterior prediction uses time-scale medians as point summaries
  because inverse-gamma variance mixtures need not have finite time-scale means.
- Added independent analytic/integration checks and survival-scale regressions.

# IntegMultiReg 0.1.2

## Native-code reliability

- Fixed premature unprotection of R return objects in the fitting and
  prediction entry points.
- Fixed an allocation leak in binary cross-validation prediction.
- Initialized cross-validation work arrays and moved input-sized buffers off
  the C stack.
- Added dedicated valgrind CI and macOS ARM vignette checks under ASAN/UBSAN.

## User interface

- Added the validated `imr_data` class for training and prediction inputs.
- Added an `imr()` formula/data method while preserving the original interface.
- Formula predictions reuse training transformations, factor levels and
  contrasts; subject IDs are excluded from dot-formula expansion.
- Unsupported no-intercept and offset formulas now fail explicitly instead
  of silently fitting a different model.
- Added fitted-object validation, fit comparison, posterior uncertainty
  summaries, credible intervals and parameter-specific trace plots.
- `predict.imr()` now returns values at full numeric precision.

## Reproducibility

- Expanded the standalone replication script to execute every R command shown
  in the manuscript and to fail when full-run reference results drift.
