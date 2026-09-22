# Paper, released code, and package conventions

The package preserves its existing defaults. The optional arguments below
express separate computational decisions through the same fit and CV engines.
`cv_method` is a set of defaults; it does not override explicit arguments.

| Decision | Paper-oriented setting | Released-code setting | Package default |
|---|---|---|---|
| Selection update | `sampler_method="paper"` | `"legacy"` | `"legacy"` |
| Prior precision block boundary | `prior_indexing="standard"` | `"code2017"` | `"standard"` |
| Laplace iteration caps: initial/selection/latent/prediction | 25/40/25/40 (numerical choice) | 25/40/25/25 | 25/40/25/40 |
| CV states | `model_set="draws"` | `"ranked_unique"`, maximum 100 | depends on CV mode |
| CV coefficient penalty | `ridge=0` | `ridge=0` | 0.001 |
| Predictive df | `df_method="fractional"` | `"legacy_integer"` | depends on CV mode |
| Scoring | state the desired tie convention | `score_method="legacy"` | depends on CV mode |
| Post-fit fold random stream | no mathematical requirement | `fold_rng="continue"` | `"reset"` |
| Historical partition replay | `folds` with `row_order` | same | generated from fit seed |
| Specified chain starting points | `initial` selection/interaction matrices | historical starting values unavailable | original random selection start and fixed interactions |

The symmetric MRF has conditional log-odds `nu + 2 * sum(theta * gamma)`.
The legacy selection update uses a factor of one and omits the Hastings ratio
when a flip crosses a boundary between empty/full and interior states. Its
stationary distribution is therefore different from the stated posterior.
The paper sampler includes these corrections and a corrected Gamma rate sign
in the diagnostic log density. Legacy behavior remains explicit and testable.

`prior_indexing="code2017"` preserves an index-boundary discrepancy in the
released precision calculation; it is a historical computational option, not
an alternative coherent prior. Conditional uncertainty from `posterior_draws()`
uses the stated pMOM priors; it inherits selection weights from the fitted
sampler but does not reproduce this historical precision discrepancy.

## Examples with ordinary argument lists

```r
paper_fit_args <- list(sampler_method = "paper", prior_indexing = "standard")
code_fit_args <- list(sampler_method = "legacy", prior_indexing = "code2017",
  laplace_max_iter = c(initial = 25, selection = 40, latent = 25, prediction = 25))
paper_cv_args <- list(cv_method = "importance", ridge = 0,
  model_set = "draws", df_method = "fractional", score_method = "standard")
code_cv_args <- list(cv_method = "legacy", ridge = 0,
  model_set = "ranked_unique", max_models = 100,
  df_method = "legacy_integer", score_method = "legacy", fold_rng = "continue")
# fit <- do.call(imr, c(list(x = platforms, outcome = outcome), paper_fit_args))
# cv <- do.call(cv_imr, c(list(object = fit), paper_cv_args))
# replay <- cv_imr(fit, folds = cv$control$folds, ridge = 0,
#                  model_set = "draws", df_method = "fractional")
```

A zero penalty requires a nonsingular training design. Singular fits stop with
round/fold/subgroup context. Supplied folds match IDs, not input row positions;
`row_order` retains the original summation order. Scoring changes do not alter
predictions. Refit CV accepts scoring and supplied folds, but rejects explicit
post-fit-only options. Its preprocessing and selection are fitted on training
subjects; tuning still needs outer validation.

## What the checks establish

* Default regression: six fit configurations and eighteen CV combinations
  compared against an independently installed frozen package revision.
* Conditional CV: archived `pred_aftcv()` compiled separately; predictions agree
  to floating-point precision for fixed states, latent responses and ordered
  partitions with ridge zero and integer df. This does not test the original
  full sampler or its top-100 model ranking.
* Mathematical reference: independent linear algebra tests ridge and df choices;
  finite-state transitions verify detailed balance; long native chains are
  checked against enumerated conditional targets, and Gamma increments against
  `dgamma()` including interactions below 0.001.
* Fold replay, caches, worker execution and old-object fallbacks are tested.

These are component and algorithm checks. They do not establish reproduction
of published Table 1 or Figure 3. That additionally requires matching inputs,
preprocessing, starting states, update order, seeds, simulation replicates and
aggregation. Continuing the package's random stream does not imply that its
preceding random draws match the historical standalone program. The supplement
has 778 gene columns, whereas the article reports 776; results must state which
source was used. Full experiments and historical-digit agreement are separate
validation levels, never inferred from the presence of an option.

`initial=list(selection=..., interaction=...)` supplies named matrices by
platform. Use `coef(fit)` as the selection-layout template (subgroup rows,
feature columns) and replace probabilities with zeros/ones. Interaction matrices
use the corresponding subgroup names on both axes, symmetric positive
off-diagonals and zero diagonal; they apply only to IMR. Omitted components use
their original initialization. Specified selection avoids the original random
initialization draws, so its subsequent RNG trajectory intentionally differs.
Refit CV reuses the specified start; the values are recorded in fit/CV controls.
Eight explicitly different starts can support an Appendix F-style diagnostic,
but the exact historical starts are not supplied by the published supplement.

### Unpenalized original-data limitation found in the larger pilot

A 20,000-draw paper-sampler IMR clinical+molecular pilot on the bundled KIRC
inputs fails unpenalized all-draws CV. In subgroup 011, 22.065% of retained states
exceed the 64 training rows in five folds, and 28.060% exceed the 62 rows in the
other five folds. An independently reconstructed failing design has 69 columns,
64 rows and rank 64. These fractions are dimensional lower bounds; remaining
designs were not exhaustively rank-tested. A ridge=.001 control completes but
changes the estimator. Do not silently add a penalty, drop states or replace the
inverse when claiming the unpenalized paper reference. The separate code2017
sampler/ranked-model configuration requires its own evidence. Details and saved
artifacts: `output/memory-scale-20260921/PROGRESS.md` in the research workspace.
