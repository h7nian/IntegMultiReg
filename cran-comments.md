## Submission

Version 0.2.0 reorganizes the public interface and adds data/formula workflows,
independent training-fold validation, and conditional regression-posterior
sampling. These are development-version breaking changes; the migration guide
lists the retired arguments and functions. No argument aliases are accepted.

## Statistical behavior

Ordinary fits use the symmetric MRF conditional, the boundary Hastings
correction, the Gamma log-density with a negative rate term, and coefficient
prior scales assigned by coefficient role. Historical unadjusted updates are
available through archived source versions, not a current public selector.
Stored historical objects can be explicitly upgraded for inspection.

Cross-validation defaults to independent training-fold refits. Post-fit
importance reweighting is selected with `cv_method = "reweight"`; the selection
state collection is specified independently with `model_set`.

`inclusion_probabilities()` returns selection probabilities. `coef()` never
labels them as regression coefficients. `sample_regression_posterior()` starts
additional conditional MCMC; its draw budget and conditional-chain effort have
separate names. The method still uses the fitted Laplace-based selection
weights. It is not the separate experimental joint posterior sampler.

## Local validation

The installed R-devel suite passes 2792 expectations with zero failures,
warnings or skips. At corresponding explicit computational settings, six fits,
eighteen CV results and six conditional regression-posterior sample sets are
bit-identical to the preceding convention-name snapshot. These comparisons do
not claim equivalence between the old defaults and the new defaults.

The source snapshot is `be38715b95568854e1a89140381ed4274f3102b2`; the built archive has
SHA-256 `d86737a6dd51ef7589af941b100969bfaf26e984becd42f1016b0fa3f89400b7`.
A complete R-devel source check passes examples (including donttest), 2764
in-check expectations, rebuilt vignettes and the PDF manual. It has zero errors
or warnings and one NOTE for unavailable local HTML Tidy/V8 tooling. CRAN
incoming feasibility is disabled for this local check. GitHub source checks
supply those tools; the website also validates equations through KaTeX.
The separate installed-suite run and the in-check run are reported with their
actual expectation counts. The final review record also binds manuscript
replay artifacts to this archive. Older CI
results are not relabelled as checks of this revision.

## Earlier computational baseline

Source `fe9ce3b7a286f2cee8f69ac0b8f25456212a79e1` passed cross-platform source
checks, GCC/Clang and macOS ARM sanitizers, and full Valgrind checks. Its
full-suite artifact contains 177 process logs and its worker artifact 36, with
zero reported errors or unfreed definite allocations. The present revision
also renames native identifiers; current native checks are recorded separately.

## Numerical reference scope

The frozen ARM Mac manuscript reference is preserved with its own source and
script checksums. Independent Linux runs agree with each other; binary AUC
can differ from the Mac reference by up to 0.002924. The archived source itself
reproduces the Linux behavior. A shared version number or a replaced reference
CSV cannot establish the cause of an older numerical discrepancy.
