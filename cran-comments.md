## Development version 0.3.0

This revision adds joint and partially marginalized pMOM/MRF samplers while
retaining the previously corrected Laplace selection algorithm and its optional
conditional regression sampling. `marginalize` selects the regression parameters
integrated out of the transition target. Priors, MCMC settings and numerical
controls use separate constructors. The fitted-object schema has changed.

Cross-validation defaults to training-fold refits with the chosen sampler.
Optional PSIS reweighting requires regression posterior draws and reports fold
weight diagnostics. GSL remains a dependency of the retained Laplace engine.

## Validation status

Current validation is in progress on this source revision. Previous 0.2.0
source checks, manuscript replay and historical artifact audits apply only to
their recorded sources. They are not evidence for the 0.3.0 engine.

The maintained checks cover installed public APIs, independent integrated
small-model posterior targets, native memory instrumentation, worker processes,
source-package checks and rendered website mathematics. The release review
must record the actual candidate commit and outputs after those checks finish.

## Scientific scope

The joint sampler changes the computational procedure. Historical numerical
outputs remain source-qualified; replacement reference files cannot establish
why an older discrepancy occurred. Small-model correctness and repeatability
do not establish adequate real-data mixing or reproduce all original scientific
conclusions. Real-data diagnostics are reported separately.
