## Development version 0.3.0

This development revision replaces Laplace selection plus conditional
coefficient sampling with a joint pMOM/MRF sampler. Coefficients, variances,
selection indicators and interactions belong to the same retained states.
The old conditional-only API and engines are removed. Prior and MCMC options
use grouped constructors, and saved selection-only fits require refitting.

Cross-validation defaults to training-fold refits. Optional reweighting uses
PSIS on the joint observed-data likelihood and reports fold weight diagnostics.
The native engine uses R allocation and RNG APIs and no longer requires GSL.

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
