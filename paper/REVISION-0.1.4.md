# JSS revision and correctness fixes — 0.1.4

The editor's comments were read from `comments.rtf`. This revision fixes the
five issues identified in the source review, while retaining the data class,
formula interface and uncertainty facilities already implemented in 0.1.3.
The package fixes and CI updates are published on the isolated GitHub branch
`revision/jss-0.1.4-api-audit`. The main branch is unchanged; no editor message
or submission has been sent.

## Correctness changes

1. `cv_imr()` now makes an independent MCMC fit per training fold. It recomputes
   standardization, formula transformations, selection models and posterior
   model weights using training subjects only. The former native CV path was
   removed: changing its weight sign alone would leave response leakage and
   full-data selection leakage. Prior settings and MCMC lengths are retained;
   groups admitted by the original size threshold remain in the training folds.
   Tuning these settings still requires another validation layer.
2. The corrected CV path predicts using training posterior weights. It no
   longer rewards large held-out errors or uses held-out outcomes as weights.
3. Concordance counts each comparable pair once, awards half credit to tied
   predictions and excludes incomparable pairs. No comparable pairs returns NA.
4. Binary model averaging now calculates `sum(w * pnorm(eta))`, instead of
   `pnorm(sum(w * eta))`. Returned values retain full numeric precision.
5. Theta posterior summaries follow the lower-triangle ordering used by C,
   fixing subgroup-pair labels for platforms present in four or more groups.

CV preserves the caller's R RNG state, provides reproducible folds and sampler
seeds, and retains subject-level out-of-fold records in the `predictions`
attribute. Single-class AUC and undefined concordance remain NA rather than
being silently discarded. Old fits without raw inputs, and old formula fits
without raw formula data, must be refitted before CV.

## Editor-comment mapping

| Editor concern | Current implementation and verification |
| --- | --- |
| Dedicated validated data type | `imr_data`, validation, print/summary, availability extraction, and prediction input; existing data-interface tests pass. |
| Inspect, validate and compare fitted models | `validate_imr`, `compare_imr`, theta/selection traces, posterior summaries and intervals; existing diagnostics and posterior tests pass, and theta labels are corrected. |
| Covariate selection | Formula covariates are explicit forced adjustment terms. A paired formula comparison and nested CV selection example now supplement descriptive fit comparisons, without a new clinical inclusion prior. |
| Executable manuscript commands and printed output | The current replication script includes all nine manuscript CodeInput blocks, including both print calls, fit summaries, predictions and KIRC summaries. |
| Prediction-table discrepancy | The script requires 0.1.4 and uses a separate versioned reference. New CV results must replace the old post-fitting table; the frozen 0.1.3 archive and checksums are preserved. |
| Broader software overview | A sourced table covers MOFA2/mofapy2, mvlearn, MATLAB/R SNF and Julia CCA, distinguishing their modelling targets and input requirements. |
| Avoid rounding predictions | Full precision is preserved; probability averaging is additionally corrected. |
| Formula/data interface | Existing formula support is retained. Raw formula data is now stored so CV can learn transformations using training subjects only. |

## Verification

- Current source compiled successfully on macOS ARM, R 4.4.2, GSL 2.8.
- `R CMD check` completed with Status OK: no errors, warnings or notes.
- All 467 test expectations passed, with no skips or warnings.
- New regressions change held-out responses and features, verify that training
  fits and held-out predictions remain independent of held-out responses,
  compare concordance with `survival::concordance`, check training-only
  polynomial bases, and test binary probability averaging and theta labels.
- Package examples, vignette rebuilding and the PDF reference manual passed.
- All 14 exported functions and 21 registered S3 methods were executed under
  instrumentation; local line coverage is 89.42% (R 91.28%, C 87.63%).
- Additional boundary tests cover integer overflow, ambiguous column names,
  identifier collisions, corrupted saved fits and one retained draw/zero burn-in.
- All five selected GitHub Actions workflows (ten jobs) passed. Linux release/
  devel, macOS and Windows package checks report Status OK. Linux/macOS ARM
  sanitizers passed; R-devel Valgrind passed 467 expectations with zero errors,
  zero definitely/indirectly/possibly lost bytes and no suppressions.
- See TEST-AUDIT-0.1.4.md and the CI-RESULTS.md output for the detailed coverage
  matrix, current GitHub run identifiers and limitations.

Logs are stored in `output/revision-0.1.4/`. The package source tarball is
`IntegMultiReg_0.1.4.tar.gz`. Two independent full default replications completed successfully. All six CSV
tables matched byte-for-byte, and the saved unrounded prediction metrics were
identical. The standalone quick workflow also passed outside the project. After the API
boundary audit, another full default run of the final source again matched all
six CSV files and unrounded metrics exactly.
The manuscript table and its interpretation have been synchronized with these
results. See `output/revision-0.1.4/comparison.json` and the replication logs.
The obsolete full-scale KIRC table has now been removed from the manuscript and
current replication workflow and preserved in a separate historical audit.

## Deliverables

- `IntegMultiReg_0.1.4.tar.gz`: installable source package.
- `output/revision-0.1.4/replication/`: standalone scripts, source, references,
  provenance and validation scope, with SHA256 checksums.
- `output/pdf/IntegMultiReg.pdf`: revised manuscript.
- `output/pdf/response-to-editor.pdf`: revised point-by-point response.

Cross-validation now costs approximately one full fit per fold and round.
The validation intentionally does not reuse full-data selection or latent draws.

## Editor-completion addendum

The three further requests are completed in the current manuscript and reply.
See `../output/editor-completion-0.1.4/COMPLETION.md` for the current source SHA,
five new successful CI jobs, fourteen-table replication evidence and submission
bundle. The earlier six-table runs above refer to the pre-addendum workflow;
the computational-scale table intentionally changed with the new manuscript scope.
