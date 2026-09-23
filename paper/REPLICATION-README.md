# IntegMultiReg 0.2.0 replication

The maintained scripts use the current API. Previous scripts and reference
values remain under `historical-audit/` and `expected-results-0.1.4/`.
`make_analysis.R` forwards its arguments to `IntegMultiReg-replication.R`;
its former implementation is archived under `historical-audit/pre-flexible-0.2.0/`.
Do not compare results from different sampler/CV conventions as if they were
repeat runs of one method.

## Install and reproduce the manuscript

```sh
R CMD INSTALL IntegMultiReg_0.2.0.tar.gz
Rscript paper/IntegMultiReg-replication.R --quick
Rscript paper/IntegMultiReg-replication.R
Rscript paper/execute-manuscript.R --out-dir command-audit
```

The manuscript main fit and predictive comparison explicitly use
`sampler_method="paper"`. Predictive comparison requests `cv_method="refit"`.
The covariate-formula example retains the documented legacy sampler in both
candidates, with paired/nested refit CV, for continuity with its separate
reference. Its purpose is to demonstrate predictive formula selection.

Normal replication retains 8,000 draws after 2,000 burn-in for figures and
4,000 after 1,000 for the six outcome/model comparisons (10 rounds of 5-fold
CV). Raw precision predictions and CV objects are saved alongside displayed
tables. `execute-manuscript.R` executes every CodeInput block and produces
`code.html`, individual code/output files, and a machine-readable command map.
`--refresh-output` updates adjacent CodeOutput blocks after successful execution.

The current `expected-results-0.2.0-paper/` reference was accepted after two
complete runs of the same candidate, independent score/fold audits, agreement
of 23 CSV/RDS artifacts at tolerance 1e-12, and pixel-identical figures. Its
provenance records the package and evidence checksums. This is a manuscript-scale
reference, not a reconstruction of the original study or historical mismatch.

To establish a new numerical reference, write to a new reference directory
using `--update-reference --reference-dir DIRECTORY`, then rerun without
`--update-reference` into a second clean output directory. Inspect the causes
of changes before promoting a reference. `--quick` cannot establish a reference.
The 2017 paper/code differences and the new sampler corrections are described
in the package's `METHOD-COVERAGE.md`; replacing an old CSV does not resolve a
historical discrepancy by itself.

## Original study experiments

```sh
Rscript paper/prepare_biom12587_table1_data.R
Rscript paper/original-experiments.R --experiment table1 --reference code2017 --quick
Rscript paper/original-experiments.R --experiment simulation --reference paper --quick
Rscript paper/original-experiments.R --experiment correlated --reference code2017 --quick
```

Omit `--quick` for original-scale settings: 350,000 retained draws after 50,000
burn-in; Table 1 uses ten rounds of ten-fold CV; main simulations
use 50 replicates per configuration (main Figure 3), and correlated simulations
use 30 (Web Appendix E, Figure 5). The default simulation prediction is one
10-fold round; the article does not specify every historical seed and folding
choice. Set `--rounds` explicitly for additional rounds. Other execution controls
are `--draws`, `--burnin`, `--k`, `--replicates`, `--workers`, `--seed`, `--data`
and `--out-dir`. Use one parallel layer and account for memory per worker.

`reference` may be `paper`, `code2017`, or `package`. These select named
argument lists passed to the same public functions. Completed model jobs are
checkpointed; a directory can be resumed only with identical settings, data
and source hashes. Store fits, folds/order, predictions, per-validation scores,
selection AUC, means, SD, SE, effective counts, timing and session information.
Table 1 uses **fold-level** averages, not pooled scores or the SD of round means.

The independent original generator is compiled from the supplied supplement.
Only the unavailable malloc header and an auxiliary RNG leak are patched;
source hashes and the changes are recorded. For simulation and correlated
experiments, `--marker-design random|fixed` selects signal identities separately
from fitting conventions. The paper reference defaults to random markers;
code2017 and package default to the released fixed indices. The random design
applies one seeded column permutation per platform before the archived generator,
then restores feature order and truth labels. This preserves signal counts and
100%/50%/0% overlap without a second response-generating algorithm. Actual
permutations, the R RNG convention and seed are saved with the generated data;
caller RNG state is preserved. The original auxiliary GSL seed is retained.
Exact historical replicate seeds are unavailable. These
facts limit claims about reproducing published digits. The 778/776 gene-count
inconsistency is also retained in provenance. Full CPH and L1-CPH benchmarks
are separate from IMR/BMS; Uni-CPH uses Hommel-adjusted per-platform tests.
Cox convergence/separation warnings remain in run logs and must be reviewed.
Each model and Cox job also writes a `.status.rds` and numbered `.attempt-N.rds`
record with start/end times, elapsed CPU/wall time, warnings and any error.
Cox benchmark results include a `diagnostics` list identifying each subgroup,
round and fold, training sample/event counts, warning messages and (for the
unpenalized model) named coefficients. Inspect these when convergence warnings
occur; finite predictions alone do not establish a stable coefficient estimate.
L1-CPH diagnostics also record the full-data selection stage, actual inner-CV
ID/fold assignments, selected lambda, CV path and full-fit error code. Warnings
from a truncated glmnet path remain visible and require checking whether the
chosen lambda is supported by every inner fit; a completed job alone is not
evidence that every attempted lambda converged.
Only jobs with both output and a completed status are skipped on resume;
failed or interrupted attempts are retained and rerun.
Fitting is checkpointed before CV, so a failed CV attempt can reuse its saved
fit when the manifest is unchanged. The temporary fit checkpoint is removed
only after the complete model result and completed status have been saved.
CLI flags must be known, unique and supplied with a value where required. R heap measurements exclude
native allocations and worker memory; they are not total process peak RSS.

Preparation saves **raw positive months** for public fitting with
`survival_scale="log"`; an audit-only `outcome.log` stores the transformed
response. The runner recognizes older prepared objects by using `outcome.raw`,
avoiding a second logarithm. Correlated simulations preserve the transformed
original generator output with `standardize=FALSE`.

Full-scale MCMC retains a large joint-state history. The earlier unshared integer
representation required approximately 6.4 GiB for 350,000 retained states before
native copies and workspaces. The current candidate shares repeated matrices
without thinning states, and the research writer packs them losslessly for disk.
The 16 GiB machine has completed the code2017 Table 1 run and one original-scale
replicate in each main simulation scenario; these results do not establish the
resource requirements of every reference mode or correlated configuration.
Pilot and manuscript demonstration runs remain distinct from complete simulation
studies. Appendix F chain interpretation and Appendix G coefficient/network
comparisons remain separate work; proprietary IPA networks cannot be recreated
from the bundled public inputs alone.

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

Research fit checkpoints now use a versioned, lossless matrix-pool container;
`digest` is required by the runner to build short indices without deparsing
large matrices. Every mapped matrix is checked with `identical()` before saving.
All draws and their order are retained, and loading restores ordinary fitted
objects with shared matrices. Older plain-fit checkpoints remain readable but
may require the original unshared memory. The container is a runner persistence
format, not a public fit API or a thinning scheme. Final result RDS files remain
ordinary objects and may expand in memory when reloaded.

## Appendix F diagnostic export

Source `paper/appendix-chain-diagnostics.R`, then call
`appendix_chain_diagnostics(fits, out_dir)` with a list of matched fitted objects.
It exports platform/subgroup feature-mPIP correlations, classical coda
Gelman–Rubin PSRF with `autoburnin=FALSE`, and per-chain theta autocorrelations.
The matrices already contain retained draws. Undefined constant-vector
statistics remain NA. This exporter does not establish convergence: independent
long chains, trace plots and scientific review remain required. The historical
eight initial states are unknown; the short eight-start fixture checks the
interface only.

Run the eight-chain study sequentially with:

```sh
Rscript paper/run-appendix-chains.R --quick --out-dir output/appendix-smoke
Rscript paper/run-appendix-chains.R --out-dir output/appendix-full
```

The full default is 350,000 retained draws after 50,000 burn-in, on the original
KIRC matrices with clinical covariates and the paper sampler. Seeds are 100–107
by default. The eight explicit starts vary selection and positive interaction
values; they are documented new choices, not recovered historical starts.
`--data`, `--draws`, `--burnin` and `--seed` override those settings. Runs retain
lossless full-fit checkpoints, compact diagnostic projections and per-chain
initial states/status records. Resume requires unchanged settings, data and
source hashes. Keep the installed package library fixed while a run is active.

After all eight chains and diagnostic tables are complete, source
`paper/plot-appendix-chains.R` and call `plot_appendix_chains(run_dir)`.
The four PNG figures show full and retained log-posterior traces, theta ACF,
and mPIP correlations. Trace panels share axis limits within each figure;
correlation and ACF panels use a common [-1,1] scale. All log-posterior draws
are plotted. Short-run figures are explicitly labelled as workflow checks.
The plotter runs separately so its layout can be refined without changing the
frozen sampler run or its resume hashes. Reinspect figures at full-run scales
before using them in a manuscript.

For Appendix G selection rankings, source `paper/appendix-marker-rankings.R`
and call `appendix_marker_rankings(fit, kirc_full, out_dir)`. It ranks each
platform/subgroup separately, taking 10 genes, 6 miRNAs and 10 methylation
probes, and preserves original archive feature names. Ties use original column
order and are flagged at the cutoff. Counts across top lists are descriptive.
This exports the selection part only: the published coefficients are posterior
modes, and neither `coef(imr)` nor posterior coefficient means may substitute
for them. The exact historical table-specific conditioning remains under audit.

New original-experiment model results use a versioned lossless container to avoid
serializing repeated selection matrices separately. Source
`paper/original-experiment-helpers.R` and use `read_experiment_result(path)` to
restore the ordinary result list and full fit; `include_fit=FALSE` reads only
the other result fields for aggregation. The reader also accepts older ordinary
result RDS files. Do not access a new container with `readRDS(path)$fit`.
The Table 1 validator handles both formats. This changes storage only: all
retained states, their order, folds, predictions and scores are preserved.
New result containers use xz compression. Allow additional serialization time;
compression and exact restoration do not constitute scientific result acceptance.

Audit a finished simulation manifest with
`Rscript paper/validate-original-simulation.R RUN_DIRECTORY`.
The audit verifies recorded source/data hashes, fit settings, saved folds,
subject coverage, predictive scores, and independently computes Bayesian
selection AUC by positive/negative feature comparisons. It accepts both result
formats. Add `--require-complete-study` to require the original iteration budget
and all 50 main or 30 correlated replicates. Without that flag, passing validates
only the executed manifest; it does not promote a one-replicate scale check to
a completed study. Benchmark numerical warnings still require scientific review.

Also run `Rscript paper/validate-simulation-baselines.R RUN_DIRECTORY` to
reconstruct the L1-CPH selection fits using their saved inner folds and check
lambda.min and the cross-validation paths. This separate audit refits every
univariate Cox model, applies the recorded Hommel adjustment, and independently
recomputes both benchmark selection AUCs by pairwise feature comparisons.
It saves acceptance records and feature-specific warnings without altering the
original results. Matching rerun values do not establish the statistical
validity of a warned fit; review the warnings before interpreting comparisons.

After both audits pass, run
`Rscript paper/summarize-original-simulation.R RUN_DIRECTORY`. It cross-checks
exported model CSVs against their saved results and writes
`replicate-metrics.csv` plus `across-replicate-summary.csv`. Prediction folds are
averaged within each replicate first; the reported SD and SE then describe
variation across independent simulation replicates. A single-replicate scale
run therefore has `NA` for both SD and SE. Requested and original-design repeat
counts remain separate. The saved summary provenance also records score-definition
differences; in particular, the L1 overall score is a weighted subgroup score.

For an additional rank-based interaction-parameter review, run
`Rscript paper/rank-chain-diagnostics.R RUN_DIRECTORY OUTPUT_DIRECTORY` after
all eight chain statuses are complete. It uses the installed `posterior` package
for rank-normalized split/folded Rhat and bulk/tail ESS, retaining every saved
post-burn-in draw. This complements, rather than renames, the original classical
coda PSRF tables. A final integer argument explicitly requests an interim review
of the first 2–8 completed chains; fewer than eight is labelled incomplete.
The output records input/script hashes, package version and method documentation.
Undefined statistics remain undefined and are never treated as passing.
