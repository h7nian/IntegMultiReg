# Performance development evidence

Baseline source: `5bd24487ae94d151f3f7faf715115123b3cd3792`.
After the maintainer-approved attribution-only history rewrite, this baseline
is `caf9f622a65e29a71a9502994635eb7554567700`. Their source trees are identical.
Existing evidence retains its original SHA; it is not relabelled as a new run.
Keep the baseline and candidate installed in separate libraries built with the
same compiler, flags and numerical dependencies. Do not load both into one R
process. Scripts here are development tools, not published inference results.

`performance-run.R` runs one isolated scenario and saves stage timings and
numerical outputs. `performance-repeat.R` freezes a copy of that worker and
runs one warm-up followed by five serial measurements. On macOS its individual
logs include `/usr/bin/time -l` peak RSS. Other platforms currently save timings
only; their native memory/profiling collectors remain to be implemented.

From the package root:

```sh
IMR_SOURCE_SHA=5bd24487ae94d151f3f7faf715115123b3cd3792 \
Rscript .github/scripts/performance-repeat.R /absolute/new-output \
  /absolute/baseline-library sim continuous imr 2000 1000
```

Replace `sim` with the local converted supplementary `.rda` path and use
`right.censored` for KIRC. Data never need to be uploaded to CI. Set
`IMR_PERFORMANCE_STAGES=fit` to measure fitting alone; otherwise the worker also
runs repeated prediction, all three CV modes and short conditional chains.
Short conditional-chain diagnostic warnings are deliberately retained in logs.

`performance-compare.R baseline-output candidate-output` enforces exact saved
fit posterior/model/preprocessing and conditional draws, with the agreed
prediction and metric tolerances. This is an initial differential gate, not
the complete acceptance suite: RNG snapshots, native model-ranking/count
evidence, formula environment comparison, controlled cache scenarios and
worker fault tests still need dedicated coverage.

## Rejected candidate

Splitting the coordinate-ascent dot product around its diagonal preserved the
apparent source-level addition order, and all 719 existing assertions passed.
However, a 100-draw, 50-burn-in continuous IMR differential run detected changes
in `log_posterior` (reported mean relative difference about `4.1e-16`). Rebuilding
the unmodified baseline reproduced its results exactly. The candidate was
removed: exact sampler outputs are required, not merely unchanged selection
indicators. No tolerance was relaxed.

## First optimization: importance state reuse

The implementation indexes immutable, full selection states once per native
call. Each hash hit is verified against every indicator. Only state identities
are shared across folds: coefficients, predictions and weights are fold-local.
Every draw retains its original position in the numerical accumulation.
The auxiliary index/table payload is bounded at 128 MiB; a half-full table
retains its first entries and uncached states use the ordinary calculation.
Insufficient index budget disables caching. Internal adapter arguments permit
small budgets and forced collisions for tests, without a public tuning option.

On this macOS ARM machine, the same frozen KIRC fit (2,000 draws, 1,000 burn-in;
847 unique full states) gave the following 5-fold, 2-round importance timings.
Each value is the median of five repeats following a separate warm-up:

| Selection-state workload | Baseline seconds | Indexed seconds | Change |
|---|---:|---:|---:|
| Fitted KIRC states | 5.476 | 2.431 | 55.6% less time |
| Artificial all-unique states | 5.141 | 5.167 | 0.5% more time |
| Artificial all-repeated states | 5.021 | 0.176 | 96.5% less time |

These are local workload measurements, not universal speed guarantees or
published predictive results. The artificial fixtures are not posterior draws
for inference. An earlier per-fold hash implementation was rejected because it
increased the all-unique workload by about 11%; the retained implementation
hashes once per call. The final local test suite passes 737 assertions with zero
failures/warnings/skips, including cache-off, bounded-capacity and collision tests.

## Remaining acceptance and implementation

The refit path now uses a deterministic task dispatcher. Its internal PSOCK
route is tested at 2 and 3 workers against the serial result for all six
outcome/model combinations, plus formula transformations, serialization,
RNG preservation, worker errors and socket/environment cleanup. Post-fit
partition/prediction/scoring stages now preserve the full native output under
nonconsecutive fold batches, including survival's carried-over shuffle state.
Both post-fit modes pass exact 2/3-worker comparisons on all six simulated
outcome/model combinations and the frozen 448-subject KIRC fit (5 folds,
2 rounds). The public `workers` argument is now connected to all three modes.
`performance-refit-workers.R` measures the refit
dispatcher (one warm-up and five repeats at each of 1, 2 and 3 workers); use the same
frozen fit and an otherwise idle machine. Short tasks may not amortize process
startup, and each child holds its own fit and training-fold allocations.
The public worker implementation and bootstrap/cleanup fixes passed 1,049
installed assertions with zero failures, warnings or skips, including an
installation that retains R source references. Eighteen saved-baseline comparisons (all CV modes
for all outcome/model combinations) return exactly identical complete results.
The public API also validates worker counts and relays worker warnings in task
order. Parallel formula refits require serializable custom-function bindings
or package-qualified functions; the caller's global workspace is not exported.
These are correctness gates, not acceptance of the complete performance plan.

### Public worker measurements and validation

The following post-fit measurements use commit `eb9b213`, one warm-up followed
by five repeats, 5 folds and 2 rounds, and the same frozen 448-subject KIRC fit.
All saved complete CV results are exactly identical to the isolated baseline.
Times are medians in seconds; raw timings retain the individual ranges.

| Workload | Baseline serial | Optimized serial | 2 workers | 3 workers |
|---|---:|---:|---:|---:|
| KIRC legacy | 0.918 | 0.897 | 1.289 | 1.411 |
| KIRC importance | 5.366 | 2.418 | 1.866 | 1.668 |
| Artificial unique importance | 5.095 | 5.103 | 3.286 | 3.006 |
| Artificial repeated importance | 4.860 | 0.178 | 0.693 | 0.866 |

Process startup makes short workloads slower; the default remains one worker.
The artificial states are not inferential posterior samples. These timings
predate the bootstrap and error-cleanup fixes; they must not be relabelled as
measurements of the later source. Full-size refit timing remains separate.

CI exposed two worker startup issues: retained source metadata could load a
different installation before worker library initialization, and macOS PSOCK
children needed sanitizer preloading before R startup. Commit `227ada1` strips
bootstrap source references and closes each worker independently after failure.
Commit `942762d` adds a CI-only macOS launcher; production execution does not
have a sanitizer-specific path. At `942762d`, Linux/macOS/Windows checks,
GCC/Clang ASAN/UBSAN and macOS ARM ASAN/UBSAN pass. The macOS sanitizer log executes all 1,049 assertions
with zero failures, warnings or skips. No failed test was disabled.

### Completed worker-phase acceptance

Run 34670333384 at `28e72b8` passes the complete Valgrind suite: 1,049 assertions,
zero failures/warnings/skips, zero parent Memcheck errors, zero definite leaks
and zero suppressions. All 136 instrumented suite workers also pass their
individual memory reports. A separate focused job covers 36 workers across
all outcome/model/CV combinations. These commits after `942762d` change CI
tooling only, not package execution.

The initial parent-only Valgrind run had small exact-refit comparison failures
when its children ran outside Valgrind. Instrumenting both parent and children
resolved those failures without changing production code or numerical tolerances.
Startup helper forks also produce log files; launcher/process PID matching
identifies the actual worker reports, all of which must finish and pass.

Final-runtime KIRC measurements at `942762d`, again with one warm-up plus five
repeats, give the following median seconds:

| CV | Baseline serial | Candidate serial | 2 workers | 3 workers |
|---|---:|---:|---:|---:|
| Legacy | 0.874 | 0.856 | 1.273 | 1.480 |
| Importance | 5.256 | 2.279 | 1.719 | 1.482 |
| Refit | not timed in this run | 192.177 | 101.170 | 76.328 |

The full refit return and caller RNG also match an isolated original-baseline
run exactly. All 16 final post-fit benchmark cases match baseline exactly;
all-unique serial importance is 5.018 versus 5.041 seconds. Short jobs remain
slower with workers, so the default is unchanged. The complete covariate example
was rebuilt twice: eight byte-identical CSVs, preserved RNG, unique held-out
predictions and no outer-test IDs in inner folds. The manual compiles, but its
local TeX overfull-box/font warnings still require layout review.

A separate 24-process fit-only study covers 2,000 and 10,000 draws with 1,000
burn-in. Posterior/model/preprocessing/control match exactly in every run.
Median peak RSS (whole isolated process, not aggregate CV-worker memory) is
275,972,096 versus 289,193,984 bytes at 2,000 draws and 646,905,856 versus
691,830,784 bytes at 10,000 draws. Ranges overlap substantially; this is not
evidence of memory savings. Investigate allocation behavior before optimizing.

The draw-by-test prediction buffer is still present. Compression, broader
workspace reuse, prediction/conditional-posterior profiling and appropriate
matrix-kernel experiments remain. This is not the complete performance plan;
later runtime changes need fresh validation. No sampling-kernel optimization
has been retained.

## Next optimization: numeric-column validation allocations

Rprofmem attribution on a frozen KIRC fit found repeated full data-frame to
matrix copies in `.imr_check_numeric_columns`. Ordinary columns now undergo
the same numeric and finite checks without combining them into a matrix;
classed columns retain the historical matrix-coercion path. Error priority,
invisible return value and validation coverage are unchanged.

Isolated installed-engine measurements (1,000 calls per repetition; one warm-up
and five repeats) give median seconds of 3.274 to 1.317, 0.310 to 0.130 and
2.469 to 1.005 on the three KIRC platforms. The small covariate and one-feature
fixtures did not regress in these measurements. These are validation-hotspot
times, not an end-to-end fit speed claim or a demonstrated RSS reduction.

The local installed suite passes 1,077 assertions with zero failures, warnings
or skips. All six outcome/model differential cases preserve fit posterior,
model/preprocessing, predictions, all CV modes and conditional draws exactly.
Short conditional chains retain their existing R-hat warnings and are not
convergence evidence. Fresh package checks and CI for this runtime change
remain required; the worker-phase CI results do not substitute for them.

At `117ab8f`, those fresh checks now pass, including Valgrind run 34674587041:
1,077 assertions, zero failures/warnings/skips, zero root Memcheck errors and
definite leaks, 136 verified suite workers and a separate 36-worker job.
Whole-fit KIRC 2000/1000 measurements give 23.881 seconds before and 23.012 after
(medians of five repeats following warm-up); all 12 saved fits match exactly.

## Fused training-matrix preparation

Training preparation previously computed means and SDs for storage, then
recomputed them for normalization. `.imr_prepare_matrix` now computes each
column's moments once and returns the same means, safe SDs and normalized
matrix. The original mean/SD arithmetic, degenerate-scale rule, dimnames and
drop behavior are preserved. Six obsolete private helper functions were
removed after reference scanning; prediction's known-moment helpers remain.

Installed-engine measurements (100 preparations per repetition, one warm-up
and five repeats) give median seconds of 3.359 to 1.674, 0.277 to 0.127 and
2.466 to 1.142 on KIRC's platforms. Covariate/small fixtures did not regress.
A prototype allocation profile records 58,037,632 to 28,285,744 total R-heap
allocated bytes for one platform; this is allocation volume, not peak memory.
Neither result should be reported as a whole-fit speedup.

The local installed suite passes 1,979 assertions with zero failures, warnings
or skips, including 300 randomized edge cases with exact output/input/RNG
checks. All six outcome/model full-pipeline differential cases still match
the original baseline exactly. Fresh package/CI and whole-fit measurements
for this new runtime change remain required.
