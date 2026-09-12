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
dispatcher (one warm-up and five repeats at each worker count); use the same
frozen fit and an otherwise idle machine. Short tasks may not amortize process
startup, and each child holds its own fit and training-fold allocations.
The internal post-fit dispatch stage passed 999 installed assertions with zero
failures, warnings or skips. Eighteen saved-baseline comparisons (all CV modes
for all outcome/model combinations) return exactly identical complete results.
The public API also validates worker counts and relays worker warnings in task
order. Parallel formula refits require serializable custom-function bindings
or package-qualified functions; the caller's global workspace is not exported.
These are correctness gates, not final performance or cross-platform acceptance.

The existing draw-by-test prediction buffer is still present; its compression
and the full memory-growth study remain undone. This is not the complete
performance plan. Final parallel CV measurements, broader workspace reuse,
matrix-kernel experiments and final cross-platform/sanitizer/Valgrind gates
remain pending. No changes were retained in the sampling kernel.
