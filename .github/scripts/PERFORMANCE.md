# Performance development evidence

Baseline source: `5bd24487ae94d151f3f7faf715115123b3cd3792`.
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

## Current status

Only measurement infrastructure is accepted so far. No runtime improvement,
parallel CV interface or complete performance acceptance is claimed yet.
