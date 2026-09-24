# Anvil implementation, 2026-09-23

This is the active development workspace restored from the independently
verified upload. The original workspace is unchanged. The archive SHA256 is
`c092eeeca4203525bfc4767be714aa1b1e28ea1d7ff5f6921173f437a5a8cbc6`;
all 3381 payload files matched before development began.

## Implemented; validation pending

- Theta proposals that round to zero are rejected; representable Gamma densities
  retain historical arithmetic, with log-density fallback for underflow/overflow.
  NaN acceptance ratios cannot accept a proposal. IMR interaction states are
  checked for positivity and symmetry.
- Refit CV generates the original seed plan even for supplied folds, and saves
  the rounds-by-fold fitting seeds in `control$refit_seeds`.
- Fit-time Laplace calls, iteration limits, nonfinite likelihoods and failed
  factorizations are counted by subgroup and initial/selection/latent stage.
- Canonical research runners contain task selectors and explicit method/CV
  argument overrides. Strict paper auditing retains failed-fit checkpoints and
  independently checks a dimensionally singular training design.
- `hpc/` contains the active Anvil entrypoints. The transferred MSI scripts are
  retained under `provenance/msi-before-anvil/`, not maintained as a second runner.
- The task manifest has 679 rows: 10 strict audits, 660 simulation tasks,
  8 diagnostic chains and one original-scale Table 1 rerun. Twenty scientific
  tasks form the pilot subset without changing full-study repetition counts.

## Submitted jobs

Current job IDs are also saved in `runtime/current-jobs.tsv`.

- `20884522`: locked Linux dependency build, shared, 16 GiB / 9 allocated CPUs,
  two-hour limit. Returned from debug to shared after shared began scheduling.
- `20887913`: candidate-b installed suite, exact defaults, sampler/original-C
  references, checkpoint tests, real task splitting and scheduler negative tests.
  Depends on dependency installation; 16 GiB / four hours.
- `20884942`: candidate-b ASAN+UBSAN suite and allocation-failure checks;
  depends on dependency installation; 16 GiB / four hours.
- `20887909`: level-2 R-devel r90579 / Valgrind 3.27.1 runtime build and 45 locked
  validation dependencies; 16 GiB / eight hours; fresh runtime directory `native-runtime-2`.
- `20887914`: strict examples/tests/vignettes, graphics and negative controls,
  instrumented parent and PSOCK workers; depends on all validation prerequisites;
  16 GiB / twelve hours.
- `20887915`: study bootstrap, depends on strict-check success. Freezes the
  tested research sources and runtime, then schedules audit → pilot → full →
  collection through evidence-checking controllers, capped at 16 tasks.

`20884941` failed before compilation because the official R tarball did not
match the historical `ab0444...` pin. The downloaded archive is SHA256
`080d29d0791b77df9a1e856fff16160e48ec744fa931baf84afde40fe27b154c`.
Its VERSION, memory.c, Defn.h and RNG.c independently match official SVN r90579.
The audit is in `runtime/native-runtime-1/source-identity-audit.json`.
This establishes a new runtime pin, not equivalence to the former archive;
all instrumented validation must run again. The failed build remains intact.
Candidate-b package R/C sources are unchanged; only the active CI runtime
builder gained the newly verified pin and a local-archive input.

Obsolete pending jobs `20884678`, `20884868`, `20885029`, `20886175` were cancelled
and replaced. Slurm copies the main batch script at submission, so the replacement
regression job includes the final source-hash guard and scheduler tests.
Unrelated user jobs are untouched.

No scientific tasks have been submitted. Do not create `UNIT-VALIDATED`,
`NATIVE-VALIDATED`, `AUDITS-ACCEPTED` or `PILOT-ACCEPTED` manually; each requires
the corresponding successful checks. Existing historical CI does not validate
the repaired native sources.

Offline checks passed: syntax for all C sources and R orchestration files;
complete task coverage and CLI overrides; phase gates and manifest mutation;
mock Slurm submission/resource bounds; once-only OOM/timeout retry selection;
job record timing, warnings, RNG preservation and failure history. Mock scheduler
checks do not submit jobs and do not establish that the numerical suite passed.
The full installed and instrumented suites have not passed yet; see the latest scheduling update below.

Scratch quota check: 100 TB allocation, 63.9 GB used at the recorded check;
CPU allocation balance: approximately 741,927 SU. These are point-in-time
observations, not reserved resources or completion-time guarantees.

## Remaining before launch

Run and inspect candidate-b regression, real task-splitting, package check and
instrumented checks. The queued bootstrap opens scientific phases only after
those checks pass; failures stop the dependent chain. Scientific convergence,
full-study validation, final manuscript/response refresh and release packaging
remain outstanding. Native mid-chain continuation is not implemented.

## Multi-partition routing update

At the user's request, pending jobs were given multiple eligible partitions
using Slurm's native partition list. Job IDs and dependencies were preserved;
no duplicate build or study jobs were created. Details are in
`runtime/partition-routing-20260923.json`.

Dependency installation 20884522 and native runtime build 20887909 started on
`highmem`, node b000. The new R and Valgrind archive checks passed. Regression
and sanitizer jobs can use shared/highmem; short bootstrap can use shared/debug.
Strict validation remains on shared while the four-job highmem submission
limit is occupied. The highmem CPU billing multiplier is 4; whole-node and GPU
partitions were not added. Production study parameters and the 16-task cap
remain unchanged. Check live scheduler state before any subsequent rerouting.
