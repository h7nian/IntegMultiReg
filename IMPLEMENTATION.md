# Anvil implementation, 2026-09-23

## Current checkpoint: September 24, after GitHub validation

Candidate B has passed the installed regression suite, ASAN/UBSAN and all nine
new GitHub CI jobs. Final Anvil level-2 Valgrind check `20894476` is running;
study bootstrap `20887915` depends on its successful completion. Scientific
sampling has not started and `NATIVE-VALIDATED` does not yet exist. Use Slurm
and `runtime/current-jobs.tsv` for live status. Earlier progress entries below
retain their original job IDs and observations.

GitHub branch `ci/anvil-candidate-b-20260924`, commit
`7fa16bbf5afb450324b89cc67fe3d7bcc3a535ea`, matches the frozen R/C/manual/test
source. All nine jobs succeeded:

- [Cross-platform checks](https://github.com/h7nian/IntegMultiReg/actions/runs/36010339605):
  Linux release/devel, Windows release, macOS release.
- [Native sanitizer checks](https://github.com/h7nian/IntegMultiReg/actions/runs/36010342655):
  GCC and Clang, each with UBSAN and ASAN+UBSAN.
- [macOS ARM sanitizer check](https://github.com/h7nian/IntegMultiReg/actions/runs/36010345688):
  installed tests, vignette rebuild, sampler targets and fault cleanup.

All nine suites passed 2851 assertions with zero failures/warnings/skips.
Their R builds enable one existing profmem-conditional shared-address assertion;
Anvil R 4.4.1 reports 2850. Both run the same 151 test cases. Each of the five
sanitizer jobs also passed all 72 allocation failure/recovery cases. Downloaded
logs, artifacts, commit identities and SHA256 hashes are in
`runtime/github-actions-evidence/`. The current historical-regression mapping
is `runtime/candidate-b/evidence/github-regression-audit.json`. Main is unchanged.

Strict attempt `20894290` passed nine scoped graphics controls and detected the
R-heap probe, archived 7200-byte/60-block leak and historical uninitialized read.
Its package check then reported missing `checkbashisms`. Because that warning
prevents acceptance, the attempt was stopped with its logs retained and its
temporary R test-startup instrumentation restored after termination.
The missing checker is now installed from official Debian devscripts
2.25.15+deb13u1; source/checker hashes are recorded under `runtime/tools/`.

Fresh attempt `20894476` uses `runtime/candidate-b/strict-anvil-2`, reruns all
controls and the complete package check, and verifies instrumented PSOCK
children before writing `NATIVE-VALIDATED`. Bootstrap `20887915` is explicitly
dependent on this attempt. Tested research-source hashes, the frozen package
source manifest, and all 32 dependency versions/load paths were reverified.
After validation, the bootstrap freezes `runs/study-20260923` and schedules
audit → full-parameter pilot → full study → collection with concurrency 16.

The user also reported a Claude contributor avatar on the GitHub homepage.
All seven live branches, tags (none), author/committer identities and coauthor
trailers were checked. GitHub's author filter returned zero Claude commits on
every branch; its contributor API lists only h7nian. Evidence is in
`runtime/github-contributors-audit/`. No history was rewritten. Stale GitHub
statistics are the likely explanation; disappearance of the user's UI entry
has not been confirmed. GitHub documents an approximately 24-hour refresh
period and recommends Support if the stale data persists afterward.

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

## September 24 validation status

The dependency installation completed all 32 builds, but the first checker
incorrectly compared packageVersion()'s normalized spelling (for example
0.22.6) with the literal DESCRIPTION lock (0.22-6). Exact DESCRIPTION versions
and actual loaded namespace paths now pass for all 32 packages. Finalization
job 20894124 completed successfully; the library was not reinstalled.

Regression job 20887913 completed: 2850 installed assertions, zero failures,
warnings or skips; exact default comparisons for six fits and eighteen CVs;
independent theta/Gamma targets and original-C conditional predictions; eight
serial/split chain comparisons and six serial/split simulation comparisons;
checkpoint, collection and scheduler negative controls. UNIT-VALIDATED exists.

The ASAN/UBSAN suite also passed 2850 assertions and independent sampler
references. Its final fault-injection harness assumed that the package was the
Git root; that assumption failed in this nested workspace. The harness now
copies the tested source and applies the patch in a private Git root. Separate
completion job 20894153 passed all 72 allocation failure/recovery cases against
the same candidate. SANITIZER-VALIDATED exists. Prior failure logs remain intact.

Strict job 20887914 built the package and vignette, then stopped at a package-free
base-R graphics control: four unsuppressed Fontconfig/Pango leak contexts, before
IntegMultiReg was loaded. Graphics audit 20894229 ran nine controls (three each
of base, extended and vignette graphics). Their 66 raw leak records yielded
12 exact, bounded 24-frame font-only rules, with no wildcard, invalid-access,
uninitialized-read or IntegMultiReg-frame rule. These are local environment
exclusions, not evidence that those graphics libraries are leak-free.

Strict retry 20894290 must pass all nine scoped replays, the R-heap and archived
7200-byte leak negative controls, and the complete candidate/worker checks before
NATIVE-VALIDATED can be written. Bootstrap 20887915 now depends on that retry.
The scientific study has not started. No convergence or full-study acceptance
is claimed. Latest IDs are in runtime/current-jobs.tsv; graphics evidence is
under runtime/candidate-b/anvil-graphics-audit/.
