# Anvil implementation, 2026-09-23

## Research script consistency: September 26, evening

`paper/` and `hpc/` had drifted from the convention `IntegMultiReg/R` follows
without exception across 4,568 lines. The variable-shadowing fault repaired on
September 25 lived in the densest of those files, where one statement's end and
the next one's start are not visible to a reader.

Formatting was applied in two passes, each rebuilding a line from its parse
tokens rather than rewriting text, so no separator decision could be confused by
a string literal, a comment, or an expression such as `x < -1`. Every file was
accepted only when `deparse(parse())` was identical before and after, proving the
change purely lexical: 25 files for spacing and semicolon splitting, 13 more for
binary arithmetic. Indentation was widened only where the existing indent matched
the brace depth, and column-aligned continuations were restored where an earlier
attempt had flattened them.

Naming then changed in a separate pass, which is not lexical and is recorded as
such. `key()` had four definitions over three key schemas and `run()` four
incompatible signatures; both are now named for what they key or run. Identifiers
that shadowed `base::all`, `outer`, `rep`, `table` and `summary` — each of which
is genuinely called in sibling scripts — were renamed. `ridge` in
`hpc/test-plan.R` held an argument list rather than the penalty the package
argument names, and `model` in the Table 1 validator held a logical flag.

Three duplications were resolved differently by kind. `read_result()` and one
inlined copy simply repeated `read_experiment_result()` from the helpers that the
same files had already sourced; they are deleted. `auc_pairwise()` and `auc()` are
*deliberate* independent re-derivations of `selection_auc()`, as is the hand-written
aggregation in the Table 1 validator; those are kept, given one name,
`independent_pairwise_auc()`, and commented so they no longer read as accidents.

`first_folds` is now reset for each replicate and asserted before use in both
validators. It was safe only because the method order happened to put IMR first,
which is the same class of fragility as the September 25 fault.

`hpc/fetch-validation-sources.R` had no caller anywhere and is superseded by
`dependencies/` with `hpc/finalize-dependencies.sh`; it is deleted. Two other
apparent orphans are not: `paper/appendix-marker-rankings.R` is invoked by hand
per `paper/REPLICATION-README.md`, and `paper/validate-manuscript-run.R` is hashed
in the frozen manuscript provenance. Both are kept.

Twelve commented-out debug statements were removed from `src/`, including two
calls to a `SAMPLER_DEBUG` macro that no longer exists. The commented
`free(yobs[m])` in `fit.c` was *correct* and stays disabled: `yobs` is allocated
only for non-binary outcomes and is freed with them a few lines below. It now
says so instead of reading as an oversight. All 21 C files still pass a syntax
check.

Nothing prevented the original drift: there is no linter in the tree, the eight
CI workflows cover only correctness and memory safety, and `.Rbuildignore`
excludes `paper/` from `R CMD check`. `hpc/check-script-style.R` now runs first in
`hpc/validate-package.sh`. It was calibrated against `IntegMultiReg/R`, which must
pass by definition; that identified two rules that would have been wrong, and both
are encoded — `x[i, , drop = FALSE]` keeps its space, and a semicolon separating a
short paired assignment is allowed. Binary arithmetic is deliberately *not*
gated: the package sources space it in 142 of 151 uses, which is a convention but
not a rule. `tools/test-study-overlap.py`, whose fifteen tests cover the
controller's lost-submission, throttle and acceptance-forgery paths, had no runner
and now runs in the same script.

Behaviour was checked by the eight installed test scripts, the fifteen controller
tests and `hpc/test-splitting.R`, which executes real short experiments through
the rewritten validators. `paper/expected-results-0.2.0-paper/provenance.json`
records the previous and current hashes of the two pinned scripts that changed,
states that only spacing and local identifier names differ, and states explicitly
that the frozen results were **not** regenerated.

The live `research-source-hashes.rds` under `runtime/candidate-b/evidence/` no
longer matches `paper/` and `hpc/`. That file gates `hpc/prepare-study.R`, so it
must be recaptured with `hpc/check-research-source.R` before a new study starts.
The running study is unaffected: it executes the frozen copies under
`runs/study-20260923/source/`.

The work tree also gained its first off-machine copy, pushed to branch
`anvil/implementation-20260923`. The `.gitignore` allowlist had excluded the
manuscript itself, so `paper/*.tex`, the bibliography, the class files, the frozen
expected-results record, the historical audit snapshots and the pre-Anvil runner
copies are now tracked. `paper/data` stays untracked: the Biometrics/Wiley
supplement has a different license and distribution boundary from the reduced
public `kircIMR` example, and that repository is public.

## Controller watchdog: September 26, 18:51 local time

The detached controller on login01 is the only component that submits the final
`full` validation, and `overlap-state.json` still records `"validator": null`.
A temporary user crontab on login01 now runs the controller's existing
`--once` entry point every ten minutes, appending to
`runs/study-20260923/logs/watchdog-cron.log`.

No code was added. `--once` takes the same `.highmem-routing.lock`, so it exits
with `BlockingIOError` and writes nothing while the controller is alive. A
manual invocation at 18:47 confirmed that behaviour: the process exited 1, and
`overlap-state.json` kept its `2026-09-25T18:58:58Z` timestamp unchanged. The
entry point is already covered by `execution_source_hashes`, so scheduling it
perturbs no gate. Recovery procedure and the verbatim relaunch command are now
in `hpc/README.md`.

This crontab lives in login01's local spool, so it does not survive a reboot of
that node, and it is a stopgap rather than a supported scheduler. **Remove it
with `crontab -r` once `runs/study-20260923/VALIDATION-COMPLETED` exists.**

Hosting the controller inside a Slurm job was rejected: `sbatch --test-only` for
a 1-CPU, 96-hour `shared` job returned an estimated start of 2026-12-22, which
trades the login-node dependency for a worse queueing one.

At this checkpoint 531 of the 679 tasks have `TASK-COMPLETED`, the array holds
its 64-task throttle, and `blocked` and `submission_in_progress` are both null.

## Collection coverage correction: September 25, 20:55 local time

Review of the final collector found a validator variable-shadowing error:
`selection` held the requested configuration/replicate list, then was replaced
by a fitted selection-AUC table inside the loop. A combined run could therefore
validate only its first configuration. The prior six-task smoke fixture's audit
had six rows instead of the expected eighteen. Individual pilot validations each
cover one configuration, so their checks were not skipped by this issue.

The validator now uses separate `task_selection` and `selection_summary` names
and verifies complete configuration/replicate/method coverage before acceptance.
Regression job `20909668` completed in 21 seconds: it reproduced the old missed
failure, checked all eighteen expected entries with the fix, and rejected a
failed result in configuration 3 / replicate 2. Evidence is retained under
`runtime/validation-coverage-1/`. The normal serial/split regression harness now
also includes this coverage assertion and late-configuration negative control.

Corrected validation-only snapshots and an explicit collector override were
installed under the active study's `execution/` directory, with updated hashes
and archived prior policy. The scientific source tree, study manifest, running
samplers and results were not changed. Frozen scientific hashes were reverified;
prior research-source acceptance was archived before accepting the tested
postprocessing updates. Collector dispatch checks and fourteen scheduling tests
passed. The detached execution controller is healthy after the handoff.

The scheduled regression finished before a separately attempted fallback step
started. The fallback refused to overwrite its completed fixture and exited in
one second; that expected guard log is retained in EXECUTION-NOTE.txt. The
scientific batch continued, and no sampler ran in the fallback step.

The current production snapshot is 105 completed, 58 running and 486 queued.
Progress reports continue every fifteen minutes with the 64-task limit.
Final collection/validation now has a 48-hour reservation (within highmem's
limit) to cover all 660 simulation results and eight chains; the pilot check
keeps its 12-hour limit. All fifteen scheduling tests passed after this change.

## Production checkpoint: September 25, 19:50 local time

All ten audits and all twenty full-parameter pilot tasks have completed and
passed their respective acceptance gates. Pilot validator `20906255` completed
successfully in 13m15s, including prediction/score replay, baseline checks and
the original-scale Table 1 audit. Its measured resource recommendation is
12288 MiB; the already submitted production tasks retain their conservative
16384 MiB allocation.

Production array `20905186` is now at its full throttle of 64. The refreshed
snapshot records 78 completed, 46 running and 525 queued production tasks.
No task failures, resource retries or execution blockers are recorded. The
code2017 main simulation's first configuration has all 50 repetitions available;
its second configuration has 30. These task results await complete-study
collection and scientific interpretation.

The detached controller is alive. A few nonzero Slurm responses while changing
an array containing finished elements were handled by the existing refresh/retry
path; the live throttle and controller state now both confirm 64. No scientific
job was restarted for those metadata responses. The current progress CSV has
1465 score rows from 97 completed scored tasks, plus the completed diagnostic
pilot chain, which is not a score-table task.

## Concurrency increased: September 25, 12:19 local time

The user requested another increase. The active execution policy now caps the
study at 64 tasks: production array `20905186` has throttle 63 and the remaining
pilot array `20901042` has throttle 1. Production rises to 64 after pilot
acceptance. The prior cap-32 policy was archived under `execution/policy-history/`.
The detached controller was restarted with the amended policy and both live
Slurm throttles were verified. Scientific job IDs and settings were retained.

At the unchanged nine allocated CPUs per initial task, 64 tasks request at most
576 CPUs, below shared's 2048-CPU per-user limit. Highmem still permits only two
running jobs per user. Actual execution at this checkpoint is one diagnostic
pilot and one production chain; other production tasks await scheduling.
The nineteenth pilot task has completed. Raising the array cap does not itself
remove the current partition/priority constraints.

## Production submitted: September 25

At the user's explicit request, the remaining 649 tasks have been submitted and
released as array `20905186`. They are 642 simulation repetitions plus seven
diagnostic chains; the initial 20 scientific tasks count toward the same study.
Production no longer waits for the two unfinished pilot tasks before computing.
All scientific methods, seeds, draw counts, data and source hashes are unchanged.

The total study concurrency limit is now 32. The existing pilot array
`20901042` is capped at two slots, and production starts with 30; after genuine
pilot acceptance, production rises to 32. The initial 16 GiB / nine allocated
CPU / 48-hour requests are retained conservatively for this early release.
Highmem's site limits (two running/four submitted per user) still apply.

Obsolete controller `20901043` was cancelled so it cannot submit the full array
a second time. A detached, bounded concurrent controller now owns separate
pilot and full ledgers, validation-only Slurm jobs, and once-only resource
retries. Numerical/cancelled failures stop new admissions for review; final
collection requires pilot acceptance and all scientific tasks. No validation
marker was fabricated to permit the concurrent submission.

Execution policy/state and hashed controller copies are under
`runs/study-20260923/execution/`; the supervisor PID and log are in
`runtime/overlap-controller.json`. The scientific manifest remains immutable.
Fourteen scheduling tests passed, covering the shared concurrency budget,
phase-specific retry resources/reservations, loss-of-response protection,
validation evidence, collection gates and highmem admission priority. The
validation job request also passed Slurm's test-only check. Progress still
refreshes every fifteen minutes without routine chat notifications.

## Morning checkpoint: September 25, 11:22 local time

All ten strict-method audit tasks completed and `AUDITS-ACCEPTED` was written.
Two method runs retained independently verified zero-ridge rank failures:
paper simulation configuration 2 IMR, and paper Table 1 clinical/molecular IMR.
They remain audit evidence rather than successful unpenalized CV results.

Pilot array `20901042` has completed 18 of 20 full-parameter tasks: all nine
code2017 simulations, eight paper/ridge=.001 simulations, and the original-scale
Table 1 rerun. The million-draw chain (task 1) and paper main configuration 3
(task 439) are still actively computing. Task 439 has saved its fit checkpoint.
No pilot task has an abnormal Slurm exit. Full production's remaining 649 tasks
have not been submitted; they await pilot acceptance.

Completed scientific tasks have 280 score rows in
`runs/study-20260923/progress/scientific-scores.csv`. Table 1 has 100/100 valid
validation scores for every reported method/subgroup. These are current-run
results, not claims of historical digit agreement or completed scientific review.
There are 228 saved baseline warnings across six completed pilot tasks: 200
clinical Cox warnings, 25 glmnet numerical-path warnings and three univariate
Cox warnings. Their messages are preserved in
`runs/study-20260923/progress/pilot-warnings-20260925.csv`; the scheduled validators
will replay results before full release. Warnings are not silently discarded.

Pilot acceptance controller `20901043` now requests twelve hours and can use
shared/highmem, while retaining its dependency on the entire pilot array.
The detached monitor remains alive, follows phase changes, and refreshes the
live [progress report](runs/study-20260923/PROGRESS.md) every fifteen minutes.
Source, methods, seeds, draw counts and active scientific jobs are unchanged.

## Current checkpoint: September 24, after GitHub validation

Candidate B has passed the installed regression suite, ASAN/UBSAN, all nine
new GitHub CI jobs, and final Anvil strict check `20899640`. The latter completed
in 3h06m02s with exit code zero; native acceptance was written at
2026-09-25 02:59:41 UTC (September 24 local time).
At the user's request, bootstrap `20900072` now starts audits and the pilot
concurrently with that check, through the explicit `--pending-native-job` option.
The former waiting bootstrap `20887915` was cancelled. All three candidate
validation markers now exist; full production and collection still require
their scientific phase gates. Use Slurm
and `runtime/current-jobs.tsv` for live status. Earlier progress entries below
retain their original job IDs and observations.

Bootstrap `20900072` completed successfully and froze all 679 tasks into
`runs/study-20260923`. Audit array `20900090` contains tasks 670–679; task 670
has begun IMR fitting on highmem alongside native check `20899640`. Tasks 670
and 671 can use shared/highmem; the remaining audit tasks use shared. Controller
`20900091` will accept the audit and submit the 20 full-parameter pilot tasks.
The pilot's acceptance/full-release controller waits for native-check success.
The frozen source/runtime manifest and the live method status were verified.

At the September 24 22:20 check, audit tasks 670 and 671 had completed in
14m06s and 18m43s. Task 670 completed both methods; task 671 retained the
expected zero-ridge IMR rank failure and completed BMS. Eight audits remained
queued, so tasks 672–674 were given shared/highmem eligibility. Task 672 has
started. Routing job `20900652` also remained queued and was cancelled. A
lightweight metadata-only process now follows arrays recorded in this study,
adding highmem eligibility within its four-submitted-job limit. It changes
neither task resources nor scientific settings and excludes tasks exceeding
highmem's 48-hour limit. It does not route unrelated jobs, retries transient
Slurm errors, and is bounded to 24 hours. Its process/session/log details are
in `runtime/partition-router.json`; routing events are saved in the study.
Audit acceptance controller `20900091` now requests 30 minutes and can use
shared/debug/highmem; it only reads the ten audit results and submits the pilot.

The strict package check passed examples, tests and vignette rebuilding with
zero errors/warnings and two reviewed NOTEs (instrumented example timing and
unavailable local HTML tools, covered by CI). The additional unrestricted
suite passed 2850 assertions without failures/warnings/skips, and verified 176
instrumented workers with zero Memcheck errors and definite leaks. All twelve
graphics controls and both negative controls passed their acceptance checks.
Hashes and a complete summary are in
`runtime/candidate-b/evidence/native-validation-summary.json`.
At September 24 23:35 local time, audit tasks 670–678 had completed and task
679 was running. The 20-task pilot
and full production remain pending scientific audit acceptance. The manifest's
`native_validated_at_start = FALSE` is retained as an accurate historical record.

The overnight monitor is now detached from the chat's launching shell with
its own process session, a 24-hour timeout and durable redirected logs. It
polls routing metadata every 60 seconds and silently refreshes
`runs/study-20260923/PROGRESS.md` every 15 minutes. The report records phase
counts, outstanding controller failures and links to completed-task score CSVs.
Audit and scientific scores are separate; large RDS fits are never loaded by
this monitor, and it never creates acceptance markers. The current supervisor
PID, worker PID, log and schedule are in `runtime/partition-router.json`.
Scientific arrays and their acceptance controllers remain Slurm jobs.

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

Attempt `20894476` used `runtime/candidate-b/strict-anvil-2`: examples and
2804 assertions passed under R CMD check's process limit, and the worker-log
verifier accepted all 110 actual PSOCK workers with zero Memcheck errors or
definite leaks. Vignette rebuilding exposed two additional Fontconfig stacks
through strwidth(). These results are recorded in
`strict-anvil-2/strict-evidence/partial-validation.json`; they do not constitute
native acceptance. The standalone suite with its additional worker counts
remains required. Anvil also lacks tidy/V8 for HTML checks; the completed
cross-platform CI covers the manual checks separately.

Graphics job `20899613` reproduced both unresolved 24-frame stacks in each
of three package-free width-first rendering runs. The combined audit under
`runtime/candidate-b/anvil-graphics-audit-2` retains the original twelve rules
and adds three bounded rules: the two reproduced definite-leak contexts and
one language-cache context from the new control. There are still no wildcard,
invalid-access, uninitialized-read or package-frame rules.

New attempt `20899640` uses `strict-anvil-3`, replays all twelve graphics
controls, reruns both negative controls and the full package/standalone worker
checks, and writes `NATIVE-VALIDATED` only on success. The dependency on this
check now applies to the controller that accepts the pilot and releases full
production, allowing audit/pilot computation to proceed first.
Tested research-source hashes, the frozen package
source manifest, and all 32 dependency versions/load paths were reverified.
After validation, the bootstrap freezes `runs/study-20260923` and schedules
audit → full-parameter pilot → full study → collection with concurrency 16.

The concurrent-launch change passed five orchestration checks: explicit
preparation/default refusal, native and phase gates, real submitter/controller
logic with mock Slurm, bounded retries, and all 679 task settings. The source
audit proved that every prior non-orchestration input and the frozen package
source were unchanged. Original and refreshed research-source acceptance are
both retained under `runtime/candidate-b/evidence/concurrent-orchestration/`.
Package validation markers were not changed. The frozen manifest records native
validation as pending, and actual acceptance is required before full production.

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
