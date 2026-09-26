# Anvil execution

The active study now uses the user-authorized concurrent execution policy in
`runs/study-20260923/execution/overlap-policy.json`. Its remaining 649 tasks are
array `20905186`: 642 simulation repetitions and seven diagnostic chains.
Production now has 63 slots while the last pilot task runs, then 64 after
pilot acceptance. The old pilot controller was cancelled to prevent a duplicate
production submission. `manage-study-overlap.py` owns phase-specific retries
and validation-only jobs; its execution sources are copied and hashed alongside
the policy. The original scientific source and manifest are unchanged.
Do not run the serial `advance.sh` concurrently with this controller.

The active execution snapshot also contains a corrected simulation validator
and `collect-study.R`. The validation-only wrapper supplies the explicit
`--simulation-validator FILE` option to that collector. This fixes complete
configuration coverage without changing any frozen inputs used by running
tasks; the override files and prior policy are hashed and archived.

## Recovering the overlap controller

The controller is a detached login-node process, and it is the only component
that submits the final `full` validation. If it stops, array `20905186` keeps
computing under its Slurm-side `ArrayTaskThrottle`, but highmem routing, the
fifteen-minute progress refresh, bounded resource retries and that final
submission all stop with it. A stale `PROGRESS.md` timestamp is the signal.

`manage-study-overlap.py STUDY --once` performs exactly one poll and exits. It
takes the same `.highmem-routing.lock`, so it fails with `BlockingIOError` and
writes nothing while a controller is alive, and otherwise carries out whatever
that poll is due to do. It is safe to run at any time and safe to schedule.

To restore continuous operation, relaunch with its own process session so it
survives logout, and give it a new log name rather than truncating the record
of the previous run:

```
cd /anvil/scratch/x-szhang30/IntegMultReg/work/anvil-20260923
setsid timeout 604800 python3 -u \
  runs/study-20260923/execution/manage-study-overlap.py runs/study-20260923 \
  < /dev/null > runs/study-20260923/logs/overlap-controller-NAME.log 2>&1 &
```

Never `SIGKILL` a live controller. `retry()` and `validate()` record
`submission_in_progress` before calling `sbatch` and clear it afterwards; a kill
inside that window leaves the flag set, and the next poll then calls `block()`,
which holds both array jobs.

Two recorded states stop a restart. When `blocked` is non-null the process
starts, refreshes the report and exits immediately; no code path clears it, so
set `"blocked": null` in `execution/overlap-state.json`, `scontrol release` the
jobs named in `hold_requests`, and delete `execution/ATTENTION.md`. When
`submission_in_progress` is non-null, first establish from `sacct` whether the
recorded submission actually landed: if it did, append it to
`state[phase]["arrays"]` and to `submissions.tsv` before clearing the flag.
Editing `overlap-state.json` is safe; it is not among the hashed execution
sources. Editing any of those six hashed files is not, and blocks the study.

The maintained runners are in `paper/`; this directory prepares, schedules
and verifies their independent tasks. Every reference is an ordinary set of
package arguments. Explicit overrides are saved alongside the reference name.
Historical MSI entrypoints are archived under `provenance/msi-before-anvil/`.

The configured CPU allocation is `cis260924`, partition `shared`, with modules
GCC 11.2.0, R 4.4.1, GSL 2.4 and Python 3.9.5. Dependencies are installed from the
verified source lock into `runtime/dependencies`. Do not replace that library
while a study is using it. Installation, full tests and experiments use Slurm.

`setup.sh` installs the dependency lock. `tools/freeze-candidate.py NAME` creates
an immutable package-source snapshot; `validate-package.sh CANDIDATE` runs its
installed tests, baseline comparisons and independent sampler checks. Further
instrumented validation is required before `NATIVE-VALIDATED` can be created.

By default, `prepare-study.R CANDIDATE STUDY` requires completed native and
installed-suite validation. Run it in the configured R environment. It copies research sources and data, hashes
the runtime, and records 679 tasks: ten strict paper audits and 669 scientific
tasks. These contain code2017 and separately labelled paper/ridge=.001 studies
(150 main and 180 correlated tasks each), eight million-draw chains, and an
original-scale code2017 Table 1 repeat.

For the default serial policy, submit phases with `bash hpc/submit.sh STUDY audit`, then `pilot`, then `full`.
Only one array can be active for a study, with a global limit of 16 tasks. Run
`accept-phase.R STUDY audit` or `accept-phase.R STUDY pilot` after the preceding
array finishes; these check evidence before opening the next phase. Pilot
tasks retain full-study settings and count toward the final study.

`start-study.sh CANDIDATE STUDY` is the batch entrypoint for unattended execution
after candidate validation. Dependent controllers run the evidence checks and
advance audit → pilot → full → collection. They stop on unexplained failures.
Terminal OOM or timeout failures may receive one resource increase; node failures
or preemption may receive one unchanged retry. Numerical failures and cancelled
jobs are not automatically retried. All retries retain the original task identity
and scientific settings.

To overlap scientific computation with an already running native check, use
`start-study.sh CANDIDATE STUDY --pending-native-job JOB_ID`. This explicit
option requires installed-suite and sanitizer acceptance, verifies that the
native job belongs to the submitting user, and records its ID and pending
status in the frozen manifest. The original no-option launch remains gated.
Audits and the full-parameter pilot can run while the native check continues.
The pilot acceptance controller waits for both its array and native-check
success; full-task submission and final collection independently require the
actual `NATIVE-VALIDATED` marker. A failed native check stops further release.
No acceptance marker is fabricated for this mode.

Each array task has its own directory and lock. Completed tasks are skipped;
failed tasks keep their status, logs and completed fit checkpoints. A fit that
did not finish must restart from its original seed. Verify Slurm's terminal
state before resubmission; a quiet log does not establish that a job stopped.

All scientific tasks use one R process and one BLAS/OpenMP thread. Initial
requests are 16 GiB and 48 hours. Anvil's memory-per-core limit is 1896 MiB, so
the CPU allocation also accounts for memory. Pilot acceptance sets production
memory to 1.5 times the observed peak, rounded up to 2 GiB, with a 4 GiB minimum.
This is a measured planning margin, not a universal bound.

`collect.R STUDY` requires all scientific tasks, checks source/settings identity,
and invokes study and baseline validators before summarization. Strict audit
failures remain separate from ridge-study results. Inspect numerical diagnostics,
Cox warnings, rank Rhat, ESS, MCSE and biomarker stability before changing
manuscript claims; program acceptance alone does not establish convergence.
