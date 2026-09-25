# Anvil execution

The active study now uses the user-authorized concurrent execution policy in
`runs/study-20260923/execution/overlap-policy.json`. Its remaining 649 tasks are
array `20905186`: 642 simulation repetitions and seven diagnostic chains.
Production has 30 slots while the two remaining pilot tasks run, then 32 after
pilot acceptance. The old pilot controller was cancelled to prevent a duplicate
production submission. `manage-study-overlap.py` owns phase-specific retries
and validation-only jobs; its execution sources are copied and hashed alongside
the policy. The original scientific source and manifest are unchanged.
Do not run the serial `advance.sh` concurrently with this controller.

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
