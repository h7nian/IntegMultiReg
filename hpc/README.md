# Anvil execution

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

After the candidate passes both gates, run `prepare-study.R CANDIDATE STUDY`
in the configured R environment. It copies research sources and data, hashes
the runtime, and records 679 tasks: ten strict paper audits and 669 scientific
tasks. These contain code2017 and separately labelled paper/ridge=.001 studies
(150 main and 180 correlated tasks each), eight million-draw chains, and an
original-scale code2017 Table 1 repeat.

Submit phases with `bash hpc/submit.sh STUDY audit`, then `pilot`, then `full`.
Only one array can be active for a study, with a global limit of 16 tasks. Run
`accept-phase.R STUDY audit` or `accept-phase.R STUDY pilot` after the preceding
array finishes; these check evidence before opening the next phase. Pilot
tasks retain full-study settings and count toward the final study.

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
