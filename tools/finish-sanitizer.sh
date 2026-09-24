#!/bin/bash -l
# Complete the failed harness stage against the already-tested sanitizer build.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
candidate="$(realpath "$1")"
work="$candidate/sanitizer"
Rscript --vanilla - "$work" <<'RS'
work <- commandArgs(TRUE)[1]
tests <- read.csv(file.path(work, "evidence/test-evidence/tests.csv"))
stopifnot(sum(tests$passed) == 2850L, !any(tests$failed), !any(tests$error),
          !any(tests$warning), !any(tests$skipped))
stopifnot(any(grepl("Native Gamma density, prior indexing and stationary-distribution checks passed.",
                    readLines(file.path(work, "sampler.log")), fixed = TRUE)))
RS
export R_MAKEVARS_USER="$work/Makevars"
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
export SANITIZER_RUNTIME="$(gcc -print-file-name=libasan.so):$(gcc -print-file-name=libubsan.so)"
export LD_PRELOAD="$SANITIZER_RUNTIME"
export R_LIBS_USER="$work/library" R_LIBS="$work/library:$IMR_DEPENDENCIES"
export GITHUB_WORKSPACE="$work/source" RUNNER_TEMP="$work/fault-retry-1"
[[ ! -e "$RUNNER_TEMP" ]]
mkdir -p "$RUNNER_TEMP"
cp "$IMR_ROOT/IntegMultiReg/.github/scripts/verify-cv-row-failures.sh" "$RUNNER_TEMP/harness.sh"
sha256sum "$RUNNER_TEMP/harness.sh" > "$RUNNER_TEMP/harness.sha256"
bash "$RUNNER_TEMP/harness.sh" > "$RUNNER_TEMP/cleanup.log" 2>&1
if grep -E 'runtime error:|ERROR: AddressSanitizer|SUMMARY: (AddressSanitizer|UndefinedBehaviorSanitizer)' \
  "$work/install.log" "$work/tests.log" "$work/sampler.log" "$RUNNER_TEMP/cleanup.log"; then
  exit 1
fi
date -u +%FT%TZ > "$candidate/SANITIZER-VALIDATED"
