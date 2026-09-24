#!/bin/bash -l
# Collect package-free controls before defining any Anvil-specific exclusion.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
candidate="$(realpath "$1")"
native="$IMR_NATIVE_RUNTIME"
out="$candidate/anvil-graphics-audit"
[[ ! -e "$out" && -f "$native/READY" ]]
mkdir -p "$out"
source "$native/build-environment.sh"
export PATH="$native/strict-runtime/bin:$IMR_ROOT/runtime/tools/pandoc/bin:$PATH"
export R_LIBS_USER="$native/validation-library" R_LIBS="$native/validation-library"
export CLI_NO_THREAD=1
unset LD_PRELOAD
cd "$candidate/source"
for name in graphics-control graphics-extended-control graphics-vignette-control; do
  for repetition in 1 2 3; do
    log="$out/$name-$repetition.log"
    set +e
    R -d "valgrind --tool=memcheck --leak-check=full --show-leak-kinds=definite,possible --errors-for-leak-kinds=definite,possible --num-callers=40 --track-origins=yes --gen-suppressions=all --error-exitcode=97" \
      --vanilla -f ".github/strict/$name.R" > "$log" 2>&1
    status=$?
    set -e
    printf '%s\n' "$status" > "$log.exitcode"
    [[ "$status" == 0 || "$status" == 97 ]] || exit "$status"
  done
done
sha256sum "$out/"*.log > "$out/SHA256SUMS"
date -u +%FT%TZ > "$out/COLLECTED"
