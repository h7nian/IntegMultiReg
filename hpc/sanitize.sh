#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
candidate="$(realpath "$1")"
work="$candidate/sanitizer"
[[ ! -e "$work" ]]
mkdir -p "$work/library" "$work/evidence"
cp -a "$candidate/source" "$work/source"
cat > "$work/Makevars" <<'EOF'
CC = gcc
CFLAGS = -g -O1 -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer
LDFLAGS = -fsanitize=address,undefined
EOF
export R_MAKEVARS_USER="$work/Makevars"
export ASAN_OPTIONS=detect_leaks=0:halt_on_error=1
export UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1
export SANITIZER_RUNTIME="$(gcc -print-file-name=libasan.so):$(gcc -print-file-name=libubsan.so)"
export LD_PRELOAD="$SANITIZER_RUNTIME"
export R_LIBS_USER="$work/library" R_LIBS="$work/library:$IMR_DEPENDENCIES"
export RUNNER_TEMP="$work/evidence" GITHUB_WORKSPACE="$work/source"
cd "$work/source"
R CMD INSTALL --preclean --install-tests --library="$work/library" . > "$work/install.log" 2>&1
Rscript --vanilla .github/scripts/verify-tests.R > "$work/tests.log" 2>&1
Rscript --vanilla .github/scripts/verify-sampler.R "$work/evidence/sampler" > "$work/sampler.log" 2>&1
bash .github/scripts/verify-cv-row-failures.sh > "$work/cleanup.log" 2>&1
if grep -E 'runtime error:|ERROR: AddressSanitizer|SUMMARY: (AddressSanitizer|UndefinedBehaviorSanitizer)' "$work/"*.log; then
  exit 1
fi
date -u +%FT%TZ > "$candidate/SANITIZER-VALIDATED"
