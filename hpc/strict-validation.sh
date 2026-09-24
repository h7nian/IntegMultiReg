#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
module load texlive/20200406
candidate="$(realpath "$1")"
native="$IMR_NATIVE_RUNTIME"
[[ -f "$native/READY" && -f "$candidate/UNIT-VALIDATED" && -f "$candidate/SANITIZER-VALIDATED" ]]
work="$candidate/strict"
[[ ! -e "$work" ]]
mkdir -p "$work/strict-evidence" "$work/library"
cp -a "$candidate/source" "$work/source"
source "$native/build-environment.sh"
export PATH="$native/strict-runtime/bin:$IMR_ROOT/runtime/tools/pandoc/bin:$PATH"
export RUNNER_TEMP="$work" GITHUB_WORKSPACE="$work/source"
export R_LIBS_USER="$work/library" R_LIBS="$work/library:$native/validation-library"
export CLI_NO_THREAD=1 R_KEEP_PKG_SOURCE=yes R_TEXI2DVICMD=emulation
export _R_CHECK_ALWAYS_LOG_VIGNETTE_OUTPUT_=true
export _R_CHECK_CRAN_INCOMING_REMOTE_=false
unset LD_PRELOAD
cp "$native/strict-runtime/provenance/"* "$work/strict-evidence/"
cp "$native/source-identity-audit.json" "$work/strict-evidence/"
cd "$work"
R CMD build source > build.log 2>&1
cp IntegMultiReg_0.2.0.tar.gz source/.github/strict/fixtures/IntegMultiReg_0.2.0.tar.gz
sha256sum IntegMultiReg_0.2.0.tar.gz > CANDIDATE-SHA256SUMS
cd source
export VALGRIND_OPTS="--tool=memcheck --leak-check=full --show-leak-kinds=definite,indirect,possible --errors-for-leak-kinds=definite,possible --num-callers=40 --track-origins=yes --error-exitcode=97"
for file in pango.supp glib.supp ubuntu-font-baseline.supp ubuntu-vignette-font-baseline.supp; do
  export VALGRIND_OPTS="$VALGRIND_OPTS --suppressions=$GITHUB_WORKSPACE/.github/strict/upstream/$file"
done
# Keep the existing scoped third-party rules and prove the controls here.
for name in graphics-control graphics-extended-control graphics-vignette-control; do
  R -d "valgrind $VALGRIND_OPTS" --vanilla -f ".github/strict/$name.R" > "$work/strict-evidence/$name.log" 2>&1
  grep -q 'ERROR SUMMARY: 0 errors from 0 contexts' "$work/strict-evidence/$name.log"
done
bash .github/strict/negative-control.sh
bash .github/strict/check-candidate.sh
mkdir -p "$work/candidate-library" "$work/strict-evidence/workers"
R CMD INSTALL --preclean --install-tests --library="$work/candidate-library" \
  .github/strict/fixtures/IntegMultiReg_0.2.0.tar.gz > "$work/strict-evidence/worker-install.log" 2>&1
export IMR_VALGRIND_LOG_DIR="$work/strict-evidence/workers"
chmod u+x .github/scripts/valgrind-worker-rscript.sh
R_LIBS="$work/candidate-library:$native/validation-library" \
  R -d "valgrind $VALGRIND_OPTS" --vanilla -f .github/scripts/verify-valgrind-suite.R \
  > "$work/strict-evidence/worker-suite.log" 2>&1
python3 .github/strict/verify-candidate-workers.py "$work/strict-evidence/worker-suite.log"
date -u +%FT%TZ > "$candidate/NATIVE-VALIDATED"
