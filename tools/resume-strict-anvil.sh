#!/bin/bash -l
# Recheck the same built package with independently audited Anvil font rules.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
module load texlive/20200406
candidate="$(realpath "$1")"
native="$IMR_NATIVE_RUNTIME"
audit="$(realpath "${3:-$candidate/anvil-graphics-audit}")"
attempt="${2:-strict-anvil-1}"
[[ "$attempt" =~ ^strict-anvil-[1-9][0-9]*$ ]] || exit 2
work="$candidate/$attempt"
[[ -f "$candidate/UNIT-VALIDATED" && -f "$candidate/SANITIZER-VALIDATED" && ! -e "$work" ]]
python3 - "$audit" <<'PY'
from pathlib import Path
import hashlib,json,sys
p=Path(sys.argv[1]); audit=json.loads((p/'rules-audit.json').read_text())
assert hashlib.sha256((p/'anvil-font-baseline.supp').read_bytes()).hexdigest()==audit['rule_sha256']
for name,expected in audit['source_logs'].items():
    assert hashlib.sha256((p/name).read_bytes()).hexdigest()==expected
PY
(cd "$candidate/strict" && sha256sum --check CANDIDATE-SHA256SUMS)
mkdir -p "$work/strict-evidence" "$work/library"
cp -a "$candidate/strict/source" "$work/source"
cp "$audit/anvil-font-baseline.supp" "$audit/rules-audit.json" "$work/strict-evidence/"
cp "$native/strict-runtime/provenance/"* "$work/strict-evidence/"
cp "$candidate/strict/CANDIDATE-SHA256SUMS" "$native/source-identity-audit.json" "$work/strict-evidence/"
source "$native/build-environment.sh"
export PATH="$native/strict-runtime/bin:$IMR_ROOT/runtime/tools/bin:$IMR_ROOT/runtime/tools/pandoc/bin:$PATH"
command -v checkbashisms > "$work/strict-evidence/checkbashisms-path.txt"
checkbashisms --version > "$work/strict-evidence/checkbashisms-version.txt"
cp "$IMR_ROOT/runtime/tools/devscripts-2.25.15-deb13u1/source.json" \
  "$work/strict-evidence/checkbashisms-source.json"
export RUNNER_TEMP="$work" GITHUB_WORKSPACE="$work/source"
export R_LIBS_USER="$work/library" R_LIBS="$work/library:$native/validation-library"
export CLI_NO_THREAD=1 R_KEEP_PKG_SOURCE=yes R_TEXI2DVICMD=emulation
export IMR_GRAPHICS_PROBE=default
export _R_CHECK_ALWAYS_LOG_VIGNETTE_OUTPUT_=true _R_CHECK_CRAN_INCOMING_REMOTE_=false
unset LD_PRELOAD
export VALGRIND_OPTS="--tool=memcheck --leak-check=full --show-leak-kinds=definite,indirect,possible --errors-for-leak-kinds=definite,possible --num-callers=40 --track-origins=yes --error-exitcode=97 --suppressions=$work/strict-evidence/anvil-font-baseline.supp"
cd "$work/source"
for name in graphics-control graphics-extended-control graphics-vignette-control; do
  for repetition in 1 2 3; do
    log="$work/strict-evidence/$name-$repetition.log"
    R -d "valgrind $VALGRIND_OPTS" --vanilla -f ".github/strict/$name.R" > "$log" 2>&1
    grep -q 'ERROR SUMMARY: 0 errors from 0 contexts' "$log"
  done
done
if [[ -f "$audit/graphics-vignette-width-control-1.log" ]]; then
  for repetition in 1 2 3; do
    log="$work/strict-evidence/graphics-vignette-width-control-$repetition.log"
    IMR_GRAPHICS_PROBE=width R -d "valgrind $VALGRIND_OPTS" --vanilla \
      -f .github/strict/graphics-vignette-control.R > "$log" 2>&1
    grep -q 'ERROR SUMMARY: 0 errors from 0 contexts' "$log"
  done
fi
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
