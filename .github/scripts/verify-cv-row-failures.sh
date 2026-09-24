#!/usr/bin/env bash
# Keep injected source and its library separate from the normally tested build.
set -euo pipefail
: "${GITHUB_WORKSPACE:?The tested package source directory is required}"
: "${RUNNER_TEMP:?An isolated runner output directory is required}"
: "${SANITIZER_RUNTIME:?The matching sanitizer runtime is required}"

fault_source="${RUNNER_TEMP}/row-fault-source"
fault_library="${RUNNER_TEMP}/row-fault-library"
mkdir "${fault_source}" "${fault_library}"
test -f "${GITHUB_WORKSPACE}/DESCRIPTION"
# Use the tested working source, including uncommitted changes. The package
# may be nested in a research repository or extracted without Git metadata.
tar -C "${GITHUB_WORKSPACE}" --exclude=.git --exclude='*.o' --exclude='*.so' \
  --exclude='*.dll' --exclude='*.dylib' --exclude=src/Makevars -cf - . | \
  tar -xf - -C "${fault_source}"
test -f "${fault_source}/DESCRIPTION"
# Isolate git apply from any ancestor repository so it cannot skip the patch
# as outside the current working-directory prefix.
git -C "${fault_source}" init -q
git -C "${fault_source}" apply --check \
  "${GITHUB_WORKSPACE}/.github/scripts/cv-row-fault-injection.patch"
git -C "${fault_source}" apply \
  "${GITHUB_WORKSPACE}/.github/scripts/cv-row-fault-injection.patch"

install_log="${RUNNER_TEMP}/sanitizer-row-fault-install.log"
test_log="${RUNNER_TEMP}/sanitizer-row-fault-tests.log"
pattern='runtime error:|ERROR: AddressSanitizer|SUMMARY: (AddressSanitizer|UndefinedBehaviorSanitizer)|UndefinedBehaviorSanitizer:'
if [[ "$(uname -s)" != "Darwin" ]]; then
  export LD_PRELOAD="${SANITIZER_RUNTIME}"
fi
R CMD INSTALL --preclean --no-test-load --library="${fault_library}" \
  "${fault_source}" 2>&1 | tee "${install_log}"
if grep -Eq "${pattern}" "${install_log}"; then
  echo "::error::Sanitizer error while building the isolated fault library"
  exit 1
fi

test_script="${GITHUB_WORKSPACE}/.github/scripts/verify-cv-row-failures.R"
if [[ "$(uname -s)" == "Darwin" ]]; then
  r_home="$(R RHOME)"
  # Match the ordinary macOS sanitizer parent/PSOCK launch mechanism.
  chmod +x "${GITHUB_WORKSPACE}/.github/scripts/macos-sanitizer-rscript.sh"
  R_HOME="${r_home}" DYLD_INSERT_LIBRARIES="${SANITIZER_RUNTIME}" \
    R_LIBS="${fault_library}" "${r_home}/bin/exec/R" --vanilla --no-save \
    -f "${test_script}" 2>&1 | tee "${test_log}"
else
  R_LIBS="${fault_library}" Rscript "${test_script}" 2>&1 | tee "${test_log}"
fi
if grep -Eq "${pattern}" "${test_log}"; then
  echo "::error::Sanitizer error during prediction allocation failure tests"
  exit 1
fi
