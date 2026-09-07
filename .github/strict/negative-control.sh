#!/usr/bin/env bash
set -euo pipefail
evidence="${RUNNER_TEMP}/strict-evidence"
oldlib="${RUNNER_TEMP}/old-library"
mkdir -p "${oldlib}"
# Deliberately faulty audit-only code. Level 2 must expose an uninitialized
# R-heap read; this is separate from reproducing the historical package report.
cp .github/strict/uninitialized-probe.c "${evidence}/"
(cd "${evidence}" && R CMD SHLIB uninitialized-probe.c) > "${evidence}/probe-build.log" 2>&1
set +e
R -d "valgrind ${VALGRIND_OPTS}" --vanilla -e \
  "dyn.load('${evidence}/uninitialized-probe.so'); .Call('imr_uninitialized_probe')" \
  > "${evidence}/uninitialized-probe.log" 2>&1
probe_status=$?
set -e
echo "${probe_status}" > "${evidence}/probe-exit-status.txt"
test "${probe_status}" -eq 97
grep -q 'Conditional jump or move depends on uninitialised value' "${evidence}/uninitialized-probe.log"
grep -q 'imr_uninitialized_probe' "${evidence}/uninitialized-probe.log"
R CMD INSTALL --preclean --install-tests --library="${oldlib}" \
  .github/strict/fixtures/IntegMultiReg_0.1.1.tar.gz > "${evidence}/old-install.log" 2>&1
cd "${oldlib}/IntegMultiReg/tests"
set +e
R_LIBS_USER="${oldlib}:${R_LIBS_USER}" R -d "valgrind ${VALGRIND_OPTS}" \
  --vanilla -f testthat.R > "${evidence}/old-tests.log" 2>&1
status=$?
set -e
echo "${status}" > "${evidence}/old-exit-status.txt"
# A crash or unrelated test failure is not sufficient evidence of sensitivity.
test "${status}" -eq 97
grep -q 'main_function_prediction' "${evidence}/old-tests.log"
grep -q 'predict_cv_fold' "${evidence}/old-tests.log"
grep -q 'definitely lost: 7,200 bytes in 60 blocks' "${evidence}/old-tests.log"
echo 'R-heap uninitialized-read probe detected; old package reproduced its exact 7,200-byte/60-block leak.' | tee "${evidence}/negative-control-result.txt"
if grep -q 'Conditional jump or move depends on uninitialised value' "${evidence}/old-tests.log"; then
  echo 'Old-package uninitialized read: reproduced.' | tee -a "${evidence}/negative-control-result.txt"
else
  echo 'LIMITATION: old-package uninitialized read NOT reproduced in this compiler/environment; the probe is not a reproduction of that historical call stack.' | tee -a "${evidence}/negative-control-result.txt"
fi
