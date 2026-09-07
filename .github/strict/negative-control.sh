#!/usr/bin/env bash
set -euo pipefail
evidence="${RUNNER_TEMP}/strict-evidence"
oldlib="${RUNNER_TEMP}/old-library"
mkdir -p "${oldlib}"
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
grep -q 'Conditional jump or move depends on uninitialised value' "${evidence}/old-tests.log"
grep -q 'main_function_prediction' "${evidence}/old-tests.log"
grep -q 'predict_cv_fold' "${evidence}/old-tests.log"
grep -Eq 'definitely lost: [1-9][0-9,]* bytes' "${evidence}/old-tests.log"
echo 'Old release reproduced uninitialized-memory and definite-leak failures.' | tee "${evidence}/negative-control-result.txt"
