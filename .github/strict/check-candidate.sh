#!/usr/bin/env bash
set -euo pipefail
evidence="${RUNNER_TEMP}/strict-evidence"
mkdir -p "${evidence}/candidate" "${evidence}/check-workers"
export IMR_VALGRIND_LOG_DIR="${evidence}/check-workers"
chmod +x .github/scripts/valgrind-worker-rscript.sh
# R's check runner copies this runtime startup file into tests/startup.Rs.
# Retain and restore the original; the pinned package artifact is unchanged.
startup="$(R RHOME)/share/R/tests-startup.R"
cp "${startup}" "${evidence}/tests-startup-original.R"
trap 'cp "${evidence}/tests-startup-original.R" "${startup}"' EXIT
printf '\nsource(file.path(Sys.getenv("GITHUB_WORKSPACE"), ".github", "scripts", "configure-valgrind-workers.R"))\n' >> "${startup}"
cp "${startup}" "${evidence}/tests-startup-instrumented.R"
set +e
R CMD check --as-cran --use-valgrind --run-donttest \
  --output="${evidence}/candidate" .github/strict/fixtures/IntegMultiReg_0.2.0.tar.gz \
  2>&1 | tee "${evidence}/candidate-console.log"
status=${PIPESTATUS[0]}
set -e
echo "${status}" > "${evidence}/candidate-exit-status.txt"
python3 .github/strict/verify-check.py "${evidence}"
exit "${status}"
