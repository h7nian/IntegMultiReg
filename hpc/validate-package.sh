#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
candidate="$(realpath "$1")"
[[ -f "$candidate/SOURCE-SHA256SUMS" && ! -e "$candidate/library" ]]
cd "$candidate"
sha256sum --check --quiet SOURCE-SHA256SUMS
mkdir library baseline evidence
export R_LIBS_USER="$candidate/library" R_LIBS="$candidate/library:$IMR_DEPENDENCIES"
export RUNNER_TEMP="$candidate/evidence"
Rscript --vanilla "$IMR_ROOT/hpc/check-research-source.R" capture "$candidate/evidence/research-source-hashes.rds"
R CMD INSTALL --preclean --install-tests --library="$candidate/library" source
R CMD INSTALL --library="$candidate/baseline" "$IMR_ROOT/IntegMultiReg_0.2.0.tar.gz"
Rscript --vanilla source/.github/scripts/verify-tests.R
Rscript --vanilla "$IMR_ROOT/hpc/compare-defaults.R" "$candidate/baseline" "$candidate/evidence/baseline.rds" capture
Rscript --vanilla "$IMR_ROOT/hpc/compare-defaults.R" "$candidate/library" "$candidate/evidence/baseline.rds" compare
cd source
Rscript --vanilla .github/scripts/verify-sampler.R "$candidate/evidence/sampler"
Rscript --vanilla .github/scripts/verify-original-code.R "$IMR_ROOT/Ref/biom12587-sup-0002-suppdata_code.zip" "$candidate/evidence/original-code"
cd "$IMR_ROOT"
Rscript --vanilla paper/test-experiment-cli.R
Rscript --vanilla paper/test-experiment-jobs.R
Rscript --vanilla paper/test-experiment-checkpoint.R
Rscript --vanilla hpc/test-plan.R
Rscript --vanilla hpc/test-phase-gates.R
Rscript --vanilla hpc/test-submit.R
Rscript --vanilla hpc/test-resource-retry.R
Rscript --vanilla hpc/test-splitting.R "$IMR_ROOT" "$candidate/evidence/task-splitting"
Rscript --vanilla "$IMR_ROOT/hpc/check-research-source.R" verify "$candidate/evidence/research-source-hashes.rds"
date -u +%FT%TZ > "$candidate/UNIT-VALIDATED"
