#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
[[ -d "$IMR_DEPENDENCIES" && ! -f "$IMR_DEPENDENCIES/READY" ]]
Rscript --vanilla "$IMR_ROOT/hpc/check-dependencies.R"
Rscript --vanilla -e 'saveRDS(R.version, file.path(Sys.getenv("IMR_DEPENDENCIES"), "R-version.rds")); sessionInfo()' > "$IMR_DEPENDENCIES/sessionInfo.txt"
gsl-config --version > "$IMR_DEPENDENCIES/gsl-version.txt"
date -u +%FT%TZ > "$IMR_DEPENDENCIES/READY"
