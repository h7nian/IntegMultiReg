#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
source "$PWD/hpc/environment.sh"
module load curl/7.76.1
export RUNNER_TEMP="$IMR_ROOT/runtime/native-runtime-1"
[[ ! -e "$RUNNER_TEMP" ]]
mkdir -p "$RUNNER_TEMP"
python3 "$IMR_ROOT/hpc/native-build-environment.py" > "$RUNNER_TEMP/build-environment.sh"
source "$RUNNER_TEMP/build-environment.sh"
pkg-config --modversion libcurl libpcre2-8 liblzma cairo pangocairo > "$RUNNER_TEMP/system-versions.txt"
bash "$IMR_ROOT/IntegMultiReg/.github/strict/build-candidate-runtime.sh"
prefix="$RUNNER_TEMP/strict-runtime"
export PATH="$prefix/bin:$PATH"
export R_LIBS_USER="$RUNNER_TEMP/validation-library" R_LIBS="$RUNNER_TEMP/validation-library"
mkdir -p "$R_LIBS_USER"
while IFS=$'\t' read -r package version archive checksum; do
  [[ "$package" == Package ]] && continue
  printf '%s  %s\n' "$checksum" "$IMR_ROOT/runtime/validation-sources/$archive" | md5sum --check --quiet
  R CMD INSTALL --library="$R_LIBS_USER" "$IMR_ROOT/runtime/validation-sources/$archive"
done < "$IMR_ROOT/runtime/validation-sources/dependencies.tsv"
Rscript --vanilla -e 'library(testthat); library(knitr); library(rmarkdown); library(survival); sessionInfo()' > "$RUNNER_TEMP/validation-session.txt"
date -u +%FT%TZ > "$RUNNER_TEMP/READY"
