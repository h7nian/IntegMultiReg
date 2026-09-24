#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ ! -f "$ROOT/msi/config.sh" ]]; then
  echo 'Copy msi/config.example.sh to msi/config.sh and fill account/modules first.' >&2; exit 2
fi
source "$ROOT/msi/config.sh"
if [[ "$MSI_ACCOUNT ${MSI_MODULES[*]-}" == *REPLACE_WITH* ]]; then
  echo 'Unfilled MSI account/module placeholders in msi/config.sh' >&2; exit 2
fi
[[ "$MSI_RESULTS" = /* && "$MSI_LIBRARY" = /* ]] || { echo 'Use absolute paths' >&2; exit 2; }
[[ "$MSI_CONCURRENCY" =~ ^[1-9][0-9]*$ ]] || exit 2
export ROOT MSI_RESULTS MSI_LIBRARY
if [[ "${MSI_LOAD_MODULES:-0}" == 1 ]]; then
  module purge
  if [[ -n "${MSI_MODULES[*]-}" ]]; then module load "${MSI_MODULES[@]}"; fi
  command -v Rscript >/dev/null
  command -v gsl-config >/dev/null
  export R_LIBS_USER="$MSI_LIBRARY" R_LIBS="$MSI_LIBRARY"
  export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
  export MAKEFLAGS=-j1
fi
