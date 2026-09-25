#!/bin/bash -l
# Validation only: the execution controller owns submissions and bounded retries.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" && $# == 2 ]] || exit 2
study="$(realpath "$1")"
phase="$2"
source "$study/source/hpc/environment.sh"
case "$phase" in
  pilot) Rscript --vanilla "$study/source/hpc/accept-phase.R" "$study" pilot ;;
  full)
    [[ -f "$study/PILOT-ACCEPTED" ]]
    Rscript --vanilla "$study/source/hpc/collect.R" "$study" ;;
  *) exit 2 ;;
esac
