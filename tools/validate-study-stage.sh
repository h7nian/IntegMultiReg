#!/bin/bash -l
# Validation only: the execution controller owns submissions and bounded retries.
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" && $# == 2 ]] || exit 2
study="$(realpath "$1")"
phase="$2"
validation_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$study/source/hpc/environment.sh"
case "$phase" in
  pilot) Rscript --vanilla "$study/source/hpc/accept-phase.R" "$study" pilot ;;
  full)
    [[ -f "$study/PILOT-ACCEPTED" ]]
    # These validated postprocessing snapshots are hashed by the controller.
    # Running task sources under study/source remain immutable.
    [[ -f "$validation_root/collect-study.R" && -f "$validation_root/validate-original-simulation.R" ]]
    Rscript --vanilla "$validation_root/collect-study.R" "$study" \
      --simulation-validator "$validation_root/validate-original-simulation.R" ;;
  *) exit 2 ;;
esac
