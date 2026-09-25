#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" && ( $# == 2 || $# == 4 ) ]] || exit 2
source "$PWD/hpc/environment.sh"
candidate="$(realpath "$1")"
study="$2"
[[ -f "$IMR_DEPENDENCIES/READY" && -f "$candidate/UNIT-VALIDATED" ]]
if [[ $# == 4 ]]; then
  [[ "$3" == --pending-native-job && "$4" =~ ^[1-9][0-9]*$ ]] || exit 2
  if [[ ! -f "$candidate/NATIVE-VALIDATED" ]]; then
    [[ "$(squeue -h -j "$4" -o '%u')" == "$USER" ]] || exit 2
  fi
fi
Rscript --vanilla "$IMR_ROOT/hpc/check-dependencies.R"
Rscript --vanilla "$IMR_ROOT/hpc/prepare-study.R" "$candidate" "$study" "${@:3}"
study="$(realpath "$study")"
job=$(bash "$study/source/hpc/submit.sh" "$study" audit)
[[ "$job" =~ ^[0-9]+$ ]] || exit 2
controller=$(sbatch --parsable --account="$IMR_ACCOUNT" --partition="$IMR_PARTITION" \
  --dependency="afterany:$job" --nodes=1 --ntasks=1 --cpus-per-task=9 --mem=16G \
  --time=96:00:00 --job-name=imr-accept-audit --chdir="$IMR_ROOT" \
  --output="$study/logs/accept-audit-%j.log" "$study/source/hpc/advance.sh" "$study" audit)
printf '%s\taudit\t%s\n' "$controller" "$job" >> "$study/controllers.tsv"
