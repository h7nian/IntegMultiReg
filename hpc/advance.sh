#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" && $# == 2 ]] || exit 2
study="$(realpath "$1")"
phase="$2"
source "$study/source/hpc/environment.sh"
exec 7>"$study/.advance.lock"
flock -n 7 || { echo 'A phase controller is already active.' >&2; exit 3; }
while true; do
  if [[ "$phase" == full ]]; then
    checker="$study/source/hpc/collect.R"
    arguments=("$study")
  else
    checker="$study/source/hpc/accept-phase.R"
    arguments=("$study" "$phase")
  fi
  if Rscript --vanilla "$checker" "${arguments[@]}"; then
    case "$phase" in
      audit) next=pilot ;;
      pilot) next=full ;;
      full) exit 0 ;;
      *) exit 2 ;;
    esac
    remaining=$(Rscript --vanilla "$study/source/hpc/select-tasks.R" "$study" "$next")
    if [[ -z "$remaining" ]]; then phase="$next"; continue; fi
    job=$(bash "$study/source/hpc/submit.sh" "$study" "$next")
  else
    retry=$(Rscript --vanilla "$study/source/hpc/resource-retry.R" "$study" "$phase")
    mapfile -t request <<< "$retry"
    [[ ${#request[@]} == 3 ]] || exit 2
    job=$(IMR_TASK_IDS="${request[0]}" IMR_MEMORY_MIB="${request[1]}" IMR_WALL_HOURS="${request[2]}" \
      bash "$study/source/hpc/submit.sh" "$study" "$phase")
    next="$phase"
  fi
  break
done
[[ "$job" =~ ^[0-9]+$ ]] || { echo "$job"; exit 2; }
controller=$(sbatch --parsable --account="$IMR_ACCOUNT" --partition="$IMR_PARTITION" \
  --dependency="afterany:$job" --nodes=1 --ntasks=1 --cpus-per-task=9 --mem=16G \
  --time=96:00:00 --job-name="imr-accept-$next" --chdir="$IMR_ROOT" \
  --output="$study/logs/accept-$next-%j.log" "$study/source/hpc/advance.sh" "$study" "$next")
printf '%s\t%s\t%s\n' "$controller" "$next" "$job" >> "$study/controllers.tsv"
