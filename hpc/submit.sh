#!/bin/bash -l
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: bash hpc/submit.sh STUDY audit|pilot|full' >&2; exit 2; }
source "$(dirname "$0")/environment.sh"
study="$(realpath "$1")"
phase="$2"
exec 8>"$study/.submission.lock"
flock -n 8 || { echo 'Another submitter is active.' >&2; exit 3; }
# There is one array at a time for this study, so %16 is a global cap.
if [[ -f "$study/last-array-id" ]]; then
  previous="$(cat "$study/last-array-id")"
  active="$(squeue -h -u "$USER" -o '%A' | awk -v id="$previous" '$1 == id')"
  [[ -z "$active" ]] || { echo 'Previous study array is still active.' >&2; exit 3; }
fi
ids="$(Rscript --vanilla "$study/source/hpc/select-tasks.R" "$study" "$phase")"
[[ -n "$ids" ]] || { echo 'No unfinished tasks in this phase.'; exit 0; }
memory_mib=16384
if [[ "$phase" == full ]]; then memory_mib="$(cat "$study/production-memory-mib")"; fi
[[ "$memory_mib" =~ ^[1-9][0-9]*$ ]] || exit 2
cpus=$(( (memory_mib + 1895) / 1896 ))
mkdir -p "$study/logs"
job=$(sbatch --parsable --account="$IMR_ACCOUNT" --partition="$IMR_PARTITION" \
  --nodes=1 --ntasks=1 --cpus-per-task="$cpus" --mem="${memory_mib}M" --time=48:00:00 \
  --array="${ids}%${IMR_CONCURRENCY}" --job-name="imr-$phase" --chdir="$IMR_ROOT" \
  --output="$study/logs/%x-%A_%a.log" "$study/source/hpc/task.sh" "$study")
job="${job%%;*}"
printf '%s\n' "$job" > "$study/last-array-id"
printf '%s\t%s\t%s\t%s\n' "$job" "$phase" "$memory_mib" "$ids" >> "$study/submissions.tsv"
printf '%s\n' "$job"
