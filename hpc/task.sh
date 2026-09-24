#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" && -n "${SLURM_ARRAY_TASK_ID:-}" ]] || exit 2
study="$(realpath "$1")"
source "$study/source/hpc/environment.sh"
id="$SLURM_ARRAY_TASK_ID"
out="$study/tasks/$(printf '%04d' "$id")"
mkdir -p "$out"
exec 9>"$out/.lock"
flock -n 9 || { echo 'Task is already running.' >&2; exit 3; }
event="$out/slurm-${SLURM_JOB_ID}"
date -u +%FT%TZ > "$event.started"
trap 'rc=$?; printf "%s\n" "$rc" > "$event.exitcode"; date -u +%FT%TZ > "$event.ended"' EXIT
/usr/bin/time -v -o "$event.resources" Rscript --vanilla "$study/source/hpc/run-task.R" "$study" "$id"
