#!/bin/bash
# Run from the MSI login node. Heavy work always goes through Slurm.
set -euo pipefail
source "$(dirname "$0")/environment.sh"
mode="${1:-}"
array=()
case "$mode" in
  setup|smoke) [[ $# == 1 ]] || exit 2 ;;
  chains) [[ $# == 1 ]] || exit 2; array=(--array="1-8%$MSI_CONCURRENCY") ;;
  simulation) [[ $# == 1 ]] || exit 2; array=(--array="1-150%$MSI_CONCURRENCY") ;;
  correlated) [[ $# == 1 ]] || exit 2; array=(--array="1-180%$MSI_CONCURRENCY") ;;
  collect) [[ $# == 2 && "$2" =~ ^(chains|simulation|correlated)$ ]] || exit 2 ;;
  *) echo 'Usage: bash msi/submit.sh setup|smoke|chains|simulation|correlated|collect TYPE' >&2; exit 2 ;;
esac
if [[ "$mode" != setup && ! -f "$MSI_LIBRARY/READY" ]]; then
  echo 'Setup has not finished successfully; inspect its Slurm log.' >&2; exit 2
fi
if [[ "$mode" != setup ]]; then
  cmp "$ROOT/msi/config.sh" "$MSI_LIBRARY/config.frozen.sh"
fi
if [[ "$mode" != setup && "$mode" != smoke && ! -f "$MSI_RESULTS/smoke/PASS" ]]; then
  echo 'Run the smoke job successfully before large calculations.' >&2; exit 2
fi
mkdir -p "$MSI_RESULTS/logs"
sbatch --parsable --account="$MSI_ACCOUNT" --partition="$MSI_PARTITION" \
  --nodes=1 --ntasks=1 --cpus-per-task=1 --mem="${MSI_MEMORY_OVERRIDE:-$MSI_MEMORY}" --time="${MSI_TIME_OVERRIDE:-$MSI_TIME}" \
  --job-name="imr-$mode" --chdir="$ROOT" \
  --output="$MSI_RESULTS/logs/%x-%A_%a.out" --error="$MSI_RESULTS/logs/%x-%A_%a.err" \
  ${array[@]+"${array[@]}"} "$ROOT/msi/job.sh" "$@"
