#!/bin/bash -l
set -euo pipefail
export MSI_LOAD_MODULES=1
# sbatch copies this script to a spool path; --chdir from submit.sh locates the bundle.
source "$PWD/msi/environment.sh"
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit with sbatch; do not run on a login node.' >&2; exit 2; }
cd "$ROOT"
mode="$1"
mkdir -p "$MSI_RESULTS"
if [[ "$mode" == setup ]]; then
  bash "$ROOT/msi/setup.sh"
  exit
fi
cmp "$ROOT/msi/config.sh" "$MSI_LIBRARY/config.frozen.sh"
[[ -f "$MSI_LIBRARY/READY" ]]
sha256sum --check --quiet RUNTIME-SHA256SUMS
Rscript --vanilla "$ROOT/msi/assert-runtime.R"
module list 2>&1
Rscript --vanilla -e 'sessionInfo()'
gsl-config --version
id="${SLURM_ARRAY_TASK_ID:-1}"
case "$mode" in
  smoke) out="$MSI_RESULTS/smoke" ;;
  chains|simulation|correlated) out="$MSI_RESULTS/$mode/tasks/$(printf '%03d' "$id")" ;;
  collect) out="$MSI_RESULTS/$2/collection-control" ;;
  *) echo 'Unknown mode' >&2; exit 2 ;;
esac
mkdir -p "$out"
# Slurm/Linux lock is released even when the process is killed; file may remain.
exec 9>"$out/.lock"
flock -n 9 || { echo "Another process owns $out" >&2; exit 3; }
event="$out/slurm-${SLURM_JOB_ID}-${id}"
date -u +%FT%TZ > "$event.started"
trap 'rc=$?; printf "%s\n" "$rc" > "$event.exitcode"; date -u +%FT%TZ > "$event.ended"' EXIT
case "$mode" in
  smoke)
    Rscript --vanilla "$ROOT/paper/run-appendix-chains.R" --quick --chain 1 --out-dir "$out/chain"
    Rscript --vanilla "$ROOT/paper/original-experiments.R" --quick --experiment simulation \
      --reference code2017 --configuration 2 --replicate 1 --out-dir "$out/simulation"
    Rscript --vanilla "$ROOT/msi/check-smoke.R" "$out"
    date -u +%FT%TZ > "$out/PASS"
    ;;
  chains)
    Rscript --vanilla "$ROOT/paper/run-appendix-chains.R" --chain "$id" \
      --draws 1000000 --burnin 50000 --seed 100 --out-dir "$out"
    ;;
  simulation|correlated)
    replicates=50; [[ "$mode" == correlated ]] && replicates=30
    configuration=$(( (id-1)/replicates+1 )); replicate=$(( (id-1)%replicates+1 ))
    Rscript --vanilla "$ROOT/paper/original-experiments.R" --experiment "$mode" \
      --reference code2017 --marker-design fixed --draws 350000 --burnin 50000 \
      --k 10 --rounds 1 --replicates "$replicates" --workers 1 --seed 100 \
      --configuration "$configuration" --replicate "$replicate" --out-dir "$out"
    ;;
  collect)
    type="$2"
    Rscript --vanilla "$ROOT/msi/collect.R" "$type" "$MSI_RESULTS/$type/tasks" "$MSI_RESULTS/$type/combined"
    if [[ "$type" == chains ]]; then
      Rscript --vanilla "$ROOT/msi/diagnose-chains.R" "$MSI_RESULTS/$type/combined"
      Rscript --vanilla "$ROOT/paper/rank-chain-diagnostics.R" "$MSI_RESULTS/$type/combined" "$MSI_RESULTS/$type/combined/diagnostics"
      Rscript --vanilla "$ROOT/paper/chain-ranking-stability.R" "$MSI_RESULTS/$type/combined" "$MSI_RESULTS/$type/combined/diagnostics"
      Rscript --vanilla "$ROOT/msi/check-halves.R" "$MSI_RESULTS/$type/combined" "$MSI_RESULTS/$type/combined/half-drift"
    else
      Rscript --vanilla "$ROOT/paper/validate-original-simulation.R" "$MSI_RESULTS/$type/combined" --require-complete-study
      Rscript --vanilla "$ROOT/paper/validate-simulation-baselines.R" "$MSI_RESULTS/$type/combined"
      Rscript --vanilla "$ROOT/paper/summarize-original-simulation.R" "$MSI_RESULTS/$type/combined"
    fi
    date -u +%FT%TZ > "$MSI_RESULTS/$type/combined/VALIDATION-COMPLETED"
    ;;
esac
