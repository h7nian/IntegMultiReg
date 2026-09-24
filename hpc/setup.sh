#!/bin/bash -l
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo 'Submit setup through Slurm.' >&2; exit 2; }
source "$PWD/hpc/environment.sh"
[[ ! -e "$IMR_DEPENDENCIES" ]] || { echo 'Choose a fresh dependency library.' >&2; exit 2; }
python3 - <<'PY'
import hashlib, json, os
from pathlib import Path
root = Path(os.environ['IMR_ROOT'])
for row in json.loads((root/'dependencies/sources.json').read_text()):
    path = root/'dependencies'/row['file']
    assert path.stat().st_size == row['bytes'], row['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row['sha256'], row['file']
print('All locked dependency archives verified.')
PY
mkdir -p "$IMR_DEPENDENCIES"
cp "$IMR_ROOT/hpc/config.sh" "$IMR_DEPENDENCIES/config.frozen.sh"
while IFS=$'\t' read -r package version; do
  [[ "$package" == Package ]] && continue
  R CMD INSTALL --library="$IMR_DEPENDENCIES" "$IMR_ROOT/dependencies/${package}_${version}.tar.gz"
done < "$IMR_ROOT/dependencies/dependencies.tsv"
bash "$IMR_ROOT/hpc/finalize-dependencies.sh"
