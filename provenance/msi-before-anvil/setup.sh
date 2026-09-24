#!/bin/bash
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 2
cd "$ROOT"
sha256sum --check --quiet SHA256SUMS
if [[ -e "$MSI_LIBRARY" ]]; then
  echo 'Library already exists. Never reinstall into a study library. Select a new MSI_LIBRARY path for a fresh setup.' >&2
  exit 2
fi
mkdir -p "$MSI_LIBRARY"
cp "$ROOT/msi/config.sh" "$MSI_LIBRARY/config.frozen.sh"
Rscript --vanilla -e 'stopifnot(getRversion() >= "4.4.0")'
gsl-config --version
# All non-base dependencies are bundled as sources in topological order.
while IFS=$'\t' read -r package version; do
  [[ "$package" == Package ]] && continue
  R CMD INSTALL --library="$MSI_LIBRARY" "$ROOT/dependencies/${package}_${version}.tar.gz"
done < "$ROOT/dependencies/dependencies.tsv"
R CMD INSTALL --library="$MSI_LIBRARY" "$ROOT/IntegMultiReg_0.2.0.tar.gz"
Rscript --vanilla "$ROOT/msi/assert-runtime.R"
Rscript --vanilla -e 'saveRDS(R.version, file.path(Sys.getenv("R_LIBS_USER"),"R-version.rds")); sessionInfo()' > "$MSI_LIBRARY/sessionInfo.txt"
gsl-config --version > "$MSI_LIBRARY/gsl-version.txt"
Rscript --vanilla "$ROOT/msi/freeze-runtime.R"
date -u +%FT%TZ > "$MSI_LIBRARY/READY"
