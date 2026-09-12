#!/bin/sh
# macOS's protected shell clears DYLD_* on entry. Restore it here, then
# directly exec R (not its shell wrapper) for base parallel's PSOCK command.
set -eu
: "${R_HOME:?R_HOME must identify the tested R installation}"
: "${SANITIZER_RUNTIME:?SANITIZER_RUNTIME must identify the sanitizer dylibs}"
export DYLD_INSERT_LIBRARIES="${SANITIZER_RUNTIME}"

case "${1-}" in
  --default-packages=*)
    R_DEFAULT_PACKAGES=${1#--default-packages=}
    export R_DEFAULT_PACKAGES
    shift
    ;;
esac
if [ "${1-}" != "-e" ] || [ "$#" -lt 2 ]; then
  echo "Expected base parallel's Rscript -e worker command" >&2
  exit 2
fi
worker_expression=$2
shift 2
exec "${R_HOME}/bin/exec/R" --no-echo --no-restore \
  -e "${worker_expression}" --args "$@"
