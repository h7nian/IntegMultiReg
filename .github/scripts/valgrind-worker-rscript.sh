#!/bin/sh
# CI-only PSOCK launcher: instrument the actual R child, not just its parent.
set -eu
: "${R_HOME:?R_HOME must identify the tested R installation}"
: "${IMR_VALGRIND_LOG_DIR:?A worker log directory is required}"
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
exec valgrind --tool=memcheck --leak-check=full --show-leak-kinds=definite \
  --errors-for-leak-kinds=definite --track-origins=yes --error-exitcode=1 \
  --log-file="${IMR_VALGRIND_LOG_DIR}/worker-%p.log" \
  "${R_HOME}/bin/exec/R" --no-echo --no-restore \
  -e "${worker_expression}" --args "$@"
