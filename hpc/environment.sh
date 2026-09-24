#!/bin/bash
# Source from a batch job or a short diagnostic command.
IMR_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$IMR_ROOT/hpc/config.sh"
module purge
module load "${IMR_MODULES[@]}"
export IMR_ROOT IMR_DEPENDENCIES IMR_RUNS
export R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
unset R_MAKEVARS_USER R_MAKEVARS_SITE
export R_LIBS_USER="$IMR_DEPENDENCIES" R_LIBS="$IMR_DEPENDENCIES"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 BLIS_NUM_THREADS=1
export MAKEFLAGS=-j2
