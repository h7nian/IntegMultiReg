#!/usr/bin/env bash
set -euo pipefail
work="${RUNNER_TEMP}/strict-build"
evidence="${RUNNER_TEMP}/strict-evidence"
prefix="${RUNNER_TEMP}/strict-runtime"
mkdir -p "${work}" "${prefix}"
cd "${work}"
curl --fail --location --retry 3 https://sourceware.org/pub/valgrind/valgrind-3.27.1.tar.bz2 -o valgrind.tar.bz2
echo '5d589152eb8071c02feab8ce6ab719e431a1fbc3e2b1700f5432632a8b9264dc  valgrind.tar.bz2' | sha256sum --check
curl --fail --location --retry 3 https://cran.r-project.org/src/base-prerelease/R-devel_2026-09-06_r90498.tar.gz -o R-devel.tar.gz
sha256sum valgrind.tar.bz2 R-devel.tar.gz > "${evidence}/runtime-source-SHA256SUMS.txt"
tar -xf valgrind.tar.bz2
cd valgrind-3.27.1
./configure --prefix="${prefix}" --enable-only64bit > "${evidence}/valgrind-configure.log" 2>&1
make -j2 > "${evidence}/valgrind-build.log" 2>&1
make install >> "${evidence}/valgrind-build.log" 2>&1
cd "${work}"
tar -xf R-devel.tar.gz
mkdir R-build
cd R-build
CPPFLAGS="-I${prefix}/include" CFLAGS='-g -O2' \
  ../R-devel/configure --prefix="${prefix}" --with-x=no \
    --enable-R-shlib --with-blas=no --with-lapack=no \
    --with-valgrind-instrumentation=2 > "${evidence}/R-configure.log" 2>&1
# Do not silently fall back to an uninstrumented R if headers were missed.
grep -E '^#define (HAVE_VALGRIND_MEMCHECK_H|VALGRIND_LEVEL)' src/include/config.h | tee "${evidence}/instrumentation.txt"
grep -Eq '^#define VALGRIND_LEVEL 2$' src/include/config.h
make -j2 > "${evidence}/R-build.log" 2>&1
make install >> "${evidence}/R-build.log" 2>&1
echo "${prefix}/bin" >> "${GITHUB_PATH}"
