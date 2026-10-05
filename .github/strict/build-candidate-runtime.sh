#!/usr/bin/env bash
set -euo pipefail
work="${RUNNER_TEMP}/strict-build"
evidence="${RUNNER_TEMP}/strict-evidence"
prefix="${RUNNER_TEMP}/strict-runtime"
mkdir -p "${work}" "${prefix}" "${evidence}"
exec > >(tee -a "${evidence}/runtime-build-console.log") 2>&1
cd "${work}"
curl --fail --location --retry 3 https://sourceware.org/pub/valgrind/valgrind-3.27.1.tar.bz2 -o valgrind.tar.bz2
echo '5d589152eb8071c02feab8ce6ab719e431a1fbc3e2b1700f5432632a8b9264dc  valgrind.tar.bz2' | sha256sum --check
if [[ -n "${IMR_R_SOURCE_ARCHIVE:-}" ]]; then
  cp "${IMR_R_SOURCE_ARCHIVE}" R-devel.tar.gz
else
  # Pinned deliberately: this audit must rebuild the same runtime. CRAN
  # removes older prerelease snapshots, so when this URL stops resolving the
  # pin and the checksum below have to be refreshed together, or
  # IMR_R_SOURCE_ARCHIVE set to a retained copy.
  curl --fail --location --retry 3 https://cran.r-project.org/src/base-prerelease/R-devel_2026-10-05_r90640.tar.gz -o R-devel.tar.gz
fi
# Pin refreshed to r90640 after CRAN withdrew the r90579 snapshot. The checksum
# below was taken from this archive as downloaded; this pin has not been
# re-audited on Anvil and it inherits acceptance of no earlier artifact.
echo 'd1a7da13ddb7f3c59a041c00d5bb842a806bce9b78a84931477152c24d82002d  R-devel.tar.gz' | sha256sum --check
sha256sum valgrind.tar.bz2 R-devel.tar.gz > "${evidence}/runtime-source-SHA256SUMS.txt"
if test -x "${prefix}/bin/valgrind" && test "$("${prefix}/bin/valgrind" --version)" = 'valgrind-3.27.1'; then
  cp "${prefix}/provenance/valgrind-configure.log" "${evidence}/"
else
  tar -xf valgrind.tar.bz2
  cd valgrind-3.27.1
  ./configure --prefix="${prefix}" --enable-only64bit > "${evidence}/valgrind-configure.log" 2>&1
  make -j2 > "${evidence}/valgrind-build.log" 2>&1
  make install >> "${evidence}/valgrind-build.log" 2>&1
fi
cd "${work}"
tar -xf R-devel.tar.gz
mkdir R-build
cd R-build
CPPFLAGS="-I${prefix}/include" CFLAGS='-g -O2 -Wall -pedantic -mtune=native' \
  CXXFLAGS='-g -O2 -Wall -pedantic -mtune=native' \
  FFLAGS='-g -O2 -mtune=native' FCFLAGS='-g -O2 -mtune=native' \
  ../R-devel/configure --prefix="${prefix}" --with-x=no \
    --enable-R-shlib --with-blas=no --with-lapack=no \
    --with-valgrind-instrumentation=2 > "${evidence}/R-configure.log" 2>&1
# Do not silently fall back to an uninstrumented R if headers were missed.
grep -E '^#define (HAVE_VALGRIND_MEMCHECK_H|VALGRIND_LEVEL|HAVE_PANGOCAIRO)' src/include/config.h | tee "${evidence}/instrumentation.txt"
grep -Eq '^#define VALGRIND_LEVEL 2$' src/include/config.h
# R Installation and Administration recommends Pango; do not silently use
# the cairo-ft fallback with its independent FreeType allocation leak.
grep -Eq '^#define HAVE_PANGOCAIRO 1$' src/include/config.h
make -j2 > "${evidence}/R-build.log" 2>&1
make install >> "${evidence}/R-build.log" 2>&1
mkdir -p "${prefix}/provenance"
cp "${evidence}/runtime-source-SHA256SUMS.txt" "${evidence}/instrumentation.txt" \
  "${evidence}/R-configure.log" "${evidence}/valgrind-configure.log" "${prefix}/provenance/"
