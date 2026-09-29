#!/usr/bin/env bash
# =============================================================================
# setup_pedsim.sh — build ped-sim v1.4.7 (optional)
# =============================================================================
# ped-sim is NOT needed to make the gametes. The pipeline uses it only to
# cross-check gamete-sim against the reference implementation on the same hg38 map.
# The Docker image always includes it. For local runs this needs make, a C++
# compiler, zlib headers and Boost headers (Debian/Ubuntu: build-essential
# zlib1g-dev libboost-dev).
#
# Usage:
#   bash setup_pedsim.sh [prefix]     # installs <prefix>/bin/ped-sim (default: ./tools)
# =============================================================================
set -euo pipefail

PEDSIM_COMMIT="44d0c7abf6206479d5d6df926ecca086a8bbee32"  # v1.4.7, 25 Jun 2026
PREFIX="${1:-$(cd "$(dirname "$0")" && pwd)/tools}"
BIN="${PREFIX}/bin/ped-sim"

if [[ -x "${BIN}" ]]; then
    echo "  [skip] ped-sim already built at ${BIN}"
    exit 0
fi

mkdir -p "${PREFIX}/src" "${PREFIX}/bin"
echo "  [download] ped-sim ${PEDSIM_COMMIT:0:7}"
curl -sSfL --retry 3 "https://github.com/williamslab/ped-sim/archive/${PEDSIM_COMMIT}.tar.gz" \
    | tar xz -C "${PREFIX}/src"
echo "  [build] ped-sim"
JOBS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 2)"
make -C "${PREFIX}/src/ped-sim-${PEDSIM_COMMIT}" distribute -j"${JOBS}" >/dev/null 2>&1
cp "${PREFIX}/src/ped-sim-${PEDSIM_COMMIT}/ped-sim" "${BIN}"
echo "  [installed] ${BIN}"
