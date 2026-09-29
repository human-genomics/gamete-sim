#!/usr/bin/env bash
# =============================================================================
# download.sh — fetch the three pipeline inputs and verify their SHA-256
# =============================================================================
# Usage:
#   bash download.sh <downloads_dir>
#
# Files (~39 MB in total):
#   Refined_EUR_genetic_map_b37.tar.gz  Bherer et al. 2017 sex-specific EUR maps (hg19)
#   hg19ToHg38.over.chain.gz            UCSC hg19 -> hg38 liftOver chain
#   nu_p_campbell.tsv                   Campbell et al. 2015 interference parameters (from ped-sim)
#
# A file that exists and passes its checksum is not downloaded again.
# =============================================================================
set -euo pipefail

DEST="${1:?usage: bash download.sh <downloads_dir>}"
mkdir -p "${DEST}"

BHERER_COMMIT="0456f74ab8a2dce4d5172914bb799ff3d26d0b43"
PEDSIM_COMMIT="44d0c7abf6206479d5d6df926ecca086a8bbee32"  # ped-sim v1.4.7

if command -v sha256sum >/dev/null 2>&1; then
    sha256() { sha256sum "$1" | awk '{print $1}'; }
else
    sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
fi

download() {
    local name="$1" url="$2" expected="$3"
    local path="${DEST}/${name}"
    if [[ -f "${path}" && "$(sha256 "${path}")" == "${expected}" ]]; then
        echo "  [skip] ${name} (checksum ok)"
        return 0
    fi
    echo "  [download] ${name}"
    curl -sSfL --retry 3 -o "${path}.part" "${url}"
    local actual
    actual="$(sha256 "${path}.part")"
    if [[ "${actual}" != "${expected}" ]]; then
        echo "  [CHECKSUM FAILED] ${name}: expected ${expected}, got ${actual}" >&2
        exit 1
    fi
    mv "${path}.part" "${path}"
    echo "  [verified] ${name}"
}

download "Refined_EUR_genetic_map_b37.tar.gz" \
    "https://raw.githubusercontent.com/cbherer/Bherer_etal_SexualDimorphismRecombination/${BHERER_COMMIT}/Refined_EUR_genetic_map_b37.tar.gz" \
    "eae6191c837a8dc3f9e418b590d69e92a44cb0a3281d7fc5341b4987265b29ec"

download "hg19ToHg38.over.chain.gz" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz" \
    "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"

download "nu_p_campbell.tsv" \
    "https://raw.githubusercontent.com/williamslab/ped-sim/${PEDSIM_COMMIT}/interfere/nu_p_campbell.tsv" \
    "825648443c400c89103295dc0a29804c969164c54e2570a3084403bb778ca1a1"
