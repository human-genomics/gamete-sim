#!/usr/bin/env bash
# =============================================================================
# gamete-sim — simulate 10,000 human gametes on GRCh38 (hg38)
# =============================================================================
# Single entry point. Local runs need only python3 (3.10+) and curl.
#
# Usage:
#   bash main.sh
#
# Settings (environment variables):
#   N_GAMETES     gametes per sex                    (default 10000)
#   SEED          random seed                        (default 42)
#   GAMETE_SIM_DATA_DIR  where downloads/ output/ logs/ go  (default: this directory)
#   STRICT_CHECKSUMS=1  fail if the default outputs differ from checksums.sha256
#
# Steps:
#   1. Download inputs (Bherer 2017 map, UCSC chain, Campbell 2015 parameters)
#   2. Lift the sex-specific genetic map from hg19 to hg38 (once, marker by marker)
#   3. Simulate N_GAMETES male and N_GAMETES female gametes on hg38
#   4. Validate outputs, cross-check with ped-sim if available, draw figures
#   5. Compare outputs with the published checksums
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${GAMETE_SIM_DATA_DIR:-${SCRIPT_DIR}}"
N_GAMETES="${N_GAMETES:-10000}"
SEED="${SEED:-42}"
DOWNLOADS="${DATA_DIR}/downloads"
OUTPUT="${DATA_DIR}/output"
MAP_DIR="${OUTPUT}/genetic_map"
export PYTHONPATH="${SCRIPT_DIR}${PYTHONPATH:+:${PYTHONPATH}}"

# ---- Python (Docker sets PYTHON; local runs get a venv in tools/venv) --------
if [[ -z "${PYTHON:-}" ]]; then
    VENV_DIR="${SCRIPT_DIR}/tools/venv"
    PYTHON="${VENV_DIR}/bin/python3"
    if [[ -x "${PYTHON}" ]]; then
        echo "[skip] Python venv"
    else
        echo "[setup] Creating Python venv in tools/venv/ ..."
        python3 -m venv "${VENV_DIR}"
        "${PYTHON}" -m pip install --quiet --upgrade pip
        "${PYTHON}" -m pip install --quiet --no-cache-dir -r "${SCRIPT_DIR}/requirements.txt"
        echo "[done] Python venv ready"
    fi
fi

# ---- Logging: everything below goes to the terminal and to logs/ -------------
mkdir -p "${DATA_DIR}/logs" "${OUTPUT}"
LOG_FILE="${DATA_DIR}/logs/run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1
echo "gamete-sim: N_GAMETES=${N_GAMETES} SEED=${SEED} DATA_DIR=${DATA_DIR}"

echo "==> [1/5] Download inputs"
bash "${SCRIPT_DIR}/download.sh" "${DOWNLOADS}"

echo "==> [2/5] Lift the sex-specific genetic map to hg38"
"${PYTHON}" -m gamete_sim build-map \
    --bherer "${DOWNLOADS}/Refined_EUR_genetic_map_b37.tar.gz" \
    --chain "${DOWNLOADS}/hg19ToHg38.over.chain.gz" \
    --out-dir "${MAP_DIR}"

echo "==> [3/5] Simulate gametes"
for sex in male female; do
    "${PYTHON}" -m gamete_sim simulate \
        --map "${MAP_DIR}/genetic_map_hg38.tsv" \
        --interference "${DOWNLOADS}/nu_p_campbell.tsv" \
        --sex "${sex}" --n "${N_GAMETES}" --seed "${SEED}" \
        --out "${OUTPUT}/gametes_${sex}_hg38.csv"
done

echo "==> [4/5] Validate, cross-check with ped-sim, draw figures"
PEDSIM="$(command -v ped-sim || true)"
if [[ -z "${PEDSIM}" ]]; then
    PEDSIM="${SCRIPT_DIR}/tools/bin/ped-sim"
    if [[ ! -x "${PEDSIM}" ]]; then
        if command -v make >/dev/null && command -v g++ >/dev/null \
            && bash "${SCRIPT_DIR}/setup_pedsim.sh" "${SCRIPT_DIR}/tools"; then
            :
        else
            echo "  [skip] could not build ped-sim (optional; see setup_pedsim.sh)"
            PEDSIM=""
        fi
    fi
fi
"${PYTHON}" -m gamete_sim report \
    --map "${MAP_DIR}/genetic_map_hg38.tsv" \
    --interference "${DOWNLOADS}/nu_p_campbell.tsv" \
    --gametes-dir "${OUTPUT}" --n "${N_GAMETES}" --seed "${SEED}" \
    ${PEDSIM:+--pedsim-bin "${PEDSIM}"} \
    --out-dir "${OUTPUT}/qc"

echo "==> [5/5] Compare with published checksums"
if [[ "${N_GAMETES}" == "10000" && "${SEED}" == "42" ]]; then
    if command -v sha256sum >/dev/null 2>&1; then HASH=(sha256sum); else HASH=(shasum -a 256); fi
    MISMATCH=0
    while read -r expected file; do
        actual="$("${HASH[@]}" "${OUTPUT}/${file}" | awk '{print $1}')"
        if [[ "${actual}" == "${expected}" ]]; then
            echo "  [identical] ${file}"
        else
            echo "  [DIFFERS] ${file}: ${actual} (published ${expected})"
            MISMATCH=1
        fi
    done < "${SCRIPT_DIR}/checksums.sha256"
    if [[ "${MISMATCH}" == "1" && "${STRICT_CHECKSUMS:-0}" == "1" ]]; then
        exit 1
    fi
else
    echo "  [skip] non-default N_GAMETES or SEED"
fi

echo "Done. Main outputs:"
echo "  ${OUTPUT}/gametes_male_hg38.csv"
echo "  ${OUTPUT}/gametes_female_hg38.csv"
echo "Log: ${LOG_FILE}"
