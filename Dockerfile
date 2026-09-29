# =============================================================================
# gamete-sim — Docker image
# =============================================================================
# Build:
#   docker build -t gamete-sim .
#
# Run (outputs land in ./gamete-sim-data/output; downloads are cached in ./gamete-sim-data):
#   docker run --rm --user "$(id -u):$(id -g)" -v "$(pwd)/gamete-sim-data:/data" gamete-sim
#
# Settings: -e N_GAMETES=10000 -e SEED=42
# =============================================================================

# ---- Stage 1: build ped-sim (used only for the concordance check) ------------
FROM python:3.11-slim-bookworm AS pedsim
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl ca-certificates make g++ zlib1g-dev libboost-dev \
    && rm -rf /var/lib/apt/lists/*
COPY setup_pedsim.sh /build/
RUN bash /build/setup_pedsim.sh /build/tools

# ---- Stage 2: runtime --------------------------------------------------------
FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/human-genomics/gamete-sim"
LABEL org.opencontainers.image.description="hg38-native simulation of 10,000 human gametes with sex-specific crossover interference"
LABEL org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r /app/requirements.txt

COPY --from=pedsim /build/tools/bin/ped-sim /usr/local/bin/ped-sim
COPY main.sh download.sh setup_pedsim.sh checksums.sha256 /app/
COPY gamete_sim /app/gamete_sim

ENV PYTHON=python3 \
    GAMETE_SIM_DATA_DIR=/data \
    MPLCONFIGDIR=/tmp/matplotlib \
    PYTHONDONTWRITEBYTECODE=1

WORKDIR /data
ENTRYPOINT ["bash", "/app/main.sh"]
