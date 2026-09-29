"""Check that gamete-sim's NumPy chain lifter equals UCSC liftOver on every map marker.

Usage (from the repo root, after download.sh):
    python3 validation/liftover_concordance.py [downloads_dir]

Downloads the UCSC liftOver binary (Linux x86_64 or macOS) into tools/bin/.
"""

from __future__ import annotations

import platform
import stat
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from gamete_sim.genetic_map import read_bherer  # noqa: E402
from gamete_sim.liftover import Chain  # noqa: E402

PLATFORMS = {
    ("Linux", "x86_64"): "linux.x86_64",
    ("Darwin", "x86_64"): "macOSX.x86_64",
    ("Darwin", "arm64"): "macOSX.arm64",
}


def liftover_binary() -> Path:
    path = REPO / "tools" / "bin" / "liftOver"
    if not path.exists():
        key = (platform.system(), platform.machine())
        if key not in PLATFORMS:
            sys.exit(f"UCSC liftOver has no binary for {key}")
        path.parent.mkdir(parents=True, exist_ok=True)
        url = f"https://hgdownload.soe.ucsc.edu/admin/exe/{PLATFORMS[key]}/liftOver"
        urllib.request.urlretrieve(url, path)
        path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def main() -> int:
    downloads = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO / "downloads"
    chain_path = downloads / "hg19ToHg38.over.chain.gz"
    markers = read_bherer(downloads / "Refined_EUR_genetic_map_b37.tar.gz")
    n = len(markers)

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        pd.DataFrame({
            "chrom": "chr" + markers["chrom"].astype(str),
            "start": markers["pos_hg19"] - 1,
            "end": markers["pos_hg19"],
            "id": np.arange(n),
        }).to_csv(tmp / "in.bed", sep="\t", header=False, index=False)
        subprocess.run([str(liftover_binary()), "-minMatch=1", str(tmp / "in.bed"), str(chain_path),
                        str(tmp / "out.bed"), str(tmp / "unmapped.bed")], check=True, capture_output=True)
        ucsc = pd.read_csv(tmp / "out.bed", sep="\t", header=None, names=["chrom", "start", "end", "id"])

    ucsc_chrom = np.full(n, "", dtype=object)
    ucsc_pos = np.full(n, -1, dtype=np.int64)
    ucsc_chrom[ucsc["id"]] = ucsc["chrom"]
    ucsc_pos[ucsc["id"]] = ucsc["end"]

    chain = Chain(chain_path)
    ours_chrom = np.empty(n, dtype=object)
    ours_pos = np.empty(n, dtype=np.int64)
    for chrom, idx in markers.groupby("chrom").indices.items():
        ours_chrom[idx], ours_pos[idx], _ = chain.lift(f"chr{chrom}", markers["pos_hg19"].to_numpy()[idx])

    same = (ours_chrom == ucsc_chrom) & (ours_pos == ucsc_pos)
    print(f"markers: {n:,}; lifted by UCSC: {(ucsc_pos > 0).sum():,}; lifted by gamete-sim: {(ours_pos > 0).sum():,}")
    print(f"identical chromosome and position: {same.sum():,} / {n:,}")
    return 0 if same.all() else 1


if __name__ == "__main__":
    sys.exit(main())
