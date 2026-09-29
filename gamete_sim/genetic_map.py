"""Build the hg38 sex-specific genetic map from the Bherer et al. (2017) hg19 map.

We lift each map marker (not each simulated breakpoint) from hg19 to hg38 once.
The simulation then runs on hg38 coordinates directly, so no breakpoint can
fail to lift and no gamete is lost.

Marker filters, applied in this order per chromosome:
  1. unmapped      - the hg19 position is not in any chain block
  2. other_chrom   - the position lifts to another chromosome or an alt/random contig
  3. duplicate     - two hg19 markers lift to the same hg38 position (both dropped)
  4. out_of_order  - markers outside the longest run whose hg19 order agrees with
                     hg38 order (inversions and moved blocks between builds)

Filter 4 keeps the genetic map monotone in hg38 coordinates. Because cM is
monotone in hg19 order, keeping markers whose hg19 order is increasing in hg38
order keeps both the male and female cM columns monotone.
"""

from __future__ import annotations

import bisect
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AUTOSOMES, GRCH38_LENGTHS
from .liftover import Chain

BHERER_DIR = "Refined_EUR_genetic_map_b37"


def read_bherer(tar_path: str | Path) -> pd.DataFrame:
    """Read male and female autosomal maps from the Bherer tarball (no extraction)."""
    frames = []
    with tarfile.open(tar_path) as tar:
        for chrom in AUTOSOMES:
            sexes = {}
            for sex in ("male", "female"):
                member = tar.extractfile(f"{BHERER_DIR}/{sex}_chr{chrom}.txt")
                sexes[sex] = pd.read_csv(member, sep="\t", usecols=["pos", "cM"])
            male, female = sexes["male"], sexes["female"]
            if not np.array_equal(male["pos"].to_numpy(), female["pos"].to_numpy()):
                raise ValueError(f"chr{chrom}: male and female maps have different positions")
            frames.append(
                pd.DataFrame(
                    {
                        "chrom": chrom,
                        "pos_hg19": male["pos"].to_numpy(np.int64),
                        "male_cM": male["cM"].to_numpy(float),
                        "female_cM": female["cM"].to_numpy(float),
                    }
                )
            )
    df = pd.concat(frames, ignore_index=True)
    for chrom, sub in df.groupby("chrom"):
        pos = sub["pos_hg19"].to_numpy()
        if np.any(np.diff(pos) <= 0):
            raise ValueError(f"chr{chrom}: hg19 positions not strictly increasing")
        for col in ("male_cM", "female_cM"):
            if np.any(np.diff(sub[col].to_numpy()) < 0):
                raise ValueError(f"chr{chrom}: {col} decreases along hg19")
    return df


def longest_increasing_subsequence(values: np.ndarray) -> np.ndarray:
    """Boolean mask of one longest strictly increasing subsequence (O(n log n))."""
    n = len(values)
    tail_values: list = []
    tail_index: list[int] = []
    previous = np.full(n, -1, dtype=np.int64)
    for i, value in enumerate(values):
        k = bisect.bisect_left(tail_values, value)
        if k > 0:
            previous[i] = tail_index[k - 1]
        if k == len(tail_values):
            tail_values.append(value)
            tail_index.append(i)
        else:
            tail_values[k] = value
            tail_index[k] = i
    keep = np.zeros(n, dtype=bool)
    i = tail_index[-1] if tail_index else -1
    while i >= 0:
        keep[i] = True
        i = previous[i]
    return keep


def lift_map(hg19_map: pd.DataFrame, chain: Chain) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Lift the map to hg38. Returns (hg38 map, per-chromosome summary)."""
    kept_frames, summary = [], []
    for chrom, sub in hg19_map.groupby("chrom", sort=True):
        name = f"chr{chrom}"
        if chain.dst_sizes.get(name) != GRCH38_LENGTHS[chrom]:
            raise ValueError(f"{name}: chain hg38 size {chain.dst_sizes.get(name)} != GRCh38 {GRCH38_LENGTHS[chrom]}")

        dst_chrom, dst_pos, _ = chain.lift(name, sub["pos_hg19"].to_numpy())
        unmapped = dst_pos < 0
        other = ~unmapped & (dst_chrom != name)
        sub = sub.assign(pos=dst_pos)[~unmapped & ~other]

        duplicated = sub["pos"].duplicated(keep=False).to_numpy()
        sub = sub[~duplicated].sort_values("pos", kind="stable")

        in_order = longest_increasing_subsequence(sub["pos_hg19"].to_numpy())
        kept = sub[in_order]
        kept_frames.append(kept)

        orig = hg19_map[hg19_map["chrom"] == chrom]
        row = {
            "chrom": chrom,
            "markers_hg19": len(orig),
            "unmapped": int(unmapped.sum()),
            "other_chrom": int(other.sum()),
            "duplicate": int(duplicated.sum()),
            "out_of_order": int((~in_order).sum()),
            "markers_hg38": len(kept),
        }
        for sex in ("male", "female"):
            col = f"{sex}_cM"
            before = orig[col].iloc[-1] - orig[col].iloc[0]
            after = kept[col].iloc[-1] - kept[col].iloc[0]
            row[f"{sex}_cM_hg19"] = round(before, 4)
            row[f"{sex}_cM_hg38"] = round(after, 4)
        row["first_pos_hg38"] = int(kept["pos"].iloc[0])
        row["last_pos_hg38"] = int(kept["pos"].iloc[-1])
        row["chrom_length_hg38"] = GRCH38_LENGTHS[chrom]
        summary.append(row)

    hg38 = pd.concat(kept_frames, ignore_index=True)[["chrom", "pos", "male_cM", "female_cM", "pos_hg19"]]
    validate_map(hg38)
    return hg38, pd.DataFrame(summary)


def validate_map(genetic_map: pd.DataFrame) -> None:
    """Raise if the map is not usable for simulation."""
    if set(genetic_map["chrom"]) != set(AUTOSOMES):
        raise ValueError("map must contain exactly chromosomes 1-22")
    for chrom, sub in genetic_map.groupby("chrom"):
        pos = sub["pos"].to_numpy()
        if np.any(np.diff(pos) <= 0):
            raise ValueError(f"chr{chrom}: hg38 positions not strictly increasing")
        if pos[0] < 1 or pos[-1] > GRCH38_LENGTHS[chrom]:
            raise ValueError(f"chr{chrom}: positions outside the chromosome")
        for col in ("male_cM", "female_cM"):
            cm = sub[col].to_numpy()
            if np.any(np.diff(cm) < 0):
                raise ValueError(f"chr{chrom}: {col} decreases along hg38")
            if cm[-1] - cm[0] <= 0:
                raise ValueError(f"chr{chrom}: {col} has zero length")


def read_map(path: str | Path) -> pd.DataFrame:
    genetic_map = pd.read_csv(path, sep="\t")
    validate_map(genetic_map)
    return genetic_map


def build(bherer_tar: str | Path, chain_path: str | Path, out_dir: str | Path) -> pd.DataFrame:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    hg19_map = read_bherer(bherer_tar)
    hg38_map, summary = lift_map(hg19_map, Chain(chain_path))

    hg38_map.to_csv(out_dir / "genetic_map_hg38.tsv", sep="\t", index=False, float_format="%.10g")
    write_simmap(hg38_map, out_dir / "refined_mf_hg38.simmap")
    summary.to_csv(out_dir / "liftover_summary.tsv", sep="\t", index=False)
    return summary


def write_simmap(genetic_map: pd.DataFrame, path: str | Path) -> None:
    """Write the map in ped-sim's sex-specific format (usable with `ped-sim -m`)."""
    with open(path, "w") as handle:
        handle.write("#chr\tpos\tmale_cM\tfemale_cM\n")
        genetic_map[["chrom", "pos", "male_cM", "female_cM"]].to_csv(
            handle, sep="\t", index=False, header=False, float_format="%.10g"
        )
