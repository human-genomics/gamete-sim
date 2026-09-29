"""Simulate gametes as hg38 haplotype segments.

Output: one row per segment, columns gamete, chrom, start, end, haplotype.
Coordinates are 1-based and inclusive (VCF convention). For every gamete and
chromosome the segments tile the whole chromosome:
  * the first segment starts at 1 and the last segment ends at the GRCh38 length
  * each next segment starts at the previous end + 1
  * haplotype alternates between 1 and 2 at each crossover
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AUTOSOMES, GRCH38_LENGTHS, SEX_SEED_KEY
from .meiosis import sample_crossovers

COLUMNS = ["gamete", "chrom", "start", "end", "haplotype"]


def read_interference(path: str | Path) -> pd.DataFrame:
    """Read a ped-sim interference file (#chr male_nu male_p female_nu female_p)."""
    table = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                        names=["chrom", "male_nu", "male_p", "female_nu", "female_p"])
    table["chrom"] = table["chrom"].astype(int)
    if list(table["chrom"]) != list(AUTOSOMES):
        raise ValueError("interference file must list chromosomes 1-22 in order")
    return table.set_index("chrom")


def stream(seed: int, sex: str, chrom: int) -> np.random.Generator:
    """Independent random stream for one (sex, chromosome)."""
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(SEX_SEED_KEY[sex], chrom)))


def breakpoints_bp(genetic_pos_cM: np.ndarray, map_cM: np.ndarray, map_pos: np.ndarray) -> np.ndarray:
    """Convert crossover genetic positions to the last base before the crossover.

    Linear interpolation between the flanking map markers, floored to a base,
    as in ped-sim. The segment before the crossover ends at the returned base;
    the next segment starts one base later.
    """
    i = np.searchsorted(map_cM, genetic_pos_cM, side="right") - 1
    i = np.clip(i, 0, len(map_cM) - 2)
    frac = (genetic_pos_cM - map_cM[i]) / (map_cM[i + 1] - map_cM[i])
    return np.floor(map_pos[i] + frac * (map_pos[i + 1] - map_pos[i])).astype(np.int64)


def cancel_pairs(gamete: np.ndarray, bp: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Drop crossovers that land between the same two bases an even number of times."""
    if len(bp) == 0:
        return gamete, bp
    key = np.stack([gamete, bp], axis=1)
    unique, counts = np.unique(key, axis=0, return_counts=True)
    odd = unique[counts % 2 == 1]
    return odd[:, 0], odd[:, 1]


def segments(n: int, gamete: np.ndarray, bp: np.ndarray, first_hap: np.ndarray, chrom_length: int):
    """Turn sorted per-gamete breakpoints into tiling segments."""
    n_seg = np.bincount(gamete, minlength=n) + 1
    seg_gamete = np.repeat(np.arange(n), n_seg)
    seg_rank = np.arange(n_seg.sum()) - np.repeat(np.cumsum(n_seg) - n_seg, n_seg)
    is_first = seg_rank == 0
    is_last = seg_rank == n_seg[seg_gamete] - 1

    start = np.empty(len(seg_gamete), dtype=np.int64)
    end = np.empty(len(seg_gamete), dtype=np.int64)
    start[is_first] = 1
    start[~is_first] = bp + 1
    end[~is_last] = bp
    end[is_last] = chrom_length
    hap0 = first_hap[seg_gamete]
    haplotype = np.where(seg_rank % 2 == 0, hap0, 3 - hap0)
    return seg_gamete, start, end, haplotype


def simulate_chromosome(chrom_map: pd.DataFrame, nu: float, p: float, sex: str,
                        chrom: int, n: int, seed: int) -> pd.DataFrame:
    rng = stream(seed, sex, chrom)
    first_hap = rng.integers(1, 3, size=n)

    map_cM = chrom_map[f"{sex}_cM"].to_numpy(float)
    map_pos = chrom_map["pos"].to_numpy(np.int64)
    length_morgans = (map_cM[-1] - map_cM[0]) / 100.0

    gamete, location = sample_crossovers(rng, n, length_morgans, nu, p)
    bp = breakpoints_bp(map_cM[0] + 100.0 * location, map_cM, map_pos)
    gamete, bp = cancel_pairs(gamete, bp)

    seg_gamete, start, end, haplotype = segments(n, gamete, bp, first_hap, GRCH38_LENGTHS[chrom])
    return pd.DataFrame({
        "gamete": seg_gamete + 1,
        "chrom": chrom,
        "start": start,
        "end": end,
        "haplotype": haplotype,
    })


def simulate(genetic_map: pd.DataFrame, interference: pd.DataFrame, sex: str,
             n: int, seed: int) -> pd.DataFrame:
    """Simulate `n` gametes from a parent of the given sex."""
    if sex not in SEX_SEED_KEY:
        raise ValueError(f"sex must be one of {list(SEX_SEED_KEY)}")
    parts = []
    for chrom, chrom_map in genetic_map.groupby("chrom", sort=True):
        params = interference.loc[chrom]
        parts.append(simulate_chromosome(chrom_map, params[f"{sex}_nu"], params[f"{sex}_p"],
                                         sex, int(chrom), n, seed))
    table = pd.concat(parts, ignore_index=True)
    table = table.sort_values(["gamete", "chrom", "start"], kind="stable", ignore_index=True)
    return table[COLUMNS]


def write(table: pd.DataFrame, path: str | Path) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, index=False, lineterminator="\n")
