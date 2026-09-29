"""Cross-check against ped-sim on the same hg38 map.

Runs the reference implementation (ped-sim --intf) for one parent with `n`
children, extracts the breakpoints of the haplotype that parent transmitted,
and compares them with gamete-sim's gametes. The two tools use different random
number generators, so we compare distributions, not values.
"""

from __future__ import annotations

import math
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AUTOSOMES

# In ped-sim BP output, a child's h0 comes from the father and h1 from the mother.
TRANSMITTED_HAP = {"male": "h0", "female": "h1"}
SEX_CODE = {"male": "M", "female": "F"}


def run(pedsim_bin: str, simmap: Path, interference: Path, sex: str, n: int,
        seed: int, work_dir: Path) -> Path:
    work_dir.mkdir(parents=True, exist_ok=True)
    def_file = work_dir / f"{sex}_parent.def"
    # n copies of a one-child pedigree (generation 2: 1 child, 1 branch). One
    # family with n children would make ped-sim compute IBD between all
    # n*(n-1)/2 sibling pairs.
    def_file.write_text(f"def {sex}_parent {n} 2 {SEX_CODE[sex]}\n2 1 1\n")
    prefix = work_dir / f"{sex}_parent"
    subprocess.run(
        [pedsim_bin, "-d", str(def_file), "-m", str(simmap), "-o", str(prefix),
         "--intf", str(interference), "--bp", "--seed", str(seed)],
        check=True, stdout=subprocess.DEVNULL,
    )
    return prefix.with_suffix(".bp")


def read_breakpoints(bp_file: Path, sex: str, n: int) -> pd.DataFrame:
    """Crossovers (gamete, chrom, bp) that the simulated parent transmitted.

    `bp` is the last base of the segment before the crossover, as in gamete-sim.
    """
    records = []
    gamete = 0
    with open(bp_file) as handle:
        for line in handle:
            fields = line.split()
            if "_g2-" not in fields[0] or fields[2] != TRANSMITTED_HAP[sex]:
                continue
            gamete += 1
            chrom = None
            ends: list[int] = []
            for token in fields[3:]:
                if "|" in token:
                    if chrom is not None:
                        records.extend((gamete, chrom, e) for e in ends[:-1])
                    chrom, ends = int(token.split("|")[0]), []
                else:
                    ends.append(int(token.split(":")[1]))
            records.extend((gamete, chrom, e) for e in ends[:-1])
    if gamete != n:
        raise ValueError(f"expected {n} transmitted haplotypes in {bp_file}, found {gamete}")
    return pd.DataFrame(records, columns=["gamete", "chrom", "bp"])


def breakpoints_from_segments(gametes: pd.DataFrame) -> pd.DataFrame:
    """Crossovers from a gamete-sim gamete table: the end of every non-final segment."""
    last = gametes.groupby(["gamete", "chrom"])["end"].transform("max")
    inner = gametes[gametes["end"] != last]
    return inner[["gamete", "chrom", "end"]].rename(columns={"end": "bp"}).reset_index(drop=True)


def ks_two_sample(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Two-sample Kolmogorov-Smirnov statistic and asymptotic p-value."""
    a, b = np.sort(a), np.sort(b)
    grid = np.concatenate([a, b])
    d = float(np.max(np.abs(np.searchsorted(a, grid, side="right") / len(a)
                            - np.searchsorted(b, grid, side="right") / len(b))))
    en = math.sqrt(len(a) * len(b) / (len(a) + len(b)))
    lam = (en + 0.12 + 0.11 / en) * d
    if lam < 1e-3:
        return d, 1.0
    p = 2.0 * sum((-1) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam) for k in range(1, 101))
    return d, float(min(max(p, 0.0), 1.0))


def to_cM(bp: np.ndarray, chrom_map: pd.DataFrame, sex: str) -> np.ndarray:
    return np.interp(bp, chrom_map["pos"].to_numpy(), chrom_map[f"{sex}_cM"].to_numpy())


def compare(ours: pd.DataFrame, theirs: pd.DataFrame, genetic_map: pd.DataFrame,
            sex: str, n: int) -> pd.DataFrame:
    """Per-chromosome comparison table."""
    rows = []
    for chrom in AUTOSOMES:
        chrom_map = genetic_map[genetic_map["chrom"] == chrom]
        length = (chrom_map[f"{sex}_cM"].iloc[-1] - chrom_map[f"{sex}_cM"].iloc[0]) / 100.0
        row = {"sex": sex, "chrom": chrom, "expected_mean": length}
        spacing = {}
        for label, table in (("gamete_sim", ours), ("pedsim", theirs)):
            sub = table[table["chrom"] == chrom]
            counts = np.bincount(sub["gamete"].to_numpy() - 1, minlength=n)
            row[f"{label}_mean"] = counts.mean()
            row[f"{label}_var"] = counts.var(ddof=1)
            cm = to_cM(sub["bp"].to_numpy(), chrom_map, sex)
            same_gamete = np.diff(sub["gamete"].to_numpy()) == 0
            spacing[label] = np.diff(cm)[same_gamete]
        a = ours.loc[ours["chrom"] == chrom, "bp"].to_numpy()
        b = theirs.loc[theirs["chrom"] == chrom, "bp"].to_numpy()
        row["position_ks_d"], row["position_ks_p"] = ks_two_sample(a, b)
        row["spacing_ks_d"], row["spacing_ks_p"] = ks_two_sample(spacing["gamete_sim"], spacing["pedsim"])
        rows.append(row)
    return pd.DataFrame(rows)
