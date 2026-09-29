"""Checks and summaries for simulated gametes."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .constants import AUTOSOMES, GRCH38_LENGTHS
from .simulate import COLUMNS


def validate_gametes(table: pd.DataFrame, n: int) -> None:
    """Raise AssertionError if the table breaks any output invariant."""
    assert list(table.columns) == COLUMNS, f"columns {list(table.columns)} != {COLUMNS}"
    gametes = table["gamete"].to_numpy()
    assert gametes.min() == 1 and gametes.max() == n and len(np.unique(gametes)) == n, \
        f"expected gametes 1..{n}"
    assert set(table["haplotype"].unique()) <= {1, 2}, "haplotype must be 1 or 2"
    assert (table["start"] <= table["end"]).all(), "segment with start > end"

    order = np.lexsort((table["start"], table["chrom"], table["gamete"]))
    assert (order == np.arange(len(table))).all(), "rows not sorted by gamete, chrom, start"
    key = table[["gamete", "chrom"]].to_numpy()
    new_group = np.r_[True, (np.diff(key, axis=0) != 0).any(axis=1)]
    last_in_group = np.r_[new_group[1:], True]

    per_gamete = table.groupby("gamete")["chrom"].unique()
    assert all(sorted(c) == list(AUTOSOMES) for c in per_gamete), "each gamete needs chromosomes 1-22"

    start, end = table["start"].to_numpy(), table["end"].to_numpy()
    hap = table["haplotype"].to_numpy()
    lengths = table["chrom"].map(GRCH38_LENGTHS).to_numpy()
    assert (start[new_group] == 1).all(), "first segment must start at 1"
    assert (end[last_in_group] == lengths[last_in_group]).all(), "last segment must end at chromosome length"
    inner = ~new_group
    assert (start[inner] == end[np.r_[inner[1:], False]] + 1).all(), "segments not contiguous"
    assert (hap[inner] != hap[np.r_[inner[1:], False]]).all(), "haplotype must switch at each breakpoint"


def crossover_counts(table: pd.DataFrame) -> pd.DataFrame:
    """Crossovers per gamete (rows) and chromosome (columns)."""
    counts = table.groupby(["gamete", "chrom"]).size().unstack("chrom") - 1
    return counts[list(AUTOSOMES)]


def summarize(table: pd.DataFrame, genetic_map: pd.DataFrame, sex: str) -> pd.DataFrame:
    """Per-chromosome crossover counts against the map's genetic length."""
    counts = crossover_counts(table)
    col = f"{sex}_cM"
    length = genetic_map.groupby("chrom")[col].agg(lambda cm: (cm.iloc[-1] - cm.iloc[0]) / 100.0)
    summary = pd.DataFrame({
        "sex": sex,
        "chrom": list(AUTOSOMES),
        "map_length_M": length.loc[list(AUTOSOMES)].to_numpy(),
        "mean_crossovers": counts.mean().to_numpy(),
        "var_crossovers": counts.var(ddof=1).to_numpy(),
        "frac_haplotype1_at_start": table[table["start"] == 1].groupby("chrom")["haplotype"]
        .apply(lambda h: (h == 1).mean()).loc[list(AUTOSOMES)].to_numpy(),
    })
    total = counts.sum(axis=1)
    genome = pd.DataFrame([{
        "sex": sex,
        "chrom": "all",
        "map_length_M": summary["map_length_M"].sum(),
        "mean_crossovers": total.mean(),
        "var_crossovers": total.var(ddof=1),
        "frac_haplotype1_at_start": (table.loc[table["start"] == 1, "haplotype"] == 1).mean(),
    }])
    return pd.concat([summary, genome], ignore_index=True)
