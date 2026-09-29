"""Validation figures (static PNG, for the README and output/qc/)."""

from __future__ import annotations

from math import lgamma
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from .constants import AUTOSOMES, GRCH38_LENGTHS  # noqa: E402
from .meiosis import sample_crossovers  # noqa: E402

GAMETE_SIM = "#2a78d6"
PEDSIM = "#eb6834"
REFERENCE = "#898781"
TEXT = "#0b0b0b"
TEXT_2 = "#52514e"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SURFACE = "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "axes.edgecolor": AXIS,
    "axes.labelcolor": TEXT_2,
    "axes.titlecolor": TEXT,
    "axes.titlesize": 11,
    "axes.titleweight": "bold",
    "axes.titlelocation": "left",
    "axes.labelsize": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "xtick.color": TEXT_2,
    "ytick.color": TEXT_2,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.frameon": False,
    "legend.fontsize": 8,
    "font.family": "DejaVu Sans",
    "lines.linewidth": 1.6,
})


def _breakpoint_cM(breaks: pd.DataFrame, genetic_map: pd.DataFrame, sex: str) -> np.ndarray:
    """Genetic distances between adjacent crossovers on the same gamete and chromosome."""
    out = []
    for chrom, sub in breaks.groupby("chrom"):
        m = genetic_map[genetic_map["chrom"] == chrom]
        cm = np.interp(sub["bp"].to_numpy(), m["pos"].to_numpy(), m[f"{sex}_cM"].to_numpy())
        same = np.diff(sub["gamete"].to_numpy()) == 0
        out.append(np.diff(cm)[same])
    return np.concatenate(out)


def _poisson_spacing(genetic_map: pd.DataFrame, sex: str, n: int) -> np.ndarray:
    """Adjacent-crossover distances with no interference, same chromosome lengths."""
    rng = np.random.default_rng(0)
    out = []
    for chrom, m in genetic_map.groupby("chrom"):
        length = (m[f"{sex}_cM"].iloc[-1] - m[f"{sex}_cM"].iloc[0]) / 100.0
        gamete, loc = sample_crossovers(rng, n, length, nu=1.0, p=0.0)
        out.append((np.diff(loc) * 100.0)[np.diff(gamete) == 0])
    return np.concatenate(out)


def crossover_figure(summaries: dict[str, pd.DataFrame], counts: dict[str, pd.Series],
                     pedsim: dict[str, pd.DataFrame] | None, pedsim_totals: dict[str, np.ndarray] | None,
                     path: Path) -> None:
    """Per-chromosome mean crossovers and genome-wide crossover count distribution."""
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.6), gridspec_kw={"height_ratios": [1.2, 1]})
    for col, sex in enumerate(("male", "female")):
        s = summaries[sex]
        s = s[s["chrom"] != "all"]
        x = np.arange(len(AUTOSOMES))
        ax = axes[0, col]
        ax.bar(x, s["map_length_M"], width=0.7, color=GRID, label="Map length (expected mean)")
        ax.scatter(x - 0.15, s["mean_crossovers"], s=16, color=GAMETE_SIM, zorder=3, label="gamete-sim")
        if pedsim is not None:
            ax.scatter(x + 0.15, pedsim[sex]["pedsim_mean"], s=16, color=PEDSIM, zorder=3, label="ped-sim")
        ax.set_xticks(x, [str(c) for c in AUTOSOMES])
        ax.set_xlabel("Chromosome")
        ax.set_ylabel("Crossovers per gamete")
        ax.set_title(f"{sex.capitalize()} meiosis: mean crossovers per chromosome")
        ax.grid(axis="x", visible=False)
        if col == 0:
            ax.legend(loc="upper right")

        ax = axes[1, col]
        total = counts[sex]
        bins = np.arange(total.min() - 5, total.max() + 6)
        ax.hist(total, bins=bins, histtype="step", color=GAMETE_SIM, linewidth=1.6, density=True, label="gamete-sim")
        if pedsim_totals is not None:
            ax.hist(pedsim_totals[sex], bins=bins, histtype="step", color=PEDSIM,
                    linewidth=1.6, density=True, label="ped-sim")
        mean = total.mean()
        k = np.arange(bins[0], bins[-1])
        poisson = np.exp([ki * np.log(mean) - mean - lgamma(ki + 1) for ki in k])
        ax.plot(k + 0.5, poisson, color=REFERENCE, linestyle="--", linewidth=1.2, label="Poisson, same mean")
        ax.set_xlabel("Crossovers per gamete (genome-wide)")
        ax.set_ylabel("Fraction of gametes")
        ax.set_title(f"{sex.capitalize()}: gamete-sim mean {mean:.1f}, variance {total.var():.1f}")
        if col == 0:
            ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def interference_figure(breaks: dict[str, pd.DataFrame], pedsim_breaks: dict[str, pd.DataFrame] | None,
                        genetic_map: pd.DataFrame, n: int, path: Path) -> None:
    """Distance between adjacent crossovers: interference removes short distances."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.6))
    bins = np.arange(0, 205, 5)
    for ax, sex in zip(axes, ("male", "female")):
        ax.hist(_poisson_spacing(genetic_map, sex, n), bins=bins, histtype="step", density=True,
                color=REFERENCE, linestyle="--", linewidth=1.2, label="No interference")
        ax.hist(_breakpoint_cM(breaks[sex], genetic_map, sex), bins=bins, histtype="step",
                density=True, color=GAMETE_SIM, linewidth=1.6, label="gamete-sim")
        if pedsim_breaks is not None:
            ax.hist(_breakpoint_cM(pedsim_breaks[sex], genetic_map, sex), bins=bins, histtype="step",
                    density=True, color=PEDSIM, linewidth=1.6, label="ped-sim")
        ax.set_xlabel("Distance between adjacent crossovers (cM)")
        ax.set_ylabel("Density")
        ax.set_title(f"{sex.capitalize()}: crossover interference")
        ax.set_xlim(0, 200)
    axes[0].legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def landscape_figure(breaks: dict[str, pd.DataFrame], genetic_map: pd.DataFrame, n: int,
                     path: Path, chrom: int = 1, bin_bp: int = 500_000) -> None:
    """Crossover density along one chromosome against the hg38 map's expectation."""
    fig, axes = plt.subplots(2, 1, figsize=(11, 5.2), sharex=True)
    m = genetic_map[genetic_map["chrom"] == chrom]
    edges = np.arange(0, GRCH38_LENGTHS[chrom] + bin_bp, bin_bp)
    mid = (edges[:-1] + edges[1:]) / 2 / 1e6
    for ax, sex in zip(axes, ("male", "female")):
        cm = np.interp(edges, m["pos"].to_numpy(), m[f"{sex}_cM"].to_numpy())
        expected = n * np.diff(cm) / 100.0
        bp = breaks[sex].loc[breaks[sex]["chrom"] == chrom, "bp"].to_numpy()
        observed = np.histogram(bp, bins=edges)[0]
        ax.plot(mid, expected, color=AXIS, linewidth=3.0, label="Expected from hg38 map")
        ax.plot(mid, observed, color=GAMETE_SIM, linewidth=1.0, label="gamete-sim crossovers")
        r = np.corrcoef(expected, observed)[0, 1]
        ax.set_ylabel(f"Crossovers per {bin_bp // 1000} kb")
        ax.set_title(f"{sex.capitalize()}: chr{chrom} crossover landscape, {n:,} gametes (r = {r:.3f})")
    axes[0].legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), ncol=2)
    axes[1].set_xlabel(f"chr{chrom} position, GRCh38 (Mb)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
