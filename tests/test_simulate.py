import hashlib

import numpy as np
import pandas as pd
import pytest

from gamete_sim import qc
from gamete_sim import simulate as sim
from gamete_sim.constants import GRCH38_LENGTHS


def test_breakpoints_bp_interpolates_and_floors():
    map_cM = np.array([0.0, 1.0, 1.0, 3.0])
    map_pos = np.array([100, 200, 500, 700])
    bp = sim.breakpoints_bp(np.array([0.0, 0.5, 0.999, 1.0, 2.0]), map_cM, map_pos)
    # 1.0 cM sits on the flat stretch 200..500: the crossover goes after it.
    assert list(bp) == [100, 150, 199, 500, 600]


def test_cancel_pairs():
    gamete = np.array([0, 0, 0, 1, 1, 1, 2])
    bp = np.array([10, 10, 20, 5, 5, 5, 7])
    g, b = sim.cancel_pairs(gamete, bp)
    assert list(zip(g, b)) == [(0, 20), (1, 5), (2, 7)]


def test_segments_tile_chromosome():
    gamete = np.array([0, 0, 2])
    bp = np.array([100, 250, 40])
    first = np.array([1, 2, 2])
    g, start, end, hap = sim.segments(3, gamete, bp, first, 1000)
    rows = list(zip(g, start, end, hap))
    assert rows == [
        (0, 1, 100, 1), (0, 101, 250, 2), (0, 251, 1000, 1),
        (1, 1, 1000, 2),
        (2, 1, 40, 2), (2, 41, 1000, 1),
    ]


@pytest.mark.parametrize("sex", ["male", "female"])
def test_output_invariants(toy_map, toy_interference, sex):
    table = sim.simulate(toy_map, toy_interference, sex, 300, seed=7)
    qc.validate_gametes(table, 300)
    assert table["gamete"].nunique() == 300


def test_exact_gamete_count_even_without_crossovers(toy_map, toy_interference):
    table = sim.simulate(toy_map, toy_interference, "male", 1, seed=1)
    qc.validate_gametes(table, 1)


def test_mean_crossovers_match_map(toy_map, toy_interference):
    n = 4000
    table = sim.simulate(toy_map, toy_interference, "female", n, seed=3)
    xo = qc.crossover_counts(table).sum(axis=1)
    assert abs(xo.mean() - 22 * 2.2) < 5 * xo.std() / np.sqrt(n)


def test_no_crossovers_in_flat_stretch(toy_map, toy_interference):
    table = sim.simulate(toy_map, toy_interference, "male", 2000, seed=4)
    for chrom, m in toy_map.groupby("chrom"):
        lo, hi = m["pos"].iloc[10], m["pos"].iloc[15]
        ends = table[(table["chrom"] == chrom) & (table["end"] != GRCH38_LENGTHS[chrom])]["end"]
        assert not ((ends >= lo) & (ends < hi)).any()


def digest(table: pd.DataFrame) -> str:
    return hashlib.sha256(table.to_csv(index=False).encode()).hexdigest()


def test_deterministic_and_seed_dependent(toy_map, toy_interference):
    a = sim.simulate(toy_map, toy_interference, "male", 200, seed=11)
    b = sim.simulate(toy_map, toy_interference, "male", 200, seed=11)
    c = sim.simulate(toy_map, toy_interference, "male", 200, seed=12)
    assert digest(a) == digest(b)
    assert digest(a) != digest(c)


def test_first_haplotype_is_balanced(toy_map, toy_interference):
    table = sim.simulate(toy_map, toy_interference, "male", 3000, seed=5)
    first = table.loc[table["start"] == 1, "haplotype"]
    assert abs((first == 1).mean() - 0.5) < 0.02


def test_sexes_use_independent_streams(toy_map, toy_interference):
    male = sim.simulate(toy_map, toy_interference, "male", 50, seed=1)
    female = sim.simulate(toy_map, toy_interference, "female", 50, seed=1)
    first_male = male.loc[male["start"] == 1, "haplotype"].to_numpy()
    first_female = female.loc[female["start"] == 1, "haplotype"].to_numpy()
    assert not np.array_equal(first_male, first_female)


def test_read_interference(tmp_path):
    path = tmp_path / "intf.tsv"
    lines = ["#chr\tmale_nu\tmale_p\tfemale_nu\tfemale_p"]
    lines += [f"{c}\t8.0\t0.05\t7.0\t0.06" for c in range(1, 23)]
    path.write_text("\n".join(lines) + "\n")
    table = sim.read_interference(path)
    assert table.loc[22, "female_p"] == 0.06
