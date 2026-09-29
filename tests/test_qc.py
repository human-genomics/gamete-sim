import pytest

from gamete_sim import qc
from gamete_sim import simulate as sim


@pytest.fixture
def table(toy_map, toy_interference):
    return sim.simulate(toy_map, toy_interference, "female", 20, seed=2)


def first_inner_row(table):
    return table.index[table["start"] != 1][0]


def test_valid_table_passes(table):
    qc.validate_gametes(table, 20)


def test_wrong_count_fails(table):
    with pytest.raises(AssertionError, match="gametes"):
        qc.validate_gametes(table, 21)


def test_gap_fails(table):
    bad = table.copy()
    bad.loc[first_inner_row(bad), "start"] += 1
    with pytest.raises(AssertionError, match="contiguous"):
        qc.validate_gametes(bad, 20)


def test_same_haplotype_twice_fails(table):
    bad = table.copy()
    i = first_inner_row(bad)
    bad.loc[i, "haplotype"] = bad.loc[i - 1, "haplotype"]
    with pytest.raises(AssertionError, match="switch"):
        qc.validate_gametes(bad, 20)


def test_short_last_segment_fails(table):
    bad = table.copy()
    bad.loc[bad.index[-1], "end"] -= 1
    with pytest.raises(AssertionError, match="chromosome length"):
        qc.validate_gametes(bad, 20)


def test_start_zero_fails(table):
    bad = table.copy()
    bad.loc[0, "start"] = 0
    with pytest.raises(AssertionError, match="start at 1"):
        qc.validate_gametes(bad, 20)
