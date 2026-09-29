import itertools

import numpy as np
import pytest

from gamete_sim import genetic_map
from gamete_sim.constants import AUTOSOMES
from gamete_sim.liftover import Chain

from .conftest import need_download


def brute_force_lis_length(values):
    for size in range(len(values), 0, -1):
        for combo in itertools.combinations(values, size):
            if all(a < b for a, b in zip(combo, combo[1:])):
                return size
    return 0


@pytest.mark.parametrize("values, expected", [
    ([], []),
    ([5], [True]),
    ([1, 2, 3], [True, True, True]),
    ([3, 1, 2], [False, True, True]),
    ([1, 5, 2, 3, 4], [True, False, True, True, True]),
    ([2, 2, 2], [False, False, True]),
])
def test_lis_examples(values, expected):
    assert list(genetic_map.longest_increasing_subsequence(np.array(values))) == expected


def test_lis_is_longest_and_increasing():
    rng = np.random.default_rng(0)
    for _ in range(200):
        values = rng.integers(0, 12, size=rng.integers(1, 10))
        keep = genetic_map.longest_increasing_subsequence(values)
        kept = values[keep]
        assert np.all(np.diff(kept) > 0)
        assert keep.sum() == brute_force_lis_length(list(values))


def test_validate_map_rejects_decreasing_cM(toy_map):
    bad = toy_map.copy()
    bad.loc[5, "male_cM"] = -1.0
    with pytest.raises(ValueError, match="decreases"):
        genetic_map.validate_map(bad)


def test_validate_map_accepts_toy(toy_map):
    genetic_map.validate_map(toy_map)


def test_real_map_lift():
    tar = need_download("Refined_EUR_genetic_map_b37.tar.gz")
    chain = need_download("hg19ToHg38.over.chain.gz")
    hg19 = genetic_map.read_bherer(tar)
    hg38, summary = genetic_map.lift_map(hg19, Chain(chain))
    assert list(summary["chrom"]) == list(AUTOSOMES)
    assert summary["markers_hg38"].sum() / summary["markers_hg19"].sum() > 0.999
    for sex in ("male", "female"):
        np.testing.assert_allclose(summary[f"{sex}_cM_hg38"], summary[f"{sex}_cM_hg19"], atol=1e-6)
    genetic_map.validate_map(hg38)
