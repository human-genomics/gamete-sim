import numpy as np
import pytest

from gamete_sim.meiosis import sample_crossovers

N = 40_000


def counts(gamete, n=N):
    return np.bincount(gamete, minlength=n)


def assert_mean(c, expected):
    se = c.std(ddof=1) / np.sqrt(len(c))
    assert abs(c.mean() - expected) < 5 * se + 1e-12, (c.mean(), expected)


@pytest.mark.parametrize("length", [0.02, 0.3, 1.0, 3.3])
@pytest.mark.parametrize("nu, p", [(1.0, 0.0), (8.0, 0.0), (8.0, 0.06), (100.0, 0.05), (6.4, 0.19)])
def test_one_crossover_per_morgan(length, nu, p):
    """Stationary start: mean count equals length even on short chromosomes."""
    rng = np.random.default_rng(123)
    gamete, loc = sample_crossovers(rng, N, length, nu, p)
    assert_mean(counts(gamete), length)
    assert np.all((loc >= 0) & (loc < length))


def test_poisson_model_variance_equals_mean():
    rng = np.random.default_rng(1)
    c = counts(sample_crossovers(rng, N, 2.0, 1.0, 0.0)[0])
    assert abs(c.var() / c.mean() - 1.0) < 0.03


def test_interference_reduces_variance_and_short_gaps():
    rng = np.random.default_rng(2)
    gamete, loc = sample_crossovers(rng, N, 2.0, 9.0, 0.05)
    c = counts(gamete)
    assert c.var() / c.mean() < 0.7
    gaps = np.diff(loc)[np.diff(gamete) == 0]
    # Without interference about 18% of gaps would be under 10 cM.
    assert np.mean(gaps < 0.10) < 0.06


def test_uniform_positions_along_chromosome():
    """Stationarity: crossover density is flat along the chromosome."""
    rng = np.random.default_rng(3)
    _, loc = sample_crossovers(rng, N, 1.5, 9.0, 0.06)
    hist = np.histogram(loc, bins=15, range=(0, 1.5))[0]
    expected = len(loc) / 15
    assert np.all(np.abs(hist - expected) < 5 * np.sqrt(expected))


def test_sorted_output():
    rng = np.random.default_rng(4)
    gamete, loc = sample_crossovers(rng, 500, 3.0, 7.0, 0.07)
    order = np.lexsort((loc, gamete))
    assert np.array_equal(order, np.arange(len(loc)))


def test_deterministic():
    a = sample_crossovers(np.random.default_rng(5), 1000, 2.0, 7.0, 0.07)
    b = sample_crossovers(np.random.default_rng(5), 1000, 2.0, 7.0, 0.07)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))


def test_zero_length():
    gamete, loc = sample_crossovers(np.random.default_rng(0), 10, 0.0, 7.0, 0.07)
    assert len(gamete) == len(loc) == 0


@pytest.mark.parametrize("length, nu, p", [(-1, 7, 0.1), (1, 0, 0.1), (1, 7, 1.5)])
def test_invalid_parameters(length, nu, p):
    with pytest.raises(ValueError):
        sample_crossovers(np.random.default_rng(0), 10, length, nu, p)
