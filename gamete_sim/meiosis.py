"""Crossover model: Housworth-Stahl two-pathway gamma model.

This is the model that ped-sim (Caballero et al. 2019) uses with `--intf`,
ported from `COInterfere::simStahl` in ped-sim, which was itself ported from
`simStahl` in Karl Broman's xoi R package. Parameters nu (interference
strength) and p (fraction of crossovers that escape interference) are the
sex- and chromosome-specific estimates of Campbell et al. (2015).

On a chromosome of genetic length L Morgans, for one gamete:
  * Chiasmata from the interfering pathway form a stationary renewal process
    with Gamma(shape=nu, rate=2*nu*(1-p)) spacing. Each chiasma involves the
    gamete's chromatid with probability 1/2.
  * Crossovers from the non-interfering pathway are Poisson(p*L), uniform on [0, L).
Both pathways together give 1 crossover per Morgan on average.

One change from ped-sim: ped-sim draws the first chiasma from a 10,000-bin
approximation of the stationary start distribution. Here it is drawn exactly
as U * Gamma(nu + 1, rate): the forward recurrence time of a stationary
renewal process is a uniform fraction of a length-biased gap, and the
length-biased Gamma(nu) gap is Gamma(nu + 1).
"""

from __future__ import annotations

import math

import numpy as np

# Treat nu == 1 as the no-interference (Poisson) model, as ped-sim does.
POISSON_NU_TOLERANCE = 1e-8


def sample_crossovers(
    rng: np.random.Generator, n: int, length: float, nu: float, p: float
) -> tuple[np.ndarray, np.ndarray]:
    """Sample crossover locations for `n` gametes on one chromosome.

    `length` is in Morgans. Returns (gamete_index, location) as flat arrays,
    sorted by gamete, then by location. Locations are in [0, length).
    """
    if not (0.0 <= p <= 1.0) or nu <= 0 or length < 0:
        raise ValueError(f"invalid parameters: length={length}, nu={nu}, p={p}")
    if length == 0:
        return np.empty(0, dtype=np.int64), np.empty(0)

    if abs(nu - 1.0) < POISSON_NU_TOLERANCE:
        counts = rng.poisson(length, size=n)
        gamete = np.repeat(np.arange(n), counts)
        location = rng.uniform(0.0, length, size=counts.sum())
        return _sort(gamete, location)

    gamete_parts, location_parts = [], []

    # Interfering pathway. Gamma spacing is in Morgans on the 4-strand bundle.
    if p < 1.0:
        scale = 1.0 / (2.0 * nu * (1.0 - p))
        current = rng.uniform(size=n) * rng.gamma(nu + 1.0, scale, size=n)
        chiasmata = [current]
        # Draw spacings in blocks until every gamete has passed the chromosome end.
        block = max(4, math.ceil(2.0 * (1.0 - p) * length * 1.5))
        while np.any(current < length):
            gaps = rng.gamma(nu, scale, size=(n, block))
            positions = current[:, None] + np.cumsum(gaps, axis=1)
            chiasmata.append(positions)
            current = positions[:, -1]
        chiasmata = np.column_stack(chiasmata)
        on_chromatid = rng.random(size=chiasmata.shape) < 0.5
        keep = (chiasmata < length) & on_chromatid
        rows, _ = np.nonzero(keep)
        gamete_parts.append(rows)
        location_parts.append(chiasmata[keep])

    # Non-interfering pathway.
    if p > 0.0:
        counts = rng.poisson(p * length, size=n)
        gamete_parts.append(np.repeat(np.arange(n), counts))
        location_parts.append(rng.uniform(0.0, length, size=counts.sum()))

    return _sort(np.concatenate(gamete_parts), np.concatenate(location_parts))


def _sort(gamete: np.ndarray, location: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.lexsort((location, gamete))
    return gamete[order].astype(np.int64), location[order]
