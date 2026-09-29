from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gamete_sim.constants import AUTOSOMES, GRCH38_LENGTHS

REPO = Path(__file__).resolve().parents[1]
DOWNLOADS = REPO / "downloads"


@pytest.fixture
def toy_map() -> pd.DataFrame:
    """Small map on all 22 GRCh38 autosomes: 50 markers, uneven rates, one flat stretch."""
    rng = np.random.default_rng(1)
    frames = []
    for chrom in AUTOSOMES:
        length = GRCH38_LENGTHS[chrom]
        pos = np.linspace(1_000_000, length - 1_000_000, 50).astype(np.int64)
        steps = rng.gamma(0.5, 1.0, size=(2, 49))
        steps[:, 10:15] = 0.0  # no recombination here
        male = np.r_[0, np.cumsum(steps[0])] / steps[0].sum() * 150.0
        female = np.r_[0, np.cumsum(steps[1])] / steps[1].sum() * 220.0
        frames.append(pd.DataFrame({"chrom": chrom, "pos": pos, "male_cM": male,
                                    "female_cM": female, "pos_hg19": pos}))
    return pd.concat(frames, ignore_index=True)


@pytest.fixture
def toy_interference() -> pd.DataFrame:
    return pd.DataFrame({"chrom": list(AUTOSOMES), "male_nu": 9.0, "male_p": 0.06,
                         "female_nu": 7.0, "female_p": 0.07}).set_index("chrom")


def need_download(name: str) -> Path:
    path = DOWNLOADS / name
    if not path.exists():
        pytest.skip(f"{name} not downloaded (run download.sh)")
    return path
