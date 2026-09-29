import numpy as np
import pytest

from gamete_sim.liftover import Chain

from .conftest import need_download

# chrA:[10,30) -> chrA:[50,70), gap, chrA:[35,60) -> chrA:[80,105) on "+";
# chrA:[70,80) -> chrB reverse-strand [5,15) of a 50-bp chromosome.
TOY_CHAIN = """\
chain 1000 chrA 100 + 10 60 chrA 200 + 50 105 1
20 5 10
25

chain 500 chrA 100 + 70 80 chrB 50 - 5 15 2
10
"""


@pytest.fixture
def toy_chain(tmp_path):
    path = tmp_path / "toy.chain"
    path.write_text(TOY_CHAIN)
    return Chain(path)


def test_plus_strand_blocks(toy_chain):
    chrom, pos, minus = toy_chain.lift("chrA", np.array([11, 30, 36, 60]))
    assert list(pos) == [51, 70, 81, 105]
    assert list(chrom) == ["chrA"] * 4
    assert not minus.any()


def test_gaps_and_outside_are_unmapped(toy_chain):
    chrom, pos, _ = toy_chain.lift("chrA", np.array([1, 10, 31, 35, 61, 100]))
    assert list(pos) == [-1, -1, -1, -1, -1, -1]
    assert list(chrom) == [""] * 6


def test_minus_strand_block(toy_chain):
    chrom, pos, minus = toy_chain.lift("chrA", np.array([71, 80]))
    # 0-based reverse offsets 5 and 14 on a 50-bp chromosome -> 1-based 45 and 36
    assert list(pos) == [45, 36]
    assert list(chrom) == ["chrB", "chrB"]
    assert minus.all()


def test_unknown_chromosome(toy_chain):
    _, pos, _ = toy_chain.lift("chrZ", np.array([5]))
    assert list(pos) == [-1]


def test_destination_sizes(toy_chain):
    assert toy_chain.dst_sizes == {"chrA": 200, "chrB": 50}


def test_overlapping_blocks_rejected(tmp_path):
    path = tmp_path / "bad.chain"
    path.write_text("chain 1 chrA 100 + 0 20 chrA 100 + 0 20 1\n20\n\n"
                    "chain 1 chrA 100 + 10 30 chrA 100 + 50 70 2\n20\n")
    with pytest.raises(ValueError, match="overlapping"):
        Chain(path)


def test_real_chain_known_positions():
    chain = Chain(need_download("hg19ToHg38.over.chain.gz"))
    # First marker of the Bherer chr1 map (hg19 777122 -> hg38 841742).
    chrom, pos, _ = chain.lift("chr1", np.array([777122]))
    assert (chrom[0], pos[0]) == ("chr1", 841742)
