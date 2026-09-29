"""Vectorised point liftover with a UCSC chain file.

A chain file lists ungapped alignment blocks between a source assembly
("target" in UCSC terms, here hg19) and a destination assembly ("query", here
hg38). A 1-bp position lifts only if it falls inside a block. In
hg19ToHg38.over.chain.gz no two blocks overlap on the source side, so every
position lifts to at most one place. The result equals `liftOver -minMatch=1`
on 1-bp BED features (validation/liftover_concordance.py checks this).
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Blocks:
    """All alignment blocks whose source is one chromosome, sorted by source start."""

    src_start: np.ndarray  # 0-based, inclusive
    src_end: np.ndarray  # 0-based, exclusive
    dst_start: np.ndarray  # 0-based, on the chain's destination strand
    dst_minus: np.ndarray  # True if the destination strand is "-"
    dst_size: np.ndarray  # destination chromosome length
    dst_name: np.ndarray  # destination chromosome name (object array)


class Chain:
    """Parsed chain file, indexed by source chromosome."""

    def __init__(self, path: str):
        rows: dict[str, list[tuple]] = {}
        self.dst_sizes: dict[str, int] = {}
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as handle:
            for line in handle:
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "chain":
                    src_name, src_strand = fields[2], fields[4]
                    dst_name, dst_size, dst_strand = fields[7], int(fields[8]), fields[9]
                    if src_strand != "+":
                        raise ValueError(f"source strand must be '+': {line.strip()}")
                    src_pos, dst_pos = int(fields[5]), int(fields[10])
                    dst_minus = dst_strand == "-"
                    self.dst_sizes[dst_name] = dst_size
                    block_rows = rows.setdefault(src_name, [])
                    continue
                size = int(fields[0])
                block_rows.append((src_pos, src_pos + size, dst_pos, dst_minus, dst_size, dst_name))
                if len(fields) == 3:
                    src_pos += size + int(fields[1])
                    dst_pos += size + int(fields[2])

        self.blocks: dict[str, Blocks] = {}
        for name, block_rows in rows.items():
            block_rows.sort(key=lambda r: r[0])
            cols = list(zip(*block_rows))
            blocks = Blocks(
                src_start=np.array(cols[0], dtype=np.int64),
                src_end=np.array(cols[1], dtype=np.int64),
                dst_start=np.array(cols[2], dtype=np.int64),
                dst_minus=np.array(cols[3], dtype=bool),
                dst_size=np.array(cols[4], dtype=np.int64),
                dst_name=np.array(cols[5], dtype=object),
            )
            if np.any(blocks.src_start[1:] < blocks.src_end[:-1]):
                raise ValueError(f"overlapping source blocks on {name}: point lift would be ambiguous")
            self.blocks[name] = blocks

    def lift(self, chrom: str, pos: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Lift 1-based positions on `chrom`.

        Returns (dst_chrom, dst_pos, dst_minus). Unmapped positions get
        dst_chrom "" and dst_pos -1.
        """
        pos = np.asarray(pos, dtype=np.int64)
        out_chrom = np.full(pos.shape, "", dtype=object)
        out_pos = np.full(pos.shape, -1, dtype=np.int64)
        out_minus = np.zeros(pos.shape, dtype=bool)
        blocks = self.blocks.get(chrom)
        if blocks is None:
            return out_chrom, out_pos, out_minus

        zero_based = pos - 1
        idx = np.searchsorted(blocks.src_start, zero_based, side="right") - 1
        safe = np.clip(idx, 0, None)
        hit = (idx >= 0) & (zero_based < blocks.src_end[safe])
        b = safe[hit]
        dst = blocks.dst_start[b] + (zero_based[hit] - blocks.src_start[b])
        minus = blocks.dst_minus[b]
        dst = np.where(minus, blocks.dst_size[b] - 1 - dst, dst)

        out_chrom[hit] = blocks.dst_name[b]
        out_pos[hit] = dst + 1
        out_minus[hit] = minus
        return out_chrom, out_pos, out_minus
