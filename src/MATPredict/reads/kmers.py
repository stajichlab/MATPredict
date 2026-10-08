"""k-mer helpers and transparent FASTQ reading (plain, .gz, .zst)."""
from __future__ import annotations

import gzip
import io
import subprocess
from collections.abc import Iterator
from pathlib import Path

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def canonical_kmers(seq: str, k: int) -> set[str]:
    """Canonical (lexicographically smaller strand) k-mers of `seq`; k-mers with non-ACGT are skipped."""
    seq = seq.upper()
    out: set[str] = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if set(kmer) <= set("ACGT"):
            rc = revcomp(kmer)
            out.add(kmer if kmer <= rc else rc)
    return out


def _open_text(path: Path) -> io.TextIOBase:
    name = str(path)
    if name.endswith(".gz"):
        return gzip.open(path, "rt")
    if name.endswith(".zst"):
        proc = subprocess.Popen(["zstd", "-dc", name], stdout=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        return proc.stdout  # type: ignore[return-value]
    return open(path)


def iter_fastq_reads(path: str | Path) -> Iterator[str]:
    """Yield read sequences (upper-case) from a FASTQ file."""
    with _open_text(Path(path)) as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            seq = fh.readline().strip().upper()
            fh.readline()
            fh.readline()
            yield seq
