"""Reference panel: idiomorph-specific k-mer sets and the k-mers shared between idiomorphs."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from MATPredict.reads.kmers import canonical_kmers, revcomp


def read_fasta_seq(path: str | Path) -> str:
    """Concatenate all records of a FASTA file into one sequence."""
    return "".join(
        line.strip() for line in Path(path).read_text().splitlines() if line and not line.startswith(">")
    )


@dataclass
class Panel:
    k: int
    unique: dict[str, set[str]]
    shared: set[str]
    lengths: dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_sequences(cls, seqs: dict[str, str], k: int = 31) -> "Panel":
        """Unique k-mers of one idiomorph are absent from every other idiomorph; k-mers in two or more
        idiomorphs (shared flank ends) are the single-copy control and are never used for the call."""
        if len(seqs) < 2:
            raise ValueError("a panel needs at least two idiomorph sequences")
        sets = {name: canonical_kmers(seq, k) for name, seq in seqs.items()}
        shared: set[str] = set()
        for name, kms in sets.items():
            for other, okms in sets.items():
                if other != name:
                    shared |= kms & okms
        unique = {name: kms - shared for name, kms in sets.items()}
        for name, kms in unique.items():
            if not kms:
                raise ValueError(f"idiomorph {name} has no unique k-mers at k={k}")
        return cls(k=k, unique=unique, shared=shared, lengths={n: len(s) for n, s in seqs.items()})

    @classmethod
    def from_fastas(cls, fastas: dict[str, str | Path], k: int = 31) -> "Panel":
        return cls.from_sequences({n: read_fasta_seq(p) for n, p in fastas.items()}, k=k)

    def lookup(self) -> tuple[dict[str, tuple[str, int]], dict[str, list[str]]]:
        """Both-strand k-mer string -> (class, id), and class -> canonical k-mers (id order)."""
        classes = {**{n: sorted(s) for n, s in self.unique.items()}, "shared": sorted(self.shared)}
        table: dict[str, tuple[str, int]] = {}
        for cls_name, kms in classes.items():
            for i, kmer in enumerate(kms):
                table[kmer] = (cls_name, i)
                table[revcomp(kmer)] = (cls_name, i)
        return table, classes
