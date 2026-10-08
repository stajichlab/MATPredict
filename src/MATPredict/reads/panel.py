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


def _kmers_in_long_runs(seq: str, k: int, unique: set[str], min_run: int) -> set[str]:
    """Canonical k-mers of `unique` that lie in a run of >= min_run consecutive positions whose k-mer is in `unique`."""
    seq = seq.upper()
    flags, kmers = [], []
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        if set(kmer) <= set("ACGT"):
            rc = revcomp(kmer)
            c = kmer if kmer <= rc else rc
            kmers.append(c); flags.append(c in unique)
        else:
            kmers.append(None); flags.append(False)
    keep: set[str] = set()
    i = 0
    while i < len(flags):
        if not flags[i]:
            i += 1
            continue
        j = i
        while j < len(flags) and flags[j]:
            j += 1
        if j - i >= min_run:
            keep.update(kmers[i:j])
        i = j
    return keep


@dataclass
class Panel:
    k: int
    unique: dict[str, set[str]]
    shared: set[str]
    lengths: dict[str, int] = field(default_factory=dict)

    @classmethod
    def from_sequences(cls, seqs: dict[str, str], k: int = 31, min_run: int = 100) -> "Panel":
        """Unique k-mers of one idiomorph are absent from every other idiomorph; k-mers in two or more
        idiomorphs (shared flank ends) are the single-copy control and are never used for the call.

        A SNP between two otherwise shared flanks makes up to k k-mers "unique" to each locus, and those k-mers
        match any strain with that flank allele. `min_run` keeps a unique k-mer only if it lies in a run of at
        least `min_run` consecutive unique positions of its locus (a SNP gives a run of at most k; an idiomorph
        core gives hundreds). Set `min_run=0` to keep every unique k-mer."""
        if len(seqs) < 2:
            raise ValueError("a panel needs at least two idiomorph sequences")
        sets = {name: canonical_kmers(seq, k) for name, seq in seqs.items()}
        shared: set[str] = set()
        for name, kms in sets.items():
            for other, okms in sets.items():
                if other != name:
                    shared |= kms & okms
        unique = {}
        for name, seq in seqs.items():
            only = sets[name] - shared
            if min_run > 0:
                only = _kmers_in_long_runs(seq, k, only, min_run)
            unique[name] = only
        for name, kms in unique.items():
            if not kms:
                raise ValueError(f"idiomorph {name} has no unique k-mers at k={k}")
        return cls(k=k, unique=unique, shared=shared, lengths={n: len(s) for n, s in seqs.items()})

    @classmethod
    def from_fastas(cls, fastas: dict[str, str | Path], k: int = 31, min_run: int = 100) -> "Panel":
        return cls.from_sequences({n: read_fasta_seq(p) for n, p in fastas.items()}, k=k, min_run=min_run)

    def lookup(self) -> tuple[dict[str, tuple[str, int]], dict[str, list[str]]]:
        """Both-strand k-mer string -> (class, id), and class -> canonical k-mers (id order)."""
        classes = {**{n: sorted(s) for n, s in self.unique.items()}, "shared": sorted(self.shared)}
        table: dict[str, tuple[str, int]] = {}
        for cls_name, kms in classes.items():
            for i, kmer in enumerate(kms):
                table[kmer] = (cls_name, i)
                table[revcomp(kmer)] = (cls_name, i)
        return table, classes
