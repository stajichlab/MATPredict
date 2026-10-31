"""GenBank as a `detect` input: sniff it, and turn it into the genome and proteome `detect` reads.

A GenBank file carries both the sequence and an annotation, so one file gives `detect` its genome, its
gene calls (the CDS proteins, which take the proteome route) and an annotation track for the report.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from MATPredict.detect.annotation_export import ConvertedGenome, convert_genbank

GENBANK_SUFFIXES = (".gb", ".gbk", ".gbff", ".genbank")


def is_genbank(path: str | Path) -> bool:
    """True when the file is GenBank: a known suffix, or a first non-blank line that starts with LOCUS."""
    path = Path(path)
    if path.suffix.lower() in GENBANK_SUFFIXES:
        return True
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if line.strip():
                    return line.startswith("LOCUS")
    except OSError:
        return False
    return False


@dataclass(frozen=True)
class GenbankInput:
    source: Path
    converted: ConvertedGenome


def prepare_genbank_input(path: str | Path, work_dir: Path) -> GenbankInput:
    """Write the genome FASTA and the coordinate-bearing proteome of a GenBank file into `work_dir`.

    A file with sequence but no CDS translations (an unannotated deposit) is not an error here: only the genome
    FASTA is written and `converted.proteome_fasta` is None, so `detect` searches the genome itself."""
    return GenbankInput(Path(path), convert_genbank(Path(path), Path(work_dir), require_proteins=False))
