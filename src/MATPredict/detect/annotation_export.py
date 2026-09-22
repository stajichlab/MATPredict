"""Turn an annotation GenBank file into the (genome, proteome) pair `detect` needs.

Why this is its own tool rather than a one-liner: `search._parse_proteome_location`
REQUIRES every proteome defline to carry `contig:start-end:strand`, and raises
`ProteomeDeflineError` when it does not. That strictness is deliberate --
silently skipping unparseable deflines would turn a malformed input into a
confident "no MAT locus found", the exact class of false negative this pipeline
exists to avoid. The project handoff records that getting this format wrong
once cost an entire run.

No annotation pipeline writes those coordinates. ZygoLife LCG deflines are
`>EDD06_000001-T1 EDD06_000001`; NCBI and funannotate proteomes are no better.
So anything driving `detect` from an existing annotation has to rebuild them,
and every caller rebuilding them separately is how the format mistake recurs.

The proteome path is worth it: ~19 s/genome against ~4 min for genome-only on
the same code, because the diamond fast path skips the genome-wide tblastn
localization and most of the polish fan-out.

CAVEAT the caller must know: automatic annotation systematically UNDER-PREDICTS
MAT genes -- the pheromone precursors especially, which are short. A
proteome-driven run therefore finds the locus mainly through its flanking genes
(SLA2/APN2/COX13 style) and should not be treated as equivalent to a genome-only
run for the core genes. Measured on Exophiala: at the true locus the flanks hit
at 63-86% identity while the core MAT genes reach only 43-56%.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO


@dataclass(frozen=True)
class ConvertedGenome:
    """What `convert_genbank` wrote, and how much of it."""

    genome_fasta: Path
    proteome_fasta: Path
    contigs: int
    proteins: int


def convert_genbank(gbk_path: Path, out_dir: Path) -> ConvertedGenome:
    """Write `<stem>.fna` (contigs) and `<stem>.faa` (proteins) into `out_dir`.

    Contig names are preserved EXACTLY: the coordinates written into each
    protein defline are meaningless if the genome FASTA renames the contig they
    refer to.

    A CDS with no `/translation` is skipped rather than written as an empty
    record, which would make diamond fail on the whole file. The protein id is
    the `locus_tag`, else the `protein_id`, else a generated
    `<contig>_<n>` -- something must be there for the defline's first token.

    Raises ValueError when the file yields no proteins at all. An empty
    proteome is the failure mode that looks like success: diamond returns no
    hits and the genome is reported as having no locus.
    """
    gbk_path = Path(gbk_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    contigs = []
    proteins: list[str] = []
    for record in SeqIO.parse(str(gbk_path), "genbank"):
        contigs.append(record)
        for index, feature in enumerate(record.features):
            if feature.type != "CDS":
                continue
            translation = (feature.qualifiers.get("translation") or [None])[0]
            if not translation:
                continue
            identifier = (
                feature.qualifiers.get("locus_tag")
                or feature.qualifiers.get("protein_id")
                or [f"{record.id}_{index}"]
            )[0]
            start = int(feature.location.start) + 1  # GenBank is 0-based half-open
            end = int(feature.location.end)
            strand = "-" if feature.location.strand == -1 else "+"
            proteins.append(
                f">{identifier} {record.id}:{start}-{end}:{strand}\n{translation}\n"
            )

    if not proteins:
        raise ValueError(
            f"{gbk_path} yielded no CDS with a /translation; refusing to write an "
            f"empty proteome, which would be silently read as 'no locus found'"
        )

    stem = gbk_path.stem
    genome_fasta = out_dir / f"{stem}.fna"
    proteome_fasta = out_dir / f"{stem}.faa"
    SeqIO.write(contigs, str(genome_fasta), "fasta")
    proteome_fasta.write_text("".join(proteins))
    return ConvertedGenome(genome_fasta, proteome_fasta, len(contigs), len(proteins))
