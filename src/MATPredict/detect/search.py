"""diamond (fast path) and exonerate (fallback path) search wrappers.

Reference-protein headers are `{record_id}|gene{gene_index}|{gene_name}`
(the form `reference_fasta.build_reference_fasta` rewrites
`gff_export.write_proteins_fasta`'s real headers into) -- parsed back out
here to recover which curated record and gene a hit corresponds to.

--outfmt design (fast path): diamond blastp is a protein-vs-protein search,
so it has no native genomic-coordinate columns. The predicted proteome
passed as `proteome_fasta` is expected to carry its gene's genomic location
on the FASTA defline after the id, as `contig:start-end:strand`
(a convention already used by several gene-prediction pipelines' protein
FASTA output). Requesting diamond's `qtitle` field (everything on the query
defline after the id) recovers that location string, which we then split
apart ourselves. The outfmt is therefore:

    6 qseqid sseqid pident qtitle

Window-restriction (genomic fallback): rather than relying on an exonerate
flag to restrict the search region (exonerate has no first-class "search
only this coordinate range of this multi-contig file" option), we actually
slice the requested contig region out of `genome_fasta` into a small
temporary FASTA file and pass *that* as `--target`. Exonerate's reported
hit coordinates are then local to the slice (1-based from the slice start),
so we add the window's start offset back on to recover real genome
coordinates before returning `SearchHit`s.
"""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from Bio import SeqIO

from MATPredict.detect.family_registry import Family, FamilyKey

_DIAMOND_OUTFMT = ["6", "qseqid", "sseqid", "pident", "qtitle"]


@dataclass(frozen=True)
class SearchHit:
    family_key: FamilyKey
    gene_name: str
    role: str  # core_MAT | flanking_conserved | flanking_variable
    contig: str
    start: int  # 1-based
    end: int  # 1-based, inclusive
    strand: str  # "+" | "-"
    identity: float
    reference_record_id: str  # which curated record's protein this matched
    method: str  # "diamond_proteome" | "exonerate_genome" | "exonerate_genome_relaxed"


def _gene_role_lookup(families: list[Family]) -> dict[str, tuple[FamilyKey, str]]:
    """Map a bare gene name to (family_key, role) for every family's expected genes."""
    lookup: dict[str, tuple[FamilyKey, str]] = {}
    for family in families:
        for gene in family.genes:
            lookup[gene["name"]] = (family.key, gene["role"])
    return lookup


def _parse_reference_header(sseqid: str) -> tuple[str, str]:
    """`{record_id}|gene{N}|{gene_name}` -> (record_id, gene_name)."""
    record_id, _, gene_name = sseqid.split("|")
    return record_id, gene_name


def search_fast_path(
    proteome_fasta: Path,
    families: list[Family],
    reference_fasta: Path,
    runner: Callable = subprocess.run,
) -> list[SearchHit]:
    """Search an existing predicted proteome against curated reference proteins via diamond blastp."""
    lookup = _gene_role_lookup(families)
    cmd = [
        "diamond", "blastp",
        "--query", str(proteome_fasta),
        "--db", str(reference_fasta),
        "--outfmt", *_DIAMOND_OUTFMT,
    ]
    result = runner(cmd, capture_output=True, text=True)

    hits: list[SearchHit] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        _qseqid, sseqid, pident, qtitle = line.split("\t")
        record_id, gene_name = _parse_reference_header(sseqid)
        if gene_name not in lookup:
            continue
        family_key, role = lookup[gene_name]
        contig, coords, strand = qtitle.split(":")
        start_str, end_str = coords.split("-")
        hits.append(
            SearchHit(
                family_key=family_key,
                gene_name=gene_name,
                role=role,
                contig=contig,
                start=int(start_str),
                end=int(end_str),
                strand=strand,
                identity=float(pident),
                reference_record_id=record_id,
                method="diamond_proteome",
            )
        )
    return hits


def _extract_window(genome_fasta: Path, window: tuple[str, int, int], tmp_dir: Path) -> Path:
    """Slice (contig, start, end) (1-based, inclusive) out of genome_fasta into its own FASTA file.

    The slice keeps the original contig name as its record id, so exonerate's
    GFF output still reports the real contig name -- only start/end need
    re-basing back onto genome coordinates afterward.
    """
    contig, start, end = window
    index = SeqIO.index(str(genome_fasta), "fasta")
    try:
        record = index[contig]
    finally:
        index.close()
    sliced = record[start - 1 : end]
    sliced.id = contig
    sliced.description = ""
    window_fasta = tmp_dir / f"window_{contig}_{start}_{end}.fasta"
    SeqIO.write(sliced, window_fasta, "fasta")
    return window_fasta


def search_genomic(
    genome_fasta: Path,
    families: list[Family],
    reference_fasta: Path,
    relaxed: bool = False,
    window: tuple[str, int, int] | None = None,
    runner: Callable = subprocess.run,
) -> list[SearchHit]:
    """Spliced protein-to-genome search via exonerate --model protein2genome.

    `window` restricts the search to (contig, start, end) (1-based, inclusive)
    -- used for the flanking-anchored second pass -- by slicing that region
    out of `genome_fasta` into a temporary FASTA and searching against just
    that slice; hit coordinates are re-based back onto genome coordinates
    before being returned. `relaxed=True` loosens exonerate's scoring
    threshold for that same second pass.
    """
    lookup = _gene_role_lookup(families)

    with tempfile.TemporaryDirectory() as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        offset = 0
        target_fasta = genome_fasta
        if window is not None:
            target_fasta = _extract_window(genome_fasta, window, tmp_dir)
            offset = window[1] - 1  # window start is 1-based; local coords start at 1

        cmd = [
            "exonerate", "--model", "protein2genome",
            "--query", str(reference_fasta),
            "--target", str(target_fasta),
            "--showtargetgff", "yes",
            "--showalignment", "no",
        ]
        if relaxed:
            cmd += ["--percent", "50"]

        result = runner(cmd, capture_output=True, text=True)

        hits: list[SearchHit] = []
        for line in result.stdout.splitlines():
            if "\tgene\t" not in line:
                continue
            fields = line.split("\t")
            contig, _src, _feat, start, end, _score, strand, _frame, attrs = fields
            # exonerate GFF attrs carry the query id under "sequence <id>"
            query_id = next(
                part.split(" ")[1] for part in attrs.split(" ; ") if part.startswith("sequence ")
            )
            record_id, gene_name = _parse_reference_header(query_id)
            if gene_name not in lookup:
                continue
            family_key, role = lookup[gene_name]
            method = "exonerate_genome_relaxed" if relaxed else "exonerate_genome"
            hits.append(
                SearchHit(
                    family_key=family_key,
                    gene_name=gene_name,
                    role=role,
                    contig=contig,
                    start=int(start) + offset,
                    end=int(end) + offset,
                    strand=strand,
                    identity=0.0,
                    reference_record_id=record_id,
                    method=method,
                )
            )
        return hits
