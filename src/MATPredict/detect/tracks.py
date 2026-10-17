"""Annotation tracks: gene models from sources other than the search.

A track is a list of gene models on one locus window. Three sources feed the report:

* ``evidence``: the pipeline's own gene models (homology search, then exonerate and miniprot).
* ``helixer``: helixer-lite predictions (ab initio; no evidence input), reported alongside.
* ``genbank``: annotation supplied by the user as a GenBank file (CDS and gene features).

Each model carries its exons, strand, and whether its ORF is complete (ATG start, stop codon at the
end, no internal stop). `summary` compares a track with the evidence track: which evidence genes a
model overlaps, and how far each end of the model extends past the evidence gene.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

STOPS = {"TAA", "TAG", "TGA"}
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class TrackGene:
    name: str
    contig: str
    start: int  # 1-based inclusive, genomic
    end: int
    strand: str  # "+" or "-"
    exons: tuple[tuple[int, int], ...]
    source: str
    protein_length: int | None = None
    start_codon: str | None = None
    stop_codon: str | None = None
    complete_orf: bool | None = None


@dataclass(frozen=True)
class Track:
    source: str
    genes: tuple[TrackGene, ...] = ()


def _revcomp(dna: str) -> str:
    return dna.translate(_COMP)[::-1]


def coding_sequence(gene_exons, strand, contig_seq: str) -> str:
    """The spliced coding sequence of a model, reverse-complemented on the minus strand."""
    dna = "".join(contig_seq[a - 1:b] for a, b in sorted(gene_exons))
    return _revcomp(dna) if strand == "-" else dna


def orf_status(dna: str) -> tuple[str, str, int, bool]:
    """(start codon, stop codon, protein length, complete ORF) for a coding sequence."""
    dna = dna.upper()
    start, stop = dna[:3], dna[-3:]
    body = dna[: len(dna) // 3 * 3]
    codons = [body[i:i + 3] for i in range(0, len(body), 3)]
    internal = any(c in STOPS for c in codons[:-1])
    complete = start == "ATG" and stop in STOPS and not internal and len(dna) % 3 == 0
    return start, stop, len(codons), complete


def _with_orf(gene: TrackGene, contig_seq: str) -> TrackGene:
    dna = coding_sequence(gene.exons, gene.strand, contig_seq)
    start, stop, length, complete = orf_status(dna)
    return TrackGene(**{**gene.__dict__, "protein_length": length, "start_codon": start,
                        "stop_codon": stop, "complete_orf": complete})


def from_helixer_gff3(path: Path, sequences: dict[str, str], window_offsets: bool = False) -> Track:
    """Gene models from helixer-lite's GFF3 (gene, mRNA and CDS lines; exons taken from CDS).

    With window_offsets, the sequence names are windows written as ``contig:start-end`` (the FASTA headers
    of helixer_windows); coordinates are shifted by start-1 and the contig is the part before the colon."""
    genes: dict[str, dict] = {}
    mrna_parent: dict[str, str] = {}
    for line in Path(path).read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        contig, kind, start, end, strand, attrs = f[0], f[2], int(f[3]), int(f[4]), f[6], f[8]
        if window_offsets:
            m = re.match(r"^(.*):(\d+)-(\d+)$", contig)
            if m:
                contig, offset = m.group(1), int(m.group(2)) - 1
                start, end = start + offset, end + offset
        ident = re.search(r"ID=([^;]+)", attrs)
        parent = re.search(r"Parent=([^;]+)", attrs)
        if kind == "gene" and ident:
            genes[ident.group(1)] = {"contig": contig, "start": start, "end": end, "strand": strand,
                                     "exons": []}
        elif kind == "mRNA" and ident and parent:
            mrna_parent[ident.group(1)] = parent.group(1)
        elif kind == "CDS" and parent:
            gid = mrna_parent.get(parent.group(1), parent.group(1))
            if gid in genes:
                genes[gid]["exons"].append((start, end))
    out = []
    for gid, g in genes.items():
        if not g["exons"]:
            continue
        exons = tuple(sorted(g["exons"]))
        gene = TrackGene(name=gid, contig=g["contig"], start=exons[0][0], end=exons[-1][1],
                         strand=g["strand"], exons=exons, source="helixer")
        out.append(_with_orf(gene, sequences[g["contig"]]))
    return Track(source="helixer", genes=tuple(out))


def from_genbank(path: Path, sequences: dict[str, str] | None = None) -> Track:
    """Gene models from a GenBank file: each CDS feature (joined exons), named by /gene, /locus_tag or
    /product. A gene without a CDS is kept only when it has no CDS of its own under the same name."""
    from Bio import SeqIO

    out = []
    for rec in SeqIO.parse(str(path), "genbank"):
        seq = str(rec.seq).upper() if sequences is None else sequences.get(rec.id, str(rec.seq)).upper()
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            name = (q.get("gene") or q.get("locus_tag") or q.get("product") or ["CDS"])[0]
            exons = tuple(sorted((int(p.start) + 1, int(p.end)) for p in feat.location.parts))
            strand = "-" if feat.location.strand == -1 else "+"
            gene = TrackGene(name=str(name), contig=rec.id, start=exons[0][0], end=exons[-1][1],
                             strand=strand, exons=exons, source="genbank")
            out.append(_with_orf(gene, seq))
    return Track(source="genbank", genes=tuple(out))


@dataclass(frozen=True)
class GeneSummary:
    evidence_gene: str
    other_gene: str
    strand_agrees: bool
    overlap_bp: int
    extends_5p: int  # bases the other model extends past the evidence gene at its 5' end (negative: shorter)
    extends_3p: int
    other_complete_orf: bool | None


def summary(evidence: Track, other: Track) -> list[GeneSummary]:
    """For each evidence gene, the overlapping model(s) of the other track: how far each end extends.

    5' and 3' follow the gene's strand. A model on the other strand overlapping the gene is reported with
    strand_agrees False and no extension (it is a different gene on the other strand)."""
    rows = []
    for e in evidence.genes:
        for o in other.genes:
            if o.contig != e.contig or o.end < e.start or o.start > e.end:
                continue
            overlap = min(o.end, e.end) - max(o.start, e.start) + 1
            agrees = o.strand == e.strand
            if agrees and e.strand == "+":
                ext5, ext3 = e.start - o.start, o.end - e.end
            elif agrees:
                ext5, ext3 = o.end - e.end, e.start - o.start
            else:
                ext5 = ext3 = 0
            rows.append(GeneSummary(e.name, o.name, agrees, overlap,
                                    ext5 if agrees else 0, ext3 if agrees else 0, o.complete_orf))
    return rows
