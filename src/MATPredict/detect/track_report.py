"""Gene-model tracks for the report: the homology models, helixer-lite on the locus windows, and a GenBank
annotation when one was supplied, compared gene by gene.

Everything here is report-only: tracks never change a call. `build_tracks_doc` turns a finished
`detection_report.yaml` document plus the optional sources into `tracks.yaml`, which the HTML report reads.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import yaml

from MATPredict.detect import tracks as T
from MATPredict.detect.helixer_windows import DEFAULT_FLANK, merged_windows, run_helixer_windows

TRACKS_YAML = "tracks.yaml"
#: Two ends count as agreeing when they differ by at most this many bp (the panel metric's tolerance).
END_TOLERANCE_BP = 30


def read_contigs(fasta: Path, wanted: set[str]) -> tuple[dict[str, str], dict[str, int]]:
    """(sequences of the wanted contigs, lengths of every contig), in one pass over a plain FASTA."""
    from Bio import SeqIO

    seqs, lengths = {}, {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        lengths[rec.id] = len(rec.seq)
        if rec.id in wanted:
            seqs[rec.id] = str(rec.seq).upper()
    return seqs, lengths


def locus_spans(doc: dict) -> list[tuple[str, int, int]]:
    """(contig, start, end) of every called locus."""
    return [(d["contig"], int(d["start"]), int(d["end"])) for d in doc.get("detected") or []
            if d.get("contig") and d.get("start") is not None and d.get("end") is not None]


def evidence_track(doc: dict, sequences: dict[str, str]) -> T.Track:
    """The homology models of a detection report: each located gene that did not lose a pair resolution."""
    genes, seen = [], set()
    for call in doc.get("detected") or []:
        for g in call.get("gene_evidence") or []:
            if g.get("lost_to") or g.get("start") is None or g.get("contig") not in sequences:
                continue
            key = (g["contig"], g["start"], g["end"], g["strand"], g["gene"])
            if key in seen:
                continue
            seen.add(key)
            exons = tuple(sorted((int(x["start"]), int(x["end"])) for x in g["exons"])) if g.get("exons") \
                else ((int(g["start"]), int(g["end"])),)
            gene = T.TrackGene(name=g["gene"], contig=g["contig"], start=exons[0][0], end=exons[-1][1],
                               strand=g["strand"], exons=exons, source="evidence")
            genes.append(T._with_orf(gene, sequences[g["contig"]]))
    return T.Track(source="evidence", genes=tuple(genes))


def restrict(track: T.Track, windows) -> T.Track:
    """The models of `track` that overlap any (contig, start, end) window."""
    keep = tuple(g for g in track.genes
                 if any(g.contig == c and g.end >= s and g.start <= e for c, s, e in windows))
    return T.Track(source=track.source, genes=keep)


def _gene_doc(g: T.TrackGene) -> dict:
    d = asdict(g)
    d["exons"] = [list(x) for x in g.exons]
    return d


def _summary_doc(rows: list[T.GeneSummary], evidence: T.Track) -> list[dict]:
    """One entry per evidence gene: its best same-strand model (largest overlap), or None."""
    out = []
    for e in evidence.genes:
        mine = [r for r in rows if r.evidence_gene == e.name and r.strand_agrees]
        best = max(mine, key=lambda r: r.overlap_bp, default=None)
        entry = {"gene": e.name, "contig": e.contig, "start": e.start, "end": e.end, "strand": e.strand,
                 "model": None}
        if best:
            entry["model"] = {**asdict(best),
                              "ends_agree": abs(best.extends_5p) <= END_TOLERANCE_BP
                              and abs(best.extends_3p) <= END_TOLERANCE_BP}
        out.append(entry)
    return out


def build_tracks_doc(doc: dict, genome_fasta: Path, *, helixer_cmd: list[str] | None = None,
                     helixer_lineage: str = "fungi", helixer_cpus: int = 4, flank: int = DEFAULT_FLANK,
                     genbank: Path | None = None) -> dict | None:
    """The `tracks.yaml` document, or None when there is nothing to compare (no loci, or no other source).

    `genbank` adds the annotation track; its records are matched to the genome FASTA by record id."""
    spans = locus_spans(doc)
    if not spans or (helixer_cmd is None and genbank is None):
        return None
    wanted = {c for c, _, _ in spans}
    seqs, lengths = read_contigs(genome_fasta, wanted)
    windows = merged_windows(spans, flank, lengths)
    evidence = evidence_track(doc, seqs)
    sources = {"evidence": evidence}
    notes = []
    if helixer_cmd is not None:
        sources["helixer"] = run_helixer_windows(spans, seqs, lengths, flank=flank, lineage=helixer_lineage,
                                                 cpus=helixer_cpus, binary=helixer_cmd)
    if genbank is not None:
        gb = T.from_genbank(genbank, seqs)
        missing = {g.contig for g in gb.genes} - set(lengths)
        if missing:
            notes.append("GenBank records not in the genome FASTA: " + ", ".join(sorted(missing)))
        sources["genbank"] = restrict(gb, windows)
    return {
        "flank_bp": flank,
        "windows": [{"contig": c, "start": s, "end": e} for c, s, e in windows],
        "end_tolerance_bp": END_TOLERANCE_BP,
        "notes": notes,
        "sources": {name: {"n_models": len(t.genes), "genes": [_gene_doc(g) for g in t.genes]}
                    for name, t in sources.items()},
        "summaries": {name: _summary_doc(T.summary(evidence, t), evidence)
                      for name, t in sources.items() if name != "evidence"},
    }


def write_tracks(path: Path, tracks_doc: dict) -> None:
    Path(path).write_text(yaml.safe_dump(tracks_doc, sort_keys=False))


def load_tracks(run_dir: Path) -> dict | None:
    p = Path(run_dir) / TRACKS_YAML
    return yaml.safe_load(p.read_text()) if p.is_file() else None
