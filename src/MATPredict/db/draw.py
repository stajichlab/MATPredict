"""Render a single-locus gene-structure diagram for one curated record, via pyGenomeViz.

This is a curated-record-side (`db/`) feature, not a detection-side one: it renders
`gff_export.write_genbank`'s real `locus.gbk` output (Task 1), which carries real
segment sequence and per-gene `CDS` features with `role`/`gene_class`/
`present_in_idiomorphs` qualifiers already attached (see `gff_export.write_genbank`'s
own docstring). See `docs/notes/2026-09-18_mat-locus-visualization-research.md` section
7 for the real, hands-on trial (Task 3) this module's design is based on -- pyGenomeViz
v1.7.0's real, installed API is class-based (`pygenomeviz.GenomeViz`,
`pygenomeviz.parser.Genbank`, `FeatureTrack`/`FeatureSegment`), confirmed by reading the
actual installed package source, not guessed from older functional-API examples online.

Scope (deliberate, disclosed limitation -- see section 7.4 of the research doc): only
single-segment (non-fragmented) records are supported. A fragmented record (a MAT locus
split across multiple scaffolds/contigs, one `SeqRecord` per segment in the GenBank
file) is explicitly detected and rejected with `NotImplementedError` rather than
attempting a wrong or partial render -- pyGenomeViz's GFF3 parser has two real,
confirmed problems with MATPredict's fragmented-multi-contig GFF3 shape (silent
multi-seqid truncation, and a hard `FeatureRangeError` from non-1-based
`##sequence-region` windows), and the GenBank path used here has no direct multi-contig
"stitch these into one figure" API either -- that is future work, not attempted here.
"""
from __future__ import annotations

from pathlib import Path

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank

# A simple, distinct, colorblind-safe 3-color palette (the "Okabe-Ito" set commonly used
# for exactly this purpose), one color per `role` qualifier that gff_export.write_genbank
# already attaches to every gene/CDS feature:
#   core_MAT            -> vermillion (#D55E00) -- the mating-type-determining genes themselves
#   flanking_conserved  -> blue       (#0072B2) -- genes conserved across idiomorphs/species
#   flanking_variable   -> green      (#009E73) -- genes present in only some idiomorphs
# Any gene whose `role` qualifier is missing or not one of these three (should not happen
# for a real curated record, but never silently mis-color one) falls back to a neutral grey.
ROLE_COLORS: dict[str, str] = {
    "core_MAT": "#D55E00",
    "flanking_conserved": "#0072B2",
    "flanking_variable": "#009E73",
}
_FALLBACK_COLOR = "#999999"


def _gene_label(feature) -> str:
    """Build a display label from `gene`/`gene_class`/`present_in_idiomorphs`
    qualifiers, whichever are present -- never fabricating a value that isn't already
    a real qualifier on the feature."""
    gene_name = feature.qualifiers.get("gene", [None])[0]
    gene_class = feature.qualifiers.get("gene_class", [None])[0]
    idiomorphs = feature.qualifiers.get("present_in_idiomorphs", [None])[0]
    parts = [part for part in (gene_name, gene_class, idiomorphs) if part]
    return "/".join(parts)


def draw_locus(gbk_path: Path, out_path: Path) -> Path:
    """Render one static image of a single-segment curated record's real, curated gene
    structure (from `gbk_path`, a `locus.gbk` file produced by `gff_export.write_genbank`)
    to `out_path`, colored by `role` and labeled by `gene`/`gene_class`/
    `present_in_idiomorphs`.

    Output format is whatever `out_path`'s extension implies (matplotlib's own
    `savefig`-driven format inference, via pyGenomeViz's `GenomeViz.savefig`) -- `.png`
    and `.svg` both work; `.png` is the default recommendation for a quick static figure.

    One `CDS` feature per gene is drawn with `FeatureSegment.add_exon_features`, which
    natively handles both a single-exon gene (a plain `FeatureLocation`) and a real
    multi-exon gene (a `CompoundLocation`, one part per exon, built by
    `gff_export._cds_location`) -- drawing separate exon boxes joined by intron lines
    for the latter, rather than one unbroken arrow spanning intronic sequence. A gene
    with no available protein sequence (so `write_genbank` emitted only a `gene`
    feature, no `CDS`) is drawn from its `gene` feature instead, so it is not silently
    omitted from the figure.

    Raises `NotImplementedError` for a fragmented (multi-segment/multi-contig) record --
    see this module's docstring for why that case is out of scope for now.
    """
    gbk_path = Path(gbk_path)
    out_path = Path(out_path)

    gbk = Genbank(gbk_path)
    if len(gbk.records) != 1:
        raise NotImplementedError(
            f"draw_locus does not support a fragmented (multi-segment/multi-contig) "
            f"record: {gbk_path} contains {len(gbk.records)} segments/contigs "
            f"({[rec.id for rec in gbk.records]}). Single-segment records only, per "
            "this task's disclosed scope limitation -- see docs/notes/"
            "2026-09-18_mat-locus-visualization-research.md section 7.4."
        )

    seq_record = gbk.records[0]

    # One feature per gene, preferring the real CDS feature (which carries the real
    # exon/intron CompoundLocation when the gene is multi-exon) over the plain `gene`
    # feature, but falling back to the `gene` feature for a gene with no CDS (no
    # available protein sequence) so it is still shown rather than silently dropped.
    feature_by_gene: dict[str, object] = {}
    for feature in seq_record.features:
        if feature.type not in ("gene", "CDS"):
            continue
        gene_name = feature.qualifiers.get("gene", [None])[0]
        if gene_name is None:
            continue
        if feature.type == "CDS" or gene_name not in feature_by_gene:
            feature_by_gene[gene_name] = feature

    gv = GenomeViz(fig_width=12, fig_track_height=1.5, show_axis=True)
    track = gv.add_feature_track(seq_record.id, len(seq_record.seq))
    segment = track.get_segment()

    for feature in feature_by_gene.values():
        role = feature.qualifiers.get("role", [None])[0]
        color = ROLE_COLORS.get(role, _FALLBACK_COLOR)
        # add_exon_features reads its label straight from feature.qualifiers[label_type][0];
        # stash the composite label under a private key rather than fabricating a new
        # feature type/qualifier convention elsewhere in the codebase.
        feature.qualifiers["_draw_label"] = [_gene_label(feature)]
        segment.add_exon_features(
            feature,
            plotstyle="bigarrow",
            label_type="_draw_label",
            patch_kws={"fc": color, "ec": "black", "lw": 0.5},
        )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    gv.savefig(out_path)
    return out_path
