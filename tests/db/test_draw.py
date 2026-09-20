from __future__ import annotations

import pytest
from pygenomeviz.parser import Genbank

from MATPredict.db.draw import _gene_label, draw_locus
from MATPredict.db.gff_export import write_genbank

# A real single-segment record with one multi-exon minus-strand gene (matching the real
# COX13-shaped curated-record convention: exons stored in descending genomic order) and
# one single-exon plus-strand gene with a gene_class/present_in_idiomorphs label, one of
# each `role` this project uses -- enough to exercise the color palette and the
# CompoundLocation exon-drawing path together, without inventing any exon/gene data
# beyond what a real curated record's schema shape already looks like.
SINGLE_SEGMENT_RECORD = {
    "record_id": "555_e_MAT_combined",
    "locus": {"core": {"segments": [
        {"segment_index": 0, "start": 1, "end": 200, "sequence_source": {"seq_region": "scaffold_1"}},
    ]}},
    "genes": [
        {"gene_index": 0, "name": "APN2like", "role": "flanking_conserved", "present": True,
         "segment_index": 0, "start": 10, "end": 100, "strand": "-",
         "exons": [{"start": 80, "end": 100}, {"start": 40, "end": 60}, {"start": 10, "end": 25}],
         "gene_class": "flanking"},
        {"gene_index": 1, "name": "MAT1-1-1", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 110, "end": 180, "strand": "+",
         "present_in_idiomorphs": ["MAT1-1"]},
    ],
}

# Real B-locus biology, matching the real curated record
# db/Basidiomycota/Agaricales/5346_a43-b43-okayama-7_PR_B43/locus.gbk's shape: several
# DISTINCT genes sharing one `gene` qualifier value on purpose (pheromone-receptor gene
# duplication/multi-allele co-occurrence is normal MAT biology, not something to
# collapse -- see this project's own domain-knowledge memory note on this). A synthetic
# fixture is used (rather than reading that real file directly) so this regression test
# is not fragile to a future re-generation of that record's real locus.gbk.
DUPLICATE_GENE_NAME_RECORD = {
    "record_id": "777_g_PR_dup",
    "locus": {"core": {"segments": [
        {"segment_index": 0, "start": 1, "end": 21000, "sequence_source": {"seq_region": "scaffold_1"}},
    ]}},
    "genes": [
        {"gene_index": 0, "name": "pheromone_receptor", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 1, "end": 2612, "strand": "-"},
        {"gene_index": 1, "name": "pheromone_B44", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 5798, "end": 6523, "strand": "-"},
        {"gene_index": 2, "name": "pheromone_B43", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 8739, "end": 8933, "strand": "+"},
        {"gene_index": 3, "name": "pheromone_receptor", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 11171, "end": 13284, "strand": "-"},
        {"gene_index": 4, "name": "pheromone_receptor", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 13490, "end": 15124, "strand": "-"},
        {"gene_index": 5, "name": "pheromone_receptor", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 15801, "end": 17793, "strand": "+"},
        {"gene_index": 6, "name": "fungal_mating_type_pheromone", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 18817, "end": 19240, "strand": "+"},
        {"gene_index": 7, "name": "fungal_mating_type_pheromone", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 20032, "end": 20706, "strand": "+"},
    ],
}

FRAGMENTED_RECORD = {
    "record_id": "666_f_MAT_combined",
    "locus": {"core": {"segments": [
        {"segment_index": 0, "start": 1, "end": 100, "sequence_source": {"seq_region": "scaffold_a"}},
        {"segment_index": 1, "start": 1, "end": 100, "sequence_source": {"seq_region": "scaffold_b"}},
    ]}},
    "genes": [
        {"gene_index": 0, "name": "G1", "role": "core_MAT", "present": True,
         "segment_index": 0, "start": 10, "end": 50, "strand": "+"},
        {"gene_index": 1, "name": "G2", "role": "core_MAT", "present": True,
         "segment_index": 1, "start": 10, "end": 50, "strand": "+"},
    ],
}


def test_draw_locus_produces_a_real_nonempty_png(tmp_path):
    gbk_path = tmp_path / "locus.gbk"
    write_genbank(SINGLE_SEGMENT_RECORD, sequences={0: "M" * 20, 1: "M" * 15}, out_path=gbk_path)

    out_path = tmp_path / "locus.png"
    result = draw_locus(gbk_path, out_path)

    assert result == out_path
    assert out_path.exists()
    assert out_path.stat().st_size > 0
    # Confirm it is a genuinely well-formed PNG, not just a non-empty file -- check the
    # real PNG magic bytes rather than trusting the extension alone.
    with open(out_path, "rb") as handle:
        magic = handle.read(8)
    assert magic == b"\x89PNG\r\n\x1a\n"


def test_draw_locus_svg_output_is_well_formed(tmp_path):
    """pyGenomeViz/matplotlib infer the output format from the file extension; confirm
    the .svg path also produces real, parseable output, not just .png."""
    gbk_path = tmp_path / "locus.gbk"
    write_genbank(SINGLE_SEGMENT_RECORD, sequences={0: "M" * 20, 1: "M" * 15}, out_path=gbk_path)

    out_path = tmp_path / "locus.svg"
    draw_locus(gbk_path, out_path)

    content = out_path.read_text()
    assert content.strip().startswith("<?xml")
    assert "<svg" in content


def test_draw_locus_draws_every_gene_even_when_names_repeat(tmp_path, monkeypatch):
    """Regression test for the real bug this fixture models: `draw_locus` used to key
    features by their `gene` qualifier string in a dict, so a later feature with the
    same gene name silently overwrote an earlier one. `DUPLICATE_GENE_NAME_RECORD` has
    8 real gene features, 4 of which share the name "pheromone_receptor" -- normal
    B-locus biology (gene duplication/multi-allele co-occurrence), not a case to
    collapse. All 8 must be drawn, not just the 4 distinct names."""
    from pygenomeviz.segment.feature import FeatureSegment

    gbk_path = tmp_path / "locus.gbk"
    write_genbank(DUPLICATE_GENE_NAME_RECORD, sequences={}, out_path=gbk_path)

    # Confirm the fixture itself really has 8 gene features in the written GenBank
    # (i.e. the bug scenario is real, not accidentally collapsed before draw_locus
    # even runs).
    gbk = Genbank(gbk_path)
    assert len(gbk.records) == 1
    gene_features = [f for f in gbk.records[0].features if f.type == "gene"]
    assert len(gene_features) == 8

    calls = []
    real_add_exon_features = FeatureSegment.add_exon_features

    def spy_add_exon_features(self, features, **kwargs):
        calls.append(features)
        return real_add_exon_features(self, features, **kwargs)

    monkeypatch.setattr(FeatureSegment, "add_exon_features", spy_add_exon_features)

    out_path = tmp_path / "locus.png"
    draw_locus(gbk_path, out_path)

    # All 8 real gene features were handed to pyGenomeViz's real feature-drawing API --
    # none silently dropped because its name repeated.
    assert len(calls) == 8
    drawn_gene_names = sorted(f.qualifiers.get("gene", [None])[0] for f in calls)
    assert drawn_gene_names == sorted(gene["name"] for gene in DUPLICATE_GENE_NAME_RECORD["genes"])


def test_draw_locus_rejects_fragmented_multi_segment_record(tmp_path):
    gbk_path = tmp_path / "locus.gbk"
    write_genbank(FRAGMENTED_RECORD, sequences={0: "M" * 10, 1: "M" * 10}, out_path=gbk_path)

    out_path = tmp_path / "locus.png"
    with pytest.raises(NotImplementedError):
        draw_locus(gbk_path, out_path)
    assert not out_path.exists()


def test_gene_label_includes_the_gene_class_qualifier_written_by_write_genbank(tmp_path):
    """The two halves of the join, end to end at the label: `write_genbank` writes
    `/gene_class` from its `gene_classes` mapping, and `_gene_label` reads it back off
    the parsed CDS feature. Before the join existed the qualifier was never written, so
    every label was gene name alone."""
    gbk_path = tmp_path / "locus.gbk"
    write_genbank(
        SINGLE_SEGMENT_RECORD,
        sequences={0: "M" * 20},
        out_path=gbk_path,
        gene_classes={0: "apn2_homolog"},
    )

    cds = [
        feature
        for feature in Genbank(gbk_path).records[0].features
        if feature.type == "CDS"
    ]
    assert len(cds) == 1
    # Also pins the direction of the fix: this fixture's gene 0 still carries a stale
    # record-level "gene_class": "flanking" field, which the record schema does not
    # define and `write_genbank` no longer reads. The label shows the mapping's value,
    # not that field's.
    assert _gene_label(cds[0]) == "APN2like/apn2_homolog"
