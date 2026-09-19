from __future__ import annotations

import pytest

from MATPredict.db.draw import draw_locus
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


def test_draw_locus_rejects_fragmented_multi_segment_record(tmp_path):
    gbk_path = tmp_path / "locus.gbk"
    write_genbank(FRAGMENTED_RECORD, sequences={0: "M" * 10, 1: "M" * 10}, out_path=gbk_path)

    out_path = tmp_path / "locus.png"
    with pytest.raises(NotImplementedError):
        draw_locus(gbk_path, out_path)
    assert not out_path.exists()
