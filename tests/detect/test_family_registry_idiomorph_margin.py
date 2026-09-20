# tests/detect/test_family_registry_idiomorph_margin.py
"""`min_idiomorph_margin` is per-locus curation data, like `max_cluster_gap_bp`.

How far apart two idiomorphs' identities must be before the call is trusted
depends on how densely that phylum's idiomorphs are represented in the curated
database, and that differs sharply between phyla -- Mucoromycota has 6 records
against Ascomycota's many. It is therefore a property of the locus entry, not
a global tuning constant.

Measured basis for the 5.0 default: across the 23 ground-truth Mucoromycota
genomes the identity rule is 23/23 correct, but the margins split hard by
idiomorph. Plus calls are separated by 44.6-64.2 points; Minus calls by only
2.3-13.0. The narrowest call seen anywhere, Blakeslea trispora at 0.59, is
outside the truth set so its correctness is unknown. A 5.0 bar caps the four
calls under five points and leaves the well-separated ones at full tier.
"""
from __future__ import annotations

from pathlib import Path

from MATPredict.detect.family_registry import (
    DEFAULT_MIN_IDIOMORPH_MARGIN,
    load_all_families,
)

_ORDER_YML = """\
phylum: TestPhylum
loci:
  - locus_name: "Declared"
    vocabulary_type: "enum"
    idiomorph_values: ["Plus", "Minus"]
    taxonomic_scope: [1]
    min_idiomorph_margin: 12.5
    genes:
      - {name: g1, role: core_MAT}
  - locus_name: "Undeclared"
    vocabulary_type: "enum"
    idiomorph_values: ["Plus", "Minus"]
    taxonomic_scope: [1]
    genes:
      - {name: g1, role: core_MAT}
"""


def _db(tmp_path: Path) -> Path:
    phylum = tmp_path / "TestPhylum"
    phylum.mkdir()
    (phylum / "order.yml").write_text(_ORDER_YML)
    return tmp_path


def test_a_declared_margin_is_loaded_onto_the_family(tmp_path):
    families = {f.key.locus_name: f for f in load_all_families(_db(tmp_path))}
    assert families["Declared"].min_idiomorph_margin == 12.5


def test_a_locus_that_declares_nothing_gets_the_default(tmp_path):
    # Absence means "the default is fine here", exactly as it does for
    # max_cluster_gap_bp -- not that every locus must restate the value.
    families = {f.key.locus_name: f for f in load_all_families(_db(tmp_path))}
    assert families["Undeclared"].min_idiomorph_margin == DEFAULT_MIN_IDIOMORPH_MARGIN
    assert DEFAULT_MIN_IDIOMORPH_MARGIN == 5.0


def test_the_real_mucoromycota_locus_declares_a_margin():
    # The curator ruled 5.0 for this locus on 2026-09-20. Asserting it here
    # rather than only in order.yml means a silent edit to the curation data
    # cannot pass unnoticed.
    families = {
        f.key: f for f in load_all_families(Path(__file__).parents[2] / "db")
    }
    mucoro = next(f for k, f in families.items() if k.phylum == "Mucoromycota")
    assert mucoro.min_idiomorph_margin == 5.0
