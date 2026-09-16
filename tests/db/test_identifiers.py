from __future__ import annotations

from MATPredict.db.identifiers import build_record_id, slugify_strain


def test_slugify_known_strain():
    assert slugify_strain("NRRL 1555", known=True, existing_slugs=set()) == "nrrl-1555"


def test_slugify_unknown_strain_gets_counter():
    existing = set()
    first = slugify_strain(None, known=False, existing_slugs=existing)
    existing.add(first)
    second = slugify_strain(None, known=False, existing_slugs=existing)
    assert first == "unknown-1"
    assert second == "unknown-2"


def test_build_record_id():
    rid = build_record_id(taxid=4837, strain_slug="nrrl-1555", locus_name="MAT", idiomorph_key="Plus")
    assert rid == "4837_nrrl-1555_MAT_Plus"
