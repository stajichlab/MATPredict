"""scripts/label_loci_table.py: the `nearest_reference` column (#200, curator 2026-11-07)."""
import csv
import importlib.util
import sys
from pathlib import Path

import pytest

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT / "scripts"))
_spec = importlib.util.spec_from_file_location("label_loci_table", _ROOT / "scripts" / "label_loci_table.py")
llt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(llt)


def _record(db, phylum, order, record_id, cls):
    d = db / phylum / order / record_id
    d.mkdir(parents=True)
    (d / "metadata.yaml").write_text(
        f"record_id: {record_id}\ntaxonomy:\n  taxid: 1\n  lineage: k__Fungi;p__{phylum};c__{cls};o__X;f__Y;g__Z\n")


@pytest.fixture()
def db(tmp_path):
    root = tmp_path / "db"
    _record(root, "Ascomycota", "Sordariales", "100_strain_MAT_MAT1-1", "Sordariomycetes")
    _record(root, "Ascomycota", "Eurotiales", "200_strain_MAT_MAT1-2", "Eurotiomycetes")
    _record(root, "Ascomycota", "Orbiliales", "300_strain_MAT_MAT1-1", "Orbiliomycetes")
    return str(root)


def _row(core_reference=None, **kw):
    row = dict(genome="G", class_="Sordariomycetes", family_called="Ascomycota:MAT", idiomorph="MAT1-1", genes_found="MAT1-1-1|MAT1-2-1")
    row.update(kw)
    if core_reference is not None:
        row["core_reference"] = core_reference
    return row


def test_both_core_genes_are_listed_with_the_system_of_their_matched_record(db):
    assert llt.nearest_reference("MAT1-1-1:100_strain_MAT_MAT1-1:55.49:120|MAT1-2-1:200_strain_MAT_MAT1-2:40:80", db) == "MAT1-1-1=NeuCra|MAT1-2-1=AspFum"


def test_a_class_without_a_registry_code_is_written_as_the_record_id(db):
    assert llt.nearest_reference("MAT1-1-1:300_strain_MAT_MAT1-1:50:90", db) == "MAT1-1-1=300_strain_MAT_MAT1-1"


def test_an_unknown_record_is_written_as_its_record_id(db):
    assert llt.nearest_reference("MAT1-1-1:no_such_record:50:90", db) == "MAT1-1-1=no_such_record"


@pytest.mark.parametrize("value", ["", None])
def test_no_matched_record_is_the_explicit_none(db, value):
    assert llt.nearest_reference(value, db) == "none"


def test_every_row_gets_a_value_and_none_is_never_blank(db):
    rows = llt.label_rows([_row("MAT1-1-1:100_strain_MAT_MAT1-1:55:1"), _row(""), _row("MAT1-2-1:300_strain_MAT_MAT1-1:50:1", class_="Orbiliomycetes")], db)
    assert [r["nearest_reference"] for r in rows] == ["MAT1-1-1=NeuCra", "none", "MAT1-2-1=300_strain_MAT_MAT1-1"]
    assert all(r["nearest_reference"] for r in rows)


def test_a_table_without_the_core_reference_column_gets_a_blank_value_not_none(db):
    (r,) = llt.label_rows([_row()], db)
    assert r["nearest_reference"] == "" and r["nearest_reference_detail"] == ""


def test_the_detail_column_keeps_the_record_identity_and_bitscore(db):
    ref = "MAT1-1-1:100_strain_MAT_MAT1-1:55.49:120"
    (r,) = llt.label_rows([_row(ref)], db)
    assert r["nearest_reference_detail"] == ref


def test_identical_input_gives_identical_output_and_the_raw_columns_are_untouched(db, tmp_path):
    src = tmp_path / "loci.tsv"
    rows = [_row("MAT1-1-1:100_strain_MAT_MAT1-1:55:1|MAT1-2-1:200_strain_MAT_MAT1-2:40:1"), _row("")]
    with open(src, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    outs = []
    for name in ("a.tsv", "b.tsv"):
        out = tmp_path / name
        sys.argv = ["label_loci_table.py", "--loci", str(src), "--out", str(out), "--db", db]
        llt.main()
        outs.append(out.read_text())
    assert outs[0] == outs[1]
    got = list(csv.DictReader(outs[0].splitlines(), delimiter="\t"))
    assert all(got[i][k] == rows[i][k] for i in range(2) for k in rows[i])


def test_the_real_registry_maps_the_six_pezizomycotina_classes():
    systems = llt.L.registry()["pezizomycotina_class_systems"]
    assert set(systems.values()) == {"NeuCra", "AspFum", "BipMay", "BotCin", "XanPar", "TubMel"}
