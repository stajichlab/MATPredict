"""exonerate is given only the requested gene's reference proteins, not the whole reference set.

`_select_requested_gene` discards every alignment whose query is a different gene, and exonerate aligns each
query on its own, so the other genes' proteins cost time and change nothing. The subset keeps every record of
the gene, across families, so the selection among records of that gene is the one it always made."""
from __future__ import annotations

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import _gene_query_records, _write_gene_query, polish_with_exonerate

FAMILY = Family(FamilyKey("Ascomycota", "MAT"), "pattern", None, r"^MAT1-[12]$",
                [{"name": "MAT1-1-1", "role": "core_MAT"}, {"name": "APN2", "role": "flanking_conserved"}], [4890])

REFERENCE = (
    ">recA|gene0|MAT1-1-1 first\nMKLV\n"
    ">recA|gene1|APN2\nMSSSS\n"
    ">recB|gene0|MAT1-1-1\nMQQQ\n"
    ">recC|gene0|SLA2\nMTTTT\n"
    ">not-a-reference-header\nMXXX\n"
)

_GFF = ("c1\texonerate\tgene\t1\t400\t.\t+\t.\tgene_id 1 ; sequence recB|gene0|MAT1-1-1 ; gene_orientation . ; identity 90.00 ; similarity 91.00\n"
        "c1\texonerate\texon\t1\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n")


class _Result:
    def __init__(self, returncode=0, stdout=""):
        self.returncode, self.stdout, self.stderr = returncode, stdout, ""


def _setup(tmp_path):
    ref = tmp_path / "reference.faa"
    ref.write_text(REFERENCE)
    (tmp_path / "genome.fa").write_text(">c1\n" + "ACGT" * 150 + "\n")
    return ref


def test_the_query_holds_only_the_requested_genes_records_in_file_order(tmp_path):
    ref = _setup(tmp_path)
    out = _write_gene_query(ref, "MAT1-1-1", tmp_path)
    assert out.read_text() == ">recA|gene0|MAT1-1-1 first\nMKLV\n>recB|gene0|MAT1-1-1\nMQQQ\n"


def test_a_gene_with_no_reference_protein_has_no_query(tmp_path):
    assert _write_gene_query(_setup(tmp_path), "NOPE", tmp_path) is None


def test_a_header_that_is_not_a_reference_header_is_left_out(tmp_path):
    assert set(_gene_query_records(_setup(tmp_path))) == {"MAT1-1-1", "APN2", "SLA2"}


def test_exonerate_is_called_with_the_gene_subset(tmp_path):
    ref, seen = _setup(tmp_path), {}

    def runner(cmd, **kwargs):
        q = cmd[cmd.index("--query") + 1]
        seen["query"] = open(q).read()
        return _Result(0, _GFF)

    model = polish_with_exonerate(tmp_path / "genome.fa", FAMILY, "MAT1-1-1", ref, {"recA": FAMILY.key, "recB": FAMILY.key},
                                  ("c1", 1, 600), runner=runner)
    assert seen["query"] == ">recA|gene0|MAT1-1-1 first\nMKLV\n>recB|gene0|MAT1-1-1\nMQQQ\n"
    assert model is not None and model.reference_record_id == "recB"


def test_exonerate_is_not_run_for_a_gene_with_no_reference_protein(tmp_path):
    ref, calls = _setup(tmp_path), []
    assert polish_with_exonerate(tmp_path / "genome.fa", FAMILY, "NOPE", ref, {}, ("c1", 1, 600),
                                 runner=lambda cmd, **k: calls.append(cmd) or _Result(0, "")) is None
    assert calls == []


def test_a_changed_reference_file_is_re_read(tmp_path):
    ref = _setup(tmp_path)
    assert "MKLV" in _write_gene_query(ref, "MAT1-1-1", tmp_path).read_text()
    ref.write_text(">recZ|gene0|MAT1-1-1\nMNEW\n")
    assert _write_gene_query(ref, "MAT1-1-1", tmp_path).read_text() == ">recZ|gene0|MAT1-1-1\nMNEW\n"
