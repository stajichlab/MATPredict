"""When two curated reference records tie exactly for one gene, the record named in the model must not depend on the order in
which the polishing tool happens to list the alignments (seen on Leppa1: the same model attributed to the A1163 record in most
runs and to the Af293 record in some). Seam: polish_with_exonerate with a stubbed tool output.
"""
from __future__ import annotations

from MATPredict.detect.search import polish_with_exonerate

from tests.detect.test_search import FAMILY


def _gff(rows):
    """rows: (record id, score). One alignment per row for gene mfa1, identical coordinates and identity."""
    out = []
    for rec, score in rows:
        out.append(f"c1\texonerate\tgene\t1\t400\t{score}\t+\t.\tgene_id 1 ; sequence {rec}|gene0|mfa1 ; gene_orientation . ; identity 59.60 ; similarity 70.00\n")
        out.append("c1\texonerate\texon\t1\t150\t.\t+\t.\tinsertions 0 ; deletions 0\n")
        out.append("c1\texonerate\texon\t200\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n")
    return "".join(out)


def _record_chosen(tmp_path, rows):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")

    def runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stdout = _gff(rows)
            stderr = ""
        return Result()

    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1", reference_fasta=tmp_path / "reference.faa",
        record_families={"rec_a1163": FAMILY.key, "rec_af293": FAMILY.key}, window=("c1", 1, 500), runner=runner)
    return model.reference_record_id


def test_an_exact_tie_names_the_same_record_whatever_order_the_tool_lists_them(tmp_path):
    one, two = tmp_path / "one", tmp_path / "two"
    one.mkdir(); two.mkdir()
    first = _record_chosen(one, [("rec_a1163", 500), ("rec_af293", 500)])
    second = _record_chosen(two, [("rec_af293", 500), ("rec_a1163", 500)])
    assert first == second == "rec_a1163"          # the lower record id, by rule


def test_a_higher_score_still_beats_a_lower_record_id(tmp_path):
    assert _record_chosen(tmp_path, [("rec_a1163", 500), ("rec_af293", 501)]) == "rec_af293"
