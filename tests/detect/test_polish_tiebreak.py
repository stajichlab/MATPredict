"""When two curated reference records tie exactly for one gene, the record named in the model must not depend on the order in
which the polishing tool happens to list the alignments (seen on Leppa1: the same model attributed to the A1163 record in most
runs and to the Af293 record in some). Seam: polish_with_exonerate with a stubbed tool output.
"""
from __future__ import annotations

from MATPredict.detect.search import polish_with_exonerate, polish_with_miniprot

from tests.detect.test_search import FAMILY


def _gff(rows):
    """rows: (record id, score[, identity]). One alignment per row for gene mfa1, identical coordinates; identity 59.60 unless given."""
    out = []
    for row in rows:
        rec, score = row[0], row[1]
        ident = row[2] if len(row) > 2 else 59.60
        out.append(f"c1\texonerate\tgene\t1\t400\t{score}\t+\t.\tgene_id 1 ; sequence {rec}|gene0|mfa1 ; gene_orientation . ; identity {ident:.2f} ; similarity 70.00\n")
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
    """Equal score AND equal identity: the lower record id decides."""
    one, two = tmp_path / "one", tmp_path / "two"
    one.mkdir(); two.mkdir()
    first = _record_chosen(one, [("rec_a1163", 500), ("rec_af293", 500)])
    second = _record_chosen(two, [("rec_af293", 500), ("rec_a1163", 500)])
    assert first == second == "rec_a1163"          # the lower record id, by rule


def test_a_higher_score_still_beats_a_lower_record_id(tmp_path):
    assert _record_chosen(tmp_path, [("rec_a1163", 500), ("rec_af293", 501)]) == "rec_af293"


def test_equal_score_goes_to_the_higher_identity_whatever_the_record_ids_or_order(tmp_path):
    """Equal score, different identity: the higher identity wins over the lower record id, in both listing orders (exonerate)."""
    one, two = tmp_path / "one", tmp_path / "two"
    one.mkdir(); two.mkdir()
    a = _record_chosen(one, [("rec_a1163", 500, 27.06), ("rec_af293", 500, 31.71)])
    b = _record_chosen(two, [("rec_af293", 500, 31.71), ("rec_a1163", 500, 27.06)])
    assert a == b == "rec_af293"                    # not the lower id: identity decides first


def test_a_higher_score_still_beats_a_higher_identity(tmp_path):
    assert _record_chosen(tmp_path, [("rec_a1163", 501, 20.0), ("rec_af293", 500, 90.0)]) == "rec_a1163"


def _miniprot_record_chosen(tmp_path, rows):
    """rows: (record id, score, identity as a fraction 0-1), one mRNA per row for gene mfa1."""
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    out = "##gff-version 3\n"
    for i, (rec, score, ident) in enumerate(rows, 1):
        out += (f"c1\tminiprot\tmRNA\t1\t400\t{score}\t+\t.\tID=MP{i:06d};Rank=1;Identity={ident:.4f};Positive=0.9;Target={rec}|gene0|mfa1 1 76\n"
                f"c1\tminiprot\tCDS\t1\t150\t144\t+\t0\tParent=MP{i:06d};Rank=1;Identity={ident:.4f};Target={rec}|gene0|mfa1 1 31\n"
                f"c1\tminiprot\tCDS\t200\t400\t257\t+\t0\tParent=MP{i:06d};Rank=1;Identity={ident:.4f};Target={rec}|gene0|mfa1 32 76\n")

    def runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stdout = out
            stderr = ""
        return Result()

    return polish_with_miniprot(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1", reference_fasta=tmp_path / "reference.faa",
        record_families={"rec_a1163": FAMILY.key, "rec_af293": FAMILY.key}, window=("c1", 1, 500), runner=runner).reference_record_id


def test_miniprot_equal_score_goes_to_the_higher_identity_then_the_lower_id(tmp_path):
    one, two, three = tmp_path / "one", tmp_path / "two", tmp_path / "three"
    for d in (one, two, three):
        d.mkdir()
    assert _miniprot_record_chosen(one, [("rec_a1163", 400, 0.50), ("rec_af293", 400, 0.80)]) == "rec_af293"
    assert _miniprot_record_chosen(two, [("rec_af293", 400, 0.80), ("rec_a1163", 400, 0.50)]) == "rec_af293"
    assert _miniprot_record_chosen(three, [("rec_af293", 400, 0.70), ("rec_a1163", 400, 0.70)]) == "rec_a1163"
