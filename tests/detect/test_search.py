from __future__ import annotations
from pathlib import Path

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import SearchHit, search_fast_path, search_genomic

FAMILY = Family(
    key=FamilyKey("Basidiomycota", "aLocus"),
    vocabulary_type="pattern",
    idiomorph_values=None,
    idiomorph_pattern="^a[0-9]+$",
    genes=[{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}],
    taxonomic_scope=[5270],
)

# Real diamond blastp (protein-vs-protein) has no native genomic-coordinate
# columns. We rely on the predicted proteome's own FASTA deflines encoding
# `contig:start-end:strand` after the first whitespace-delimited token, and
# request diamond's `qtitle` field (everything after qseqid on the query
# defline) to recover it. --outfmt is therefore:
#   6 qseqid sseqid pident qtitle
DIAMOND_TSV = "query1\t5270_521_aLocus_a1|gene0|mfa1\t95.0\tcontigA:100-400:+\n"


def fake_diamond_runner(cmd, **kwargs):
    class Result:
        returncode = 0
        stdout = DIAMOND_TSV
        stderr = ""

    return Result()


def test_search_fast_path_parses_diamond_output(tmp_path):
    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        runner=fake_diamond_runner,
    )
    assert hits == [
        SearchHit(
            family_key=FamilyKey("Basidiomycota", "aLocus"),
            gene_name="mfa1",
            role="core_MAT",
            contig="contigA",
            start=100,
            end=400,
            strand="+",
            identity=95.0,
            reference_record_id="5270_521_aLocus_a1",
            method="diamond_proteome",
        )
    ]


def test_search_fast_path_ignores_hits_for_unknown_genes(tmp_path):
    tsv = "query1\t5270_521_aLocus_a1|gene9|unknown_gene\t80.0\tcontigA:1-2:+\n"

    def fake_runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stdout = tsv
            stderr = ""

        return Result()

    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        runner=fake_runner,
    )
    assert hits == []


EXONERATE_GFF = (
    "contigB\texonerate\tgene\t500\t900\t.\t-\t.\t"
    "gene_id 1 ; sequence 5270_521_aLocus_a1|gene1|pra1 ; gene_orientation -\n"
)


def fake_exonerate_runner(cmd, **kwargs):
    class Result:
        returncode = 0
        stdout = EXONERATE_GFF
        stderr = ""

    return Result()


def test_search_genomic_parses_exonerate_output_no_window(tmp_path):
    hits = search_genomic(
        genome_fasta=tmp_path / "genome.fasta",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        runner=fake_exonerate_runner,
    )
    assert hits == [
        SearchHit(
            family_key=FamilyKey("Basidiomycota", "aLocus"),
            gene_name="pra1",
            role="core_MAT",
            contig="contigB",
            start=500,
            end=900,
            strand="-",
            identity=0.0,
            reference_record_id="5270_521_aLocus_a1",
            method="exonerate_genome",
        )
    ]


def test_search_genomic_relaxed_sets_relaxed_method(tmp_path):
    hits = search_genomic(
        genome_fasta=tmp_path / "genome.fasta",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        relaxed=True,
        runner=fake_exonerate_runner,
    )
    assert hits[0].method == "exonerate_genome_relaxed"


def test_search_genomic_window_extracts_target_region_and_offsets_coordinates(tmp_path):
    # 1000bp contig; the "gene" of interest occupies genomic 1500..1900 within
    # a contig spanning 1..2000 -- write a real contig long enough that
    # slicing out a 1000..2000 window and re-numbering from 1 is meaningful.
    genome_fasta = tmp_path / "genome.fasta"
    genome_fasta.write_text(">contigC\n" + ("A" * 2000) + "\n")

    captured_cmd = {}

    def fake_runner(cmd, **kwargs):
        captured_cmd["cmd"] = cmd
        # the target passed to exonerate must be a real (sliced) FASTA file,
        # not the full genome -- assert this here, while the temp file the
        # implementation created still exists (it is cleaned up once
        # search_genomic returns).
        target_arg = cmd[cmd.index("--target") + 1]
        assert target_arg != str(genome_fasta)
        sliced = Path(target_arg).read_text()
        # window (1000, 2000) is 1-based inclusive -> 1001 bases
        assert sliced.count("A") == 1001
        # exonerate ran against the *sliced* window fasta, so its GFF
        # coordinates are local to that slice (1-based from the window start):
        # a hit at genomic 1500-1900 in a window starting at genomic 1000
        # is reported locally as 501-901.
        class Result:
            returncode = 0
            stdout = (
                "contigC\texonerate\tgene\t501\t901\t.\t+\t.\t"
                "gene_id 1 ; sequence 5270_521_aLocus_a1|gene0|mfa1 ; gene_orientation +\n"
            )
            stderr = ""

        return Result()

    hits = search_genomic(
        genome_fasta=genome_fasta,
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        window=("contigC", 1000, 2000),
        runner=fake_runner,
    )

    assert hits == [
        SearchHit(
            family_key=FamilyKey("Basidiomycota", "aLocus"),
            gene_name="mfa1",
            role="core_MAT",
            contig="contigC",
            start=1500,
            end=1900,
            strand="+",
            identity=0.0,
            reference_record_id="5270_521_aLocus_a1",
            method="exonerate_genome",
        )
    ]
    assert captured_cmd["cmd"]  # the fake runner's own assertions above ran
