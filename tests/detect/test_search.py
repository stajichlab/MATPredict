from __future__ import annotations
from pathlib import Path

import pytest

from MATPredict.detect.family_registry import Family, FamilyKey, load_record_families
from MATPredict.detect.polish import ExonSpan
from MATPredict.detect.search import (
    METHOD_TBLASTN,
    ProteomeDeflineError,
    SearchHit,
    SearchToolError,
    polish_with_exonerate,
    search_fast_path,
    search_genomic,
    search_localize,
)

FAMILY = Family(
    key=FamilyKey("Basidiomycota", "aLocus"),
    vocabulary_type="pattern",
    idiomorph_values=None,
    idiomorph_pattern="^a[0-9]+$",
    genes=[{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}],
    taxonomic_scope=[5270],
)

RECORD_FAMILIES = {"5270_521_aLocus_a1": FAMILY.key}

# Real diamond blastp (protein-vs-protein) has no native genomic-coordinate
# columns. We rely on the predicted proteome's own FASTA deflines encoding
# `contig:start-end:strand` after the first whitespace-delimited token, and
# request diamond's `qtitle` field (everything after qseqid on the query
# defline) to recover it. `scovhsp` gives coverage of the curated reference
# protein. --outfmt is therefore:
#   6 qseqid sseqid pident scovhsp qtitle
# Verified against diamond v2.2.6: qtitle is the WHOLE defline, query id included.
DIAMOND_TSV = "query1\t5270_521_aLocus_a1|gene0|mfa1\t95.0\t88.0\tquery1 contigA:100-400:+\n"


def _result(stdout: str, returncode: int = 0, stderr: str = ""):
    class Result:
        pass

    Result.returncode = returncode
    Result.stdout = stdout
    Result.stderr = stderr
    return Result()


def fake_diamond_runner(cmd, **kwargs):
    # `diamond makedb` produces no tabular output; only `blastp` does.
    if "makedb" in cmd:
        return _result("")
    return _result(DIAMOND_TSV)


def test_search_fast_path_parses_diamond_output(tmp_path):
    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
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
            coverage=88.0,
        )
    ]


def test_search_fast_path_builds_a_diamond_database_first(tmp_path):
    """`diamond blastp --db` needs a binary .dmnd database, not a plain FASTA."""
    commands = []

    def capturing_runner(cmd, **kwargs):
        commands.append(cmd)
        return fake_diamond_runner(cmd, **kwargs)

    search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
        runner=capturing_runner,
    )
    assert commands[0][:2] == ["diamond", "makedb"]
    assert commands[1][:2] == ["diamond", "blastp"]
    blastp_db = commands[1][commands[1].index("--db") + 1]
    assert blastp_db.endswith(".dmnd")


def test_search_fast_path_ignores_hits_for_unknown_genes(tmp_path):
    tsv = "query1\t5270_521_aLocus_a1|gene9|unknown_gene\t80.0\t50.0\tquery1 contigA:1-2:+\n"

    def fake_runner(cmd, **kwargs):
        return _result("" if "makedb" in cmd else tsv)

    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
        runner=fake_runner,
    )
    assert hits == []


def test_search_fast_path_raises_on_nonzero_returncode(tmp_path):
    """A failed diamond invocation must not be indistinguishable from zero hits."""

    def failing_runner(cmd, **kwargs):
        return _result("", returncode=127, stderr="diamond: command not found")

    with pytest.raises(SearchToolError) as err:
        search_fast_path(
            proteome_fasta=tmp_path / "proteome.faa",
            families=[FAMILY],
            reference_fasta=tmp_path / "reference.faa",
            record_families=RECORD_FAMILIES,
            runner=failing_runner,
        )
    assert "diamond" in str(err.value)


@pytest.mark.parametrize(
    "defline",
    [
        "lcl|ABC123.1_prot_XP_001.1_1 [gene=mfa1] [protein=x]",  # NCBI style
        "g1.t1",  # AUGUSTUS/BRAKER style
        "contigA:100-400",  # missing strand
        "contigA:abc-400:+",  # non-numeric coordinates
    ],
)
def test_search_fast_path_raises_clear_error_on_non_conforming_defline(tmp_path, defline):
    tsv = f"query1\t5270_521_aLocus_a1|gene0|mfa1\t95.0\t88.0\t{defline}\n"

    def fake_runner(cmd, **kwargs):
        return _result("" if "makedb" in cmd else tsv)

    with pytest.raises(ProteomeDeflineError) as err:
        search_fast_path(
            proteome_fasta=tmp_path / "proteome.faa",
            families=[FAMILY],
            reference_fasta=tmp_path / "reference.faa",
            record_families=RECORD_FAMILIES,
            runner=fake_runner,
        )
    assert "contig:start-end:strand" in str(err.value)
    assert defline in str(err.value)


EXONERATE_GFF = (
    "contigB\texonerate\tgene\t500\t900\t.\t-\t.\t"
    "gene_id 1 ; sequence 5270_521_aLocus_a1|gene1|pra1 ; gene_orientation - ; "
    "identity 95.50 ; similarity 96.00\n"
)


def fake_exonerate_runner(cmd, **kwargs):
    return _result(EXONERATE_GFF)


def test_search_genomic_parses_exonerate_output_no_window(tmp_path):
    hits = search_genomic(
        genome_fasta=tmp_path / "genome.fasta",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
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
            identity=95.5,
            reference_record_id="5270_521_aLocus_a1",
            method="exonerate_genome",
            coverage=None,
        )
    ]


def test_search_genomic_parses_real_identity_from_exonerate_attrs(tmp_path):
    """Regression: search_genomic previously hardcoded identity=0.0, claiming
    exonerate's GFF has no identity field. The installed exonerate 2.4.0
    binary (--model protein2genome --showtargetgff yes) actually emits a
    real `identity` field in the gene feature's attribute string, in exactly
    this format -- confirmed against the real binary's output."""

    def real_format_runner(cmd, **kwargs):
        return _result(
            "rec1\texonerate\tgene\t10\t400\t.\t+\t.\t"
            "gene_id 1 ; sequence 5270_521_aLocus_a1|gene0|mfa1 ; "
            "gene_orientation . ; identity 100.00 ; similarity 100.00\n"
        )

    hits = search_genomic(
        genome_fasta=tmp_path / "genome.fasta",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
        runner=real_format_runner,
    )
    assert len(hits) == 1
    assert hits[0].identity == 100.0


def test_search_genomic_falls_back_to_zero_identity_when_field_missing(tmp_path):
    """Older/malformed exonerate output without an `identity` field must not
    crash -- it should fall back to 0.0 rather than raising."""

    def no_identity_runner(cmd, **kwargs):
        return _result(
            "contigB\texonerate\tgene\t500\t900\t.\t-\t.\t"
            "gene_id 1 ; sequence 5270_521_aLocus_a1|gene1|pra1 ; gene_orientation -\n"
        )

    hits = search_genomic(
        genome_fasta=tmp_path / "genome.fasta",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
        runner=no_identity_runner,
    )
    assert len(hits) == 1
    assert hits[0].identity == 0.0


def test_search_genomic_relaxed_sets_relaxed_method(tmp_path):
    hits = search_genomic(
        genome_fasta=tmp_path / "genome.fasta",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
        relaxed=True,
        runner=fake_exonerate_runner,
    )
    assert hits[0].method == "exonerate_genome_relaxed"


def test_relaxed_pass_is_more_permissive_than_the_standard_pass(tmp_path):
    """Regression: the relaxed second pass previously added `--percent 50`,
    which DISCARDS matches below 50% of the query's maximal score, while the
    standard pass passed no threshold flag at all -- making the "relaxed" pass
    strictly more stringent. Relaxation is now a LOWER absolute --score."""
    captured: list[list[str]] = []

    def capturing_runner(cmd, **kwargs):
        captured.append(cmd)
        return _result(EXONERATE_GFF)

    for relaxed in (False, True):
        search_genomic(
            genome_fasta=tmp_path / "genome.fasta",
            families=[FAMILY],
            reference_fasta=tmp_path / "reference.faa",
            record_families=RECORD_FAMILIES,
            relaxed=relaxed,
            runner=capturing_runner,
        )

    standard_cmd, relaxed_cmd = captured
    standard_score = int(standard_cmd[standard_cmd.index("--score") + 1])
    relaxed_score = int(relaxed_cmd[relaxed_cmd.index("--score") + 1])
    assert relaxed_score < standard_score
    assert "--percent" not in standard_cmd
    assert "--percent" not in relaxed_cmd


def test_search_genomic_raises_on_nonzero_returncode(tmp_path):
    def failing_runner(cmd, **kwargs):
        return _result("", returncode=1, stderr="exonerate: bad target")

    with pytest.raises(SearchToolError) as err:
        search_genomic(
            genome_fasta=tmp_path / "genome.fasta",
            families=[FAMILY],
            reference_fasta=tmp_path / "reference.faa",
            record_families=RECORD_FAMILIES,
            runner=failing_runner,
        )
    assert "exonerate" in str(err.value)


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
        return _result(
            "contigC\texonerate\tgene\t501\t901\t.\t+\t.\t"
            "gene_id 1 ; sequence 5270_521_aLocus_a1|gene0|mfa1 ; gene_orientation + ; "
            "identity 100.00 ; similarity 100.00\n"
        )

    hits = search_genomic(
        genome_fasta=genome_fasta,
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
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
            identity=100.0,
            reference_record_id="5270_521_aLocus_a1",
            method="exonerate_genome",
            coverage=None,
        )
    ]
    assert captured_cmd["cmd"]  # the fake runner's own assertions above ran


# --- Finding 1 regression: gene-name collisions across REAL families in db/ ---

DB_ROOT = Path("db")

# Confirmed by reading db/Basidiomycota/order.yml: `pheromone` and
# `pheromone_receptor` are declared by THREE distinct families -- PR, Balpha
# and Bbeta -- and `Y`/`Z` by both Aalpha and Abeta. A bare-gene-name lookup
# keeps only whichever family is written into the dict last, so every other
# family sharing the name becomes permanently unreachable.
PR = FamilyKey("Basidiomycota", "PR")
BALPHA = FamilyKey("Basidiomycota", "Balpha")
BBETA = FamilyKey("Basidiomycota", "Bbeta")


def _real_family(families, key):
    return next(f for f in families if f.key == key)


@pytest.mark.skipif(not DB_ROOT.exists(), reason="requires the real curated db/ tree")
def test_real_colliding_gene_names_are_attributed_to_the_right_family(tmp_path):
    from MATPredict.detect.family_registry import load_all_families

    families = load_all_families(DB_ROOT)
    record_families = load_record_families(DB_ROOT)

    # Real curated records, one per colliding family, confirmed present in db/.
    pr_record = "5346_a43-b43-okayama-7_PR_B43"
    balpha_record = "5334_h4-8_Balpha_3"
    bbeta_record = "5334_h4-8_Bbeta_2"
    for record_id, expected in [
        (pr_record, PR), (balpha_record, BALPHA), (bbeta_record, BBETA),
    ]:
        assert record_families[record_id] == expected

    # All three families declare a gene literally named "pheromone_receptor".
    attempted = [_real_family(families, k) for k in (PR, BALPHA, BBETA)]
    for family in attempted:
        assert "pheromone_receptor" in {g["name"] for g in family.genes}

    tsv = "".join(
        f"q{i}\t{record_id}|gene1|pheromone_receptor\t90.0\t80.0\tq{i} c1:{i * 1000}-{i * 1000 + 500}:+\n"
        for i, record_id in enumerate([pr_record, balpha_record, bbeta_record], start=1)
    )

    def fake_runner(cmd, **kwargs):
        return _result("" if "makedb" in cmd else tsv)

    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=attempted,
        reference_fasta=tmp_path / "reference.faa",
        record_families=record_families,
        runner=fake_runner,
    )

    by_record = {h.reference_record_id: h.family_key for h in hits}
    assert by_record == {pr_record: PR, balpha_record: BALPHA, bbeta_record: BBETA}
    # Every family stays reachable -- none is collapsed into another.
    assert {h.family_key for h in hits} == {PR, BALPHA, BBETA}


@pytest.mark.skipif(not DB_ROOT.exists(), reason="requires the real curated db/ tree")
def test_real_colliding_Z_gene_is_attributed_to_aalpha_not_abeta(tmp_path):
    """`Z` is declared by both Basidiomycota:Aalpha and Basidiomycota:Abeta."""
    from MATPredict.detect.family_registry import load_all_families

    families = load_all_families(DB_ROOT)
    record_families = load_record_families(DB_ROOT)
    aalpha_record = "5334_h4-8_Aalpha_4"
    assert record_families[aalpha_record] == FamilyKey("Basidiomycota", "Aalpha")

    attempted = [
        _real_family(families, FamilyKey("Basidiomycota", "Aalpha")),
        _real_family(families, FamilyKey("Basidiomycota", "Abeta")),
    ]
    tsv = f"q1\t{aalpha_record}|gene0|Z\t99.0\t95.0\tq1 c1:1-500:+\n"

    def fake_runner(cmd, **kwargs):
        return _result("" if "makedb" in cmd else tsv)

    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=attempted,
        reference_fasta=tmp_path / "reference.faa",
        record_families=record_families,
        runner=fake_runner,
    )
    assert [h.family_key for h in hits] == [FamilyKey("Basidiomycota", "Aalpha")]


# --- search_localize (tblastn genome-wide localization) ---

# tblastn's query is the curated reference protein set, its database is the
# genome -- the OPPOSITE role assignment from diamond's fast path (where the
# predicted proteome is the query). So qseqid carries the reference header
# (record_id|geneN|gene_name) and sseqid is the genome's own contig name.
# Verified against the real tblastn 2.17.0 binary: `-outfmt "6 qseqid sseqid
# pident length sstart send sframe"` produces exactly these tab-separated
# columns, and a minus-strand HSP reports sstart > send with sframe -1.
TBLASTN_TSV = (
    # qseqid                sseqid  pident length sstart send sframe
    "rec1|gene0|mfa1\tcontigA\t95.0\t40\t400\t100\t-1\n"  # minus strand: sstart > send
    "rec1|gene1|pra1\tcontigA\t90.0\t300\t3600\t4914\t1\n"  # plus strand
)


def fake_tblastn_runner(cmd, **kwargs):
    if cmd[0] == "makeblastdb":
        return _result("")
    assert "-seg" in cmd and cmd[cmd.index("-seg") + 1] == "no"
    return _result(TBLASTN_TSV)


def test_search_localize_normalizes_minus_strand_and_batches_one_call(tmp_path):
    calls = []

    def counting_runner(cmd, **kwargs):
        calls.append(cmd)
        return fake_tblastn_runner(cmd, **kwargs)

    hits = search_localize(
        genome_fasta=tmp_path / "genome.fa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        runner=counting_runner,
    )
    by_gene = {h.gene_name: h for h in hits}
    assert by_gene["mfa1"].start == 100 and by_gene["mfa1"].end == 400 and by_gene["mfa1"].strand == "-"
    assert by_gene["pra1"].start == 3600 and by_gene["pra1"].end == 4914 and by_gene["pra1"].strand == "+"
    assert all(h.method == METHOD_TBLASTN for h in hits)
    # exactly one tblastn invocation regardless of how many families/genes
    tblastn_calls = [c for c in calls if c[0] == "tblastn"]
    assert len(tblastn_calls) == 1


def test_search_localize_raises_on_nonzero_returncode(tmp_path):
    def failing_runner(cmd, **kwargs):
        if cmd[0] == "makeblastdb":
            return _result("")
        return _result("", returncode=1, stderr="tblastn: bad database")

    with pytest.raises(SearchToolError) as err:
        search_localize(
            genome_fasta=tmp_path / "genome.fa",
            families=[FAMILY],
            reference_fasta=tmp_path / "reference.faa",
            record_families={"rec1": FAMILY.key},
            runner=failing_runner,
        )
    assert "tblastn" in str(err.value)


# --- polish_with_exonerate (exon-aware --refine region polishing) ---

# Verified against the real exonerate 2.4.0 binary (--model protein2genome
# --refine region --showtargetgff yes) with a synthetic two-exon gene: the
# gene-line attribute format (sequence/identity) matches search_genomic's
# existing parsing exactly, and exon lines additionally carry
# "identity ... ; similarity ..." attrs beyond "insertions"/"deletions" --
# irrelevant here since only start/end (fields[3]/fields[4]) are parsed from
# exon lines, not their attribute string.
EXONERATE_REFINE_GFF = (
    "c1\texonerate\tgene\t1\t400\t.\t+\t.\t"
    "gene_id 1 ; sequence rec1|gene0|mfa1 ; gene_orientation . ; identity 95.00 ; similarity 96.00\n"
    "c1\texonerate\texon\t1\t150\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    "c1\texonerate\texon\t200\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
)


def fake_exonerate_refine_runner(cmd, **kwargs):
    assert "--refine" in cmd and cmd[cmd.index("--refine") + 1] == "region"
    class Result:
        returncode = 0
        stdout = EXONERATE_REFINE_GFF
        stderr = ""
    return Result()


def test_polish_with_exonerate_parses_exon_structure(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa",
        family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500),
        runner=fake_exonerate_refine_runner,
    )
    assert model.exons == [ExonSpan(1, 150), ExonSpan(200, 400)]
    assert model.identity == 95.0
    assert model.method == "exonerate_refine"


def test_polish_with_exonerate_returns_none_when_no_model(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")

    def empty_runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stdout = ""
            stderr = ""
        return Result()

    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=empty_runner,
    )
    assert model is None


def test_polish_with_exonerate_returns_none_when_gene_name_mismatches(tmp_path):
    """Regression: when a padded window overlaps two adjacent genes, exonerate
    may return an alignment for the neighboring gene instead of the requested one
    (e.g., request mfa1 but exonerate finds pra1 in the same window). The guard
    against this cross-gene misattribution must reject the result."""
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")

    def mismatched_gene_runner(cmd, **kwargs):
        # GFF reports pra1 (different gene, same record)
        gff = (
            "c1\texonerate\tgene\t1\t400\t.\t+\t.\t"
            "gene_id 1 ; sequence rec1|gene1|pra1 ; gene_orientation . ; identity 95.00 ; similarity 96.00\n"
            "c1\texonerate\texon\t1\t150\t.\t+\t.\tinsertions 0 ; deletions 0\n"
            "c1\texonerate\texon\t200\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
        )
        class Result:
            returncode = 0
            stdout = gff
            stderr = ""
        return Result()

    # Request mfa1, but the fake exonerate output carries pra1
    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=mismatched_gene_runner,
    )
    assert model is None


def test_defline_location_is_found_alongside_a_free_text_description(tmp_path):
    """A real proteome defline often carries a description after the location."""
    tsv = (
        "query1\t5270_521_aLocus_a1|gene0|mfa1\t95.0\t88.0\t"
        "query1 contigA:100-400:+ pheromone precursor mfa1\n"
    )

    def fake_runner(cmd, **kwargs):
        return _result("" if "makedb" in cmd else tsv)

    hits = search_fast_path(
        proteome_fasta=tmp_path / "proteome.faa",
        families=[FAMILY],
        reference_fasta=tmp_path / "reference.faa",
        record_families=RECORD_FAMILIES,
        runner=fake_runner,
    )
    assert (hits[0].contig, hits[0].start, hits[0].end, hits[0].strand) == ("contigA", 100, 400, "+")
