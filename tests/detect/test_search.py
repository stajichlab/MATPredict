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
    polish_with_miniprot,
    search_fast_path,
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


# --- polish_with_miniprot (exon-aware miniprot --gff polishing) ---

# Verified against the real miniprot 0.18-r281 binary, invoked as
# `miniprot --gff <target.fa> <query.faa>` (target genome FIRST, query
# protein SECOND -- the opposite argument order from exonerate's
# --target/--query flags) against a synthetic two-exon gene on a 931bp
# contig with a real 300bp GT...AG intron. Real captured stdout (miniprot's
# own progress/version banner goes to stderr, not stdout):
#
#   ##gff-version 3
#   ##PAF	rec1|gene0|mfa1	76	0	76	+	chr1	931	200	728	228	228	0	AS:i:364	ms:i:401	np:i:76	fs:i:0	st:i:0	da:i:0	do:i:0	cg:Z:31M300N45M	cs:Z::31~gt300ag:45
#   chr1	miniprot	mRNA	201	731	401	+	.	ID=MP000001;Rank=1;Identity=1.0000;Positive=1.0000;Target=rec1|gene0|mfa1 1 76
#   chr1	miniprot	CDS	201	293	144	+	0	Parent=MP000001;Rank=1;Identity=1.0000;Target=rec1|gene0|mfa1 1 31
#   chr1	miniprot	CDS	594	731	257	+	0	Parent=MP000001;Rank=1;Identity=1.0000;Target=rec1|gene0|mfa1 32 76
#   chr1	miniprot	stop_codon	729	731	0	+	0	Parent=MP000001;Rank=1
#
# So: a `##PAF` comment line to skip, a GFF3 `mRNA` feature line carrying
# `ID=`, `Target=<query_id> <qstart> <qend>` (space-separated, not another
# `key=value` pair) and `Identity=<fraction 0-1>` (NOT a percentage, unlike
# exonerate's `identity <pct>`) among its `;`-separated attributes, and one
# `CDS` line per exon carrying `Parent=<mRNA ID>` linking it back. A
# `stop_codon` line is also emitted and ignored here. A query with no hit in
# the window produces only the `##gff-version 3` header, no mRNA/CDS lines,
# with exit code 0 -- confirmed by running miniprot with an unrelated
# 49-residue protein against the same genome.
MINIPROT_GFF = (
    "##gff-version 3\n"
    "##PAF\trec1|gene0|mfa1\t76\t0\t76\t+\tc1\t500\t0\t400\t228\t228\t0\tAS:i:364\n"
    "c1\tminiprot\tmRNA\t1\t400\t401\t+\t.\t"
    "ID=MP000001;Rank=1;Identity=0.9500;Positive=0.9600;Target=rec1|gene0|mfa1 1 76\n"
    "c1\tminiprot\tCDS\t1\t150\t144\t+\t0\t"
    "Parent=MP000001;Rank=1;Identity=0.9500;Target=rec1|gene0|mfa1 1 31\n"
    "c1\tminiprot\tCDS\t200\t400\t257\t+\t0\t"
    "Parent=MP000001;Rank=1;Identity=0.9500;Target=rec1|gene0|mfa1 32 76\n"
    "c1\tminiprot\tstop_codon\t398\t400\t0\t+\t0\tParent=MP000001;Rank=1\n"
)


def fake_miniprot_runner(cmd, **kwargs):
    assert cmd[0] == "miniprot" and "--gff" in cmd
    # target (window fasta) comes before the query/reference fasta.
    gff_flag_index = cmd.index("--gff")
    assert cmd[gff_flag_index + 2].endswith("reference.faa")

    class Result:
        returncode = 0
        stdout = MINIPROT_GFF
        stderr = ""
    return Result()


def test_polish_with_miniprot_parses_exon_structure(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")
    model = polish_with_miniprot(
        genome_fasta=tmp_path / "genome.fa",
        family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa",
        record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500),
        runner=fake_miniprot_runner,
    )
    assert model.exons == [ExonSpan(1, 150), ExonSpan(200, 400)]
    assert model.identity == 95.0
    assert model.method == "miniprot_refine"
    assert model.contig == "c1"
    assert (model.start, model.end, model.strand) == (1, 400, "+")


def test_polish_with_miniprot_returns_none_when_no_model(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")

    def empty_runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stdout = "##gff-version 3\n"
            stderr = ""
        return Result()

    model = polish_with_miniprot(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=empty_runner,
    )
    assert model is None


def test_polish_with_miniprot_returns_none_when_gene_name_mismatches(tmp_path):
    """Regression: when a padded window overlaps two adjacent genes, miniprot
    may return an alignment for the neighboring gene instead of the requested
    one (e.g., request mfa1 but miniprot finds pra1 in the same window). The
    guard against this cross-gene misattribution must reject the result."""
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 500 + "\n")

    def mismatched_gene_runner(cmd, **kwargs):
        # GFF reports pra1 (different gene, same record)
        gff = (
            "##gff-version 3\n"
            "c1\tminiprot\tmRNA\t1\t400\t401\t+\t.\t"
            "ID=MP000001;Rank=1;Identity=0.9500;Positive=0.9600;Target=rec1|gene1|pra1 1 76\n"
            "c1\tminiprot\tCDS\t1\t150\t144\t+\t0\t"
            "Parent=MP000001;Rank=1;Identity=0.9500;Target=rec1|gene1|pra1 1 31\n"
            "c1\tminiprot\tCDS\t200\t400\t257\t+\t0\t"
            "Parent=MP000001;Rank=1;Identity=0.9500;Target=rec1|gene1|pra1 32 76\n"
        )
        class Result:
            returncode = 0
            stdout = gff
            stderr = ""
        return Result()

    # Request mfa1, but the fake miniprot output carries pra1
    model = polish_with_miniprot(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 500), runner=mismatched_gene_runner,
    )
    assert model is None


# --- Multi-gene windows: the normal shape of a real MAT locus ---
#
# Both wrappers are handed the WHOLE reference FASTA (every curated protein of
# every routed family) as query, and a padded window routinely contains more
# than one gene. Captured from the real binaries against the real curated
# Basidiomycota:Aalpha Z and Y proteins placed in one contig:
#
#   exonerate 2.4.0 emits one `gene` line per alignment, BEST-SCORING FIRST,
#   each followed by its own `exon` lines, with `gene_id` restarting at 1 for
#   every alignment (so gene_id cannot group them -- position does).
#   miniprot 0.18-r281 emits one `mRNA` per aligning query in QUERY ORDER,
#   each with its own Parent-linked `CDS` lines.
#
# The two fixtures below mirror that, renamed onto this file's FAMILY
# (mfa1/pra1): the NON-requested gene's alignment comes first in both.

EXONERATE_TWO_GENE_GFF = (
    # pra1 first (exonerate orders by score), with its own exons
    "c1\texonerate\tgene\t2001\t2400\t4745\t+\t.\t"
    "gene_id 1 ; sequence rec1|gene1|pra1 ; gene_orientation . ; identity 88.00\n"
    "c1\texonerate\texon\t2001\t2100\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    "c1\texonerate\texon\t2200\t2400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    # mfa1 second, with its own, entirely separate exons
    "c1\texonerate\tgene\t1\t400\t4767\t+\t.\t"
    "gene_id 1 ; sequence rec1|gene0|mfa1 ; gene_orientation . ; identity 95.00\n"
    "c1\texonerate\texon\t1\t150\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    "c1\texonerate\texon\t200\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
)


def _stdout_runner(stdout):
    def runner(cmd, **kwargs):
        class Result:
            returncode = 0
            stderr = ""
        Result.stdout = stdout
        return Result()

    return runner


def test_polish_with_exonerate_keeps_only_the_requested_genes_exons_in_a_two_gene_window(tmp_path):
    """A window holding two genes must not mix their exons, and must not answer
    None for the second gene just because the first one scored higher."""
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 3000 + "\n")
    runner = _stdout_runner(EXONERATE_TWO_GENE_GFF)
    common = dict(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY,
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 3000), runner=runner,
    )

    mfa1 = polish_with_exonerate(gene_name="mfa1", **common)
    assert mfa1 is not None  # not None just because pra1's alignment came first
    assert mfa1.gene_name == "mfa1"
    assert (mfa1.start, mfa1.end) == (1, 400)
    # pra1's exon spans (2001-2100, 2200-2400) must NOT appear here
    assert mfa1.exons == [ExonSpan(1, 150), ExonSpan(200, 400)]
    assert mfa1.identity == 95.0

    pra1 = polish_with_exonerate(gene_name="pra1", **common)
    assert pra1 is not None
    assert pra1.gene_name == "pra1"
    assert (pra1.start, pra1.end) == (2001, 2400)
    assert pra1.exons == [ExonSpan(2001, 2100), ExonSpan(2200, 2400)]
    assert pra1.identity == 88.0


MINIPROT_TWO_GENE_GFF = (
    "##gff-version 3\n"
    "##PAF\trec1|gene1|pra1\t76\t0\t76\t+\tc1\t3000\t2000\t2400\t228\t228\t0\tAS:i:4745\n"
    # pra1 first (miniprot orders by query), with its own Parent-linked CDS lines
    "c1\tminiprot\tmRNA\t2001\t2400\t4745\t+\t.\t"
    "ID=MP000001;Rank=1;Identity=0.8800;Target=rec1|gene1|pra1 1 76\n"
    "c1\tminiprot\tCDS\t2001\t2100\t144\t+\t0\tParent=MP000001;Rank=1;Target=rec1|gene1|pra1 1 31\n"
    "c1\tminiprot\tCDS\t2200\t2400\t257\t+\t0\tParent=MP000001;Rank=1;Target=rec1|gene1|pra1 32 76\n"
    "c1\tminiprot\tstop_codon\t2398\t2400\t0\t+\t0\tParent=MP000001;Rank=1\n"
    "##PAF\trec1|gene0|mfa1\t76\t0\t76\t+\tc1\t3000\t0\t400\t228\t228\t0\tAS:i:4767\n"
    # mfa1 second, with its own, entirely separate CDS lines
    "c1\tminiprot\tmRNA\t1\t400\t4767\t+\t.\t"
    "ID=MP000002;Rank=1;Identity=0.9500;Target=rec1|gene0|mfa1 1 76\n"
    "c1\tminiprot\tCDS\t1\t150\t144\t+\t0\tParent=MP000002;Rank=1;Target=rec1|gene0|mfa1 1 31\n"
    "c1\tminiprot\tCDS\t200\t400\t257\t+\t0\tParent=MP000002;Rank=1;Target=rec1|gene0|mfa1 32 76\n"
)


def test_polish_with_miniprot_keeps_only_the_requested_genes_exons_in_a_two_gene_window(tmp_path):
    """miniprot answers for every aligning query, in query order. Asking for the
    second gene must return that gene's own model, not None."""
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 3000 + "\n")
    runner = _stdout_runner(MINIPROT_TWO_GENE_GFF)
    common = dict(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY,
        reference_fasta=tmp_path / "reference.faa", record_families={"rec1": FAMILY.key},
        window=("c1", 1, 3000), runner=runner,
    )

    mfa1 = polish_with_miniprot(gene_name="mfa1", **common)
    assert mfa1 is not None  # not None just because pra1's mRNA came first
    assert mfa1.gene_name == "mfa1"
    assert (mfa1.start, mfa1.end) == (1, 400)
    assert mfa1.exons == [ExonSpan(1, 150), ExonSpan(200, 400)]
    assert mfa1.identity == 95.0

    pra1 = polish_with_miniprot(gene_name="pra1", **common)
    assert pra1 is not None
    assert pra1.gene_name == "pra1"
    assert (pra1.start, pra1.end) == (2001, 2400)
    assert pra1.exons == [ExonSpan(2001, 2100), ExonSpan(2200, 2400)]
    assert pra1.identity == 88.0


# --- Spec Stage 2: select per gene per tool by the TOOL'S OWN score ---
#
# The curated database legitimately holds several reference proteins for the
# same gene (different curated records of the same family), each producing its
# own alignment in the same window. The spec requires selecting one of them by
# a named, tool-appropriate score -- explicitly NOT raw percent identity. Both
# fixtures below make the two criteria disagree: the LOWER-identity alignment
# carries the HIGHER tool score, so a score-based selection and an
# identity-based one pick different records.

TWO_RECORD_FAMILIES = {"rec1": FAMILY.key, "rec2": FAMILY.key}


def test_polish_with_exonerate_selects_the_best_by_exonerate_score_not_identity(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 3000 + "\n")
    gff = (
        "c1\texonerate\tgene\t1\t400\t1200\t+\t.\t"
        "gene_id 1 ; sequence rec1|gene0|mfa1 ; gene_orientation . ; identity 99.00\n"
        "c1\texonerate\texon\t1\t400\t.\t+\t.\tinsertions 0 ; deletions 0\n"
        "c1\texonerate\tgene\t1000\t2600\t4767\t+\t.\t"
        "gene_id 1 ; sequence rec2|gene0|mfa1 ; gene_orientation . ; identity 71.00\n"
        "c1\texonerate\texon\t1000\t2600\t.\t+\t.\tinsertions 0 ; deletions 0\n"
    )
    model = polish_with_exonerate(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families=TWO_RECORD_FAMILIES,
        window=("c1", 1, 3000), runner=_stdout_runner(gff),
    )
    # rec2 scores 4767 vs rec1's 1200, though rec1 has the higher identity.
    assert model.reference_record_id == "rec2"
    assert (model.start, model.end) == (1000, 2600)
    assert model.exons == [ExonSpan(1000, 2600)]


def test_polish_with_miniprot_selects_the_best_by_miniprot_score_not_identity(tmp_path):
    (tmp_path / "genome.fa").write_text(">c1\n" + "N" * 3000 + "\n")
    gff = (
        "##gff-version 3\n"
        "c1\tminiprot\tmRNA\t1\t400\t1200\t+\t.\t"
        "ID=MP000001;Rank=1;Identity=0.9900;Target=rec1|gene0|mfa1 1 76\n"
        "c1\tminiprot\tCDS\t1\t400\t1200\t+\t0\tParent=MP000001;Rank=1\n"
        "c1\tminiprot\tmRNA\t1000\t2600\t4767\t+\t.\t"
        "ID=MP000002;Rank=1;Identity=0.7100;Target=rec2|gene0|mfa1 1 76\n"
        "c1\tminiprot\tCDS\t1000\t2600\t4767\t+\t0\tParent=MP000002;Rank=1\n"
    )
    model = polish_with_miniprot(
        genome_fasta=tmp_path / "genome.fa", family=FAMILY, gene_name="mfa1",
        reference_fasta=tmp_path / "reference.faa", record_families=TWO_RECORD_FAMILIES,
        window=("c1", 1, 3000), runner=_stdout_runner(gff),
    )
    assert model.reference_record_id == "rec2"
    assert (model.start, model.end) == (1000, 2600)
    assert model.exons == [ExonSpan(1000, 2600)]


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
