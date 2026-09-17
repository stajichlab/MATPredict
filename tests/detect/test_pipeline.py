from __future__ import annotations
from pathlib import Path

from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.search import SearchHit
from MATPredict.detect.pipeline import run_pipeline

FAMILY = Family(FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
                [{"name": "mfa1", "role": "core_MAT"}, {"name": "pra1", "role": "core_MAT"}], [1])

ORDER_YML = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
)


def _write_order(tmp_path, text=ORDER_YML):
    (tmp_path / "P").mkdir(exist_ok=True)
    (tmp_path / "P" / "order.yml").write_text(text)


def _write_record(tmp_path, record_id="rec1", locus_name="aLocus", proteins: str = ""):
    """A curated record dir, so load_record_families can map record_id -> family."""
    record_dir = tmp_path / "P" / "Fam" / record_id
    record_dir.mkdir(parents=True, exist_ok=True)
    (record_dir / "metadata.yaml").write_text(
        f"record_id: {record_id}\n"
        f"mating_type: {{locus_name: {locus_name}, idiomorphs: [a1]}}\n"
    )
    if proteins:
        (record_dir / "proteins.faa").write_text(proteins)
    return record_dir


def test_run_pipeline_end_to_end_with_stubbed_search(tmp_path):
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
                SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa",
        proteome_fasta=tmp_path / "proteome.faa",
        taxid=None,
        db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path,
        search_genomic=fake_genomic,
    )
    assert len(outcome.results) == 1
    result = outcome.results[0]
    assert result.family_key == FAMILY.key
    assert result.confidence == "high"
    assert result.genes_missing == []
    assert outcome.not_detected == []
    assert outcome.families_attempted == [FAMILY.key]


def test_run_pipeline_triggers_genomic_second_pass_on_missing_core_gene(tmp_path):
    _write_order(tmp_path)
    _write_record(tmp_path)
    genomic_calls = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(genome_fasta, families, reference_fasta, record_families,
                     relaxed=False, window=None, runner=None):
        genomic_calls.append(relaxed)
        return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 60.0, "rec1", "exonerate_genome_relaxed")]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert True in genomic_calls  # relaxed second pass was actually invoked
    assert outcome.results[0].confidence == "medium"  # second-pass-confirmed core gene caps at medium
    assert outcome.results[0].genes_missing == []


def test_second_pass_fires_for_a_family_with_zero_proteome_hits(tmp_path):
    """Finding 3 regression. The whole point of the second pass is the
    mfa1-style blind spot: a supplied proteome annotation that misses a
    family's genes ENTIRELY. Such a family has no cluster to anchor a
    windowed search to, so gating the second pass on an existing foothold
    skipped exactly the case the spec calls unconditional."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    windows_searched = []

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return []  # the provided annotation contains nothing for this family at all

    def fake_genomic(genome_fasta, families, reference_fasta, record_families,
                     relaxed=False, window=None, runner=None):
        windows_searched.append(window)
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 0.0, "rec1", "exonerate_genome_relaxed"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 0.0, "rec1", "exonerate_genome_relaxed"),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    # a whole-genome (window=None) relaxed search actually ran for the family
    assert None in windows_searched
    assert len(outcome.results) == 1
    assert outcome.results[0].genes_found == ["mfa1", "pra1"]
    # rescued only via the relaxed second pass -> capped at medium, not high
    assert outcome.results[0].confidence == "medium"


def test_second_pass_in_one_cluster_does_not_cap_a_different_cluster_of_the_same_family(tmp_path):
    """Same family, two independent spatial clusters (e.g. gene-duplication /
    multi-allele co-occurrence). Cluster on contig c1 needs the relaxed
    genomic second pass to confirm mfa1. Cluster on contig c2 gets both core
    genes from the fast path alone -- it must reach "high", not be capped at
    "medium" just because the c1 cluster (same family) needed a second pass."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            # c1 cluster: only pra1 found by the fast path -> needs second pass.
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            # c2 cluster: both core genes found by the fast path -> no second pass needed.
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c2", 100, 200, "+", 95.0, "rec2", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec2", "diamond_proteome"),
        ]

    def fake_genomic(genome_fasta, families, reference_fasta, record_families,
                     relaxed=False, window=None, runner=None):
        contig = window[0] if window else None
        if contig == "c1":
            return [SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 60.0, "rec1",
                               "exonerate_genome_relaxed")]
        return []

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )

    by_contig = {r.contig: r for r in outcome.results}
    assert len(outcome.results) == 2
    assert by_contig["c1"].confidence == "medium"  # this cluster's own second pass caps it
    assert by_contig["c2"].confidence == "high"  # unrelated cluster, same family, must NOT be capped


def test_short_orf_gene_reported_as_not_searchable_not_missing(tmp_path):
    _write_order(tmp_path)
    _write_record(
        tmp_path,
        proteins=(
            ">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 41 + "\n"
            ">rec1|gene_index=1|name=pra1|role=core_MAT\n" + "M" * 300 + "\n"
        ),
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []  # mfa1 genuinely not found even after the relaxed second pass

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert outcome.results[0].genes_missing == []
    assert outcome.results[0].genes_not_searchable == ["mfa1"]


def test_short_orf_split_discriminates_three_buckets(tmp_path):
    """Three core genes, three distinct fates: mfa1 is short (41 aa) and
    genuinely never found -> genes_not_searchable. pra2 is normal-length and
    genuinely never found -> genes_missing. pra1 is found -> genes_found.
    A test with only one gene per bucket can't tell "correctly bucketed" from
    "coincidentally correct because nothing else is missing" -- this proves
    the split logic actually discriminates across all three cases at once."""
    three_gene_family = Family(
        FamilyKey("P", "aLocus"), "pattern", None, "^a[0-9]+$",
        [{"name": "mfa1", "role": "core_MAT"},
         {"name": "pra1", "role": "core_MAT"},
         {"name": "pra2", "role": "core_MAT"}],
        [1],
    )
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: pra1, role: core_MAT}\n"
        "      - {name: pra2, role: core_MAT}\n",
    )
    _write_record(
        tmp_path,
        proteins=(
            ">rec1|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 41 + "\n"
            ">rec1|gene_index=1|name=pra1|role=core_MAT\n" + "M" * 300 + "\n"
            ">rec1|gene_index=2|name=pra2|role=core_MAT\n" + "M" * 300 + "\n"
        ),
    )

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(three_gene_family.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1",
                           "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []  # neither mfa1 nor pra2 is found even after the relaxed second pass

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
        ambiguity_floor=0.3,  # only 1 of 3 genes is found by design -- lower the floor
        # so this single-family cluster still clears scoring and isn't dropped,
        # without needing a fourth gene just to satisfy an unrelated threshold.
    )
    assert outcome.results[0].genes_found == ["pra1"]
    assert outcome.results[0].genes_missing == ["pra2"]
    assert outcome.results[0].genes_not_searchable == ["mfa1"]


def test_short_orf_scan_is_scoped_per_family_and_uses_the_longest_curated_protein(tmp_path):
    """Finding 2 regression, modelled on the real MATsc `cha1` case: one
    curated record holds a 22-aa FRAGMENT of a gene that is normal length
    elsewhere. Scanning globally and keeping the SHORTEST length made that
    fragment condemn the gene everywhere, falsely reporting a perfectly
    searchable gene as "not searchable by this method"."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n"
        "  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT}\n      - {name: cha1, role: flanking_variable}\n"
        "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: cha1, role: flanking_variable}\n",
    )
    # aLocus: one record has a full-length cha1, another has only a 22-aa fragment.
    _write_record(tmp_path, "recFull", "aLocus", proteins=(
        ">recFull|gene_index=0|name=mfa1|role=core_MAT\n" + "M" * 300 + "\n"
        ">recFull|gene_index=1|name=cha1|role=flanking_variable\n" + "M" * 400 + "\n"
    ))
    _write_record(tmp_path, "recFrag", "aLocus", proteins=(
        ">recFrag|gene_index=0|name=cha1|role=flanking_variable\n" + "M" * 22 + "\n"
    ))
    # bLocus's own cha1 evidence is genuinely short everywhere.
    _write_record(tmp_path, "recB", "bLocus", proteins=(
        ">recB|gene_index=0|name=bE|role=core_MAT\n" + "M" * 300 + "\n"
        ">recB|gene_index=1|name=cha1|role=flanking_variable\n" + "M" * 20 + "\n"
    ))

    a_key, b_key = FamilyKey("P", "aLocus"), FamilyKey("P", "bLocus")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(a_key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "recFull", "diamond_proteome"),
            SearchHit(b_key, "bE", "core_MAT", "c9", 100, 200, "+", 95.0, "recB", "diamond_proteome"),
        ]

    def fake_genomic(*args, **kwargs):
        return []

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    by_family = {r.family_key: r for r in outcome.results}
    # aLocus has a full-length curated cha1 -> genuinely missing, not "unsearchable"
    assert by_family[a_key].genes_missing == ["cha1"]
    assert by_family[a_key].genes_not_searchable == []
    # bLocus's best curated cha1 really is short -> not searchable by this method
    assert by_family[b_key].genes_missing == []
    assert by_family[b_key].genes_not_searchable == ["cha1"]


def test_sub_floor_families_are_reported_as_not_detected_not_dropped(tmp_path):
    """Finding 6 regression: the spec requires a sub-floor result to be
    reported as "not detected" listing which families were attempted and why
    each fell short -- never silently omitted."""
    _write_order(
        tmp_path,
        "phylum: P\nloci:\n"
        "  - locus_name: aLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: g1, role: core_MAT}\n      - {name: g2, role: core_MAT}\n"
        "      - {name: g3, role: core_MAT}\n      - {name: g4, role: core_MAT}\n"
        "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: bW, role: core_MAT}\n",
    )
    _write_record(tmp_path, "recA", "aLocus")
    _write_record(tmp_path, "recB", "bLocus")
    a_key, b_key = FamilyKey("P", "aLocus"), FamilyKey("P", "bLocus")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        # 1 of aLocus's 4 genes -> fraction 0.25, below the 0.5 floor.
        return [SearchHit(a_key, "g1", "core_MAT", "c1", 100, 200, "+", 90.0, "recA", "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert outcome.results == []
    not_detected = {n.family_key: n for n in outcome.not_detected}
    assert set(not_detected) == {a_key, b_key}
    assert not_detected[a_key].best_fraction_found == 0.25
    assert "below the ambiguity floor" in not_detected[a_key].reason
    assert not_detected[a_key].genes_found == ["g1"]
    assert not_detected[b_key].best_fraction_found == 0.0
    assert "no reference-protein hits" in not_detected[b_key].reason
    assert sorted(outcome.families_attempted, key=lambda k: k.locus_name) == [a_key, b_key]


def test_isolated_single_hit_is_low_tier(tmp_path):
    """Finding 6, second half: "low" was unreachable because score_cluster
    never emits a score for a family with zero hits. Per the spec's own Low
    description, a single isolated hit with nothing else from the family
    nearby is Low."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome")]

    def fake_genomic(*args, **kwargs):
        return []

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert outcome.results[0].confidence == "low"


def test_genes_split_across_contigs_are_one_fragmented_multi_segment_call(tmp_path):
    """Finding 5: a family whose core genes land on different contigs, with no
    single cluster carrying them all, is one multi-segment locus with
    fragmented=True and a one-tier confidence downgrade."""
    _write_order(tmp_path)
    _write_record(tmp_path)
    genome = tmp_path / "genome.fa"
    genome.write_text(">c1\n" + "A" * 1000 + "\n>c2\n" + "A" * 1000 + "\n")

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    def fake_genomic(*args, **kwargs):
        return []

    outcome = run_pipeline(
        genome_fasta=genome, proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=fake_genomic,
    )
    assert len(outcome.results) == 1
    result = outcome.results[0]
    assert result.fragmented is True
    assert sorted(result.genes_found) == ["mfa1", "pra1"]
    assert [(s.contig, s.start, s.end) for s in result.segments] == [("c1", 100, 200), ("c2", 300, 400)]
    # contig_edge_distance is populated from the real genome FASTA
    assert result.segments[0].contig_edge_distance == 99
    assert result.segments[1].contig_edge_distance == 299
    # would have been "high" on one contig; multi-segment downgrades one tier
    assert result.confidence == "medium"


def test_same_family_on_two_contigs_each_complete_is_not_fragmented(tmp_path):
    """Gene duplication / multi-allele co-occurrence is normal at MAT loci.
    Two COMPLETE copies of a family on different contigs are two real loci,
    never one fragmented locus."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c2", 100, 200, "+", 95.0, "rec1", "diamond_proteome"),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c2", 300, 400, "+", 95.0, "rec1", "diamond_proteome"),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=lambda *a, **k: [],
    )
    assert len(outcome.results) == 2
    assert all(r.fragmented is False for r in outcome.results)


def test_detection_result_carries_per_gene_evidence(tmp_path):
    """Finding 7: identity, coverage, coordinates, role and the matched
    curated record must survive to the DetectionResult, not be discarded."""
    _write_order(tmp_path)
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 91.5, "rec1",
                      "diamond_proteome", coverage=77.5),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "-", 88.0, "rec1",
                      "diamond_proteome", coverage=99.0),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=lambda *a, **k: [],
    )
    evidence = {e.gene_name: e for e in outcome.results[0].gene_evidence}
    assert evidence["mfa1"].identity == 91.5
    assert evidence["mfa1"].coverage == 77.5
    assert evidence["mfa1"].reference_record_id == "rec1"
    assert evidence["pra1"].strand == "-"
    assert evidence["pra1"].start == 300
    assert outcome.results[0].reference_records == ["rec1"]


def test_pipeline_output_feeds_the_report_writers_directly(tmp_path):
    """End-to-end boundary check: whatever run_pipeline returns must be exactly
    what report.py consumes, so "not detected" entries and per-gene evidence
    actually reach the CLI's files rather than stopping at an internal type."""
    import yaml

    from MATPredict.detect.report import write_detection_gff3, write_detection_report

    _write_order(
        tmp_path,
        ORDER_YML
        + "  - locus_name: bLocus\n    vocabulary_type: pattern\n"
        "    idiomorph_pattern: \"^b[0-9]+$\"\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: bE, role: core_MAT}\n      - {name: bW, role: core_MAT}\n",
    )
    _write_record(tmp_path)

    def fake_fast_path(proteome_fasta, families, reference_fasta, record_families, runner=None):
        return [
            SearchHit(FAMILY.key, "mfa1", "core_MAT", "c1", 100, 200, "+", 95.0, "rec1",
                      "diamond_proteome", coverage=80.0),
            SearchHit(FAMILY.key, "pra1", "core_MAT", "c1", 300, 400, "+", 93.0, "rec1",
                      "diamond_proteome", coverage=90.0),
        ]

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa", taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=fake_fast_path, search_genomic=lambda *a, **k: [],
    )
    write_detection_gff3(outcome, tmp_path / "out.gff3")
    write_detection_report(outcome, tmp_path / "out.yaml")

    gff3 = (tmp_path / "out.gff3").read_text()
    assert "Name=mfa1" in gff3 and "identity=95.0" in gff3
    doc = yaml.safe_load((tmp_path / "out.yaml").read_text())
    assert doc["detected"][0]["gene_evidence"][0]["coverage"] == 80.0
    # the family that never cleared the floor is reported, not dropped
    assert [n["family"] for n in doc["not_detected"]] == ["P:bLocus"]
