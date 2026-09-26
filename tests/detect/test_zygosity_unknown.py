"""A single-idiomorph call in a taxon whose assemblies collapse MTL
heterozygosity is reported with zygosity unknown.

Curator's ruling 2026-09-26 (results/2026-09-26_calbicans_mtl_reads/NOTE.md).
In 16 C. albicans isolates, 8 of the 11 that are a/alpha by read depth had an
assembly holding only one idiomorph, and no assembly showed both when the
reads showed one. An assembly-based single-idiomorph call in this species says
nothing about zygosity. The taxa are curated data (`db/assembly_zygosity.yml`),
not code. The idiomorph calls themselves are unchanged.
"""
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult
from MATPredict.detect.rollout_aggregate import aggregate_reports
from MATPredict.detect.zygosity import (
    ZYGOSITY_FILE, genome_zygosity, load_zygosity_rules,
)

KEY = FamilyKey("Ascomycota", "MTL")

RULES_YML = """\
assemblies_collapse_heterozygosity:
  - taxid: 5476
    name: Candida albicans
    reason: assemblies collapse heterozygous MTL haplotypes
    evidence: results/2026-09-26_calbicans_mtl_reads/NOTE.md
"""


def _result(idiomorph):
    return DetectionResult(
        family_key=KEY, contig="c1", start=1, end=100, confidence="high",
        idiomorph=idiomorph, ambiguous_with=[], genes_found=["MTLA1"],
        genes_missing=[], fragmented=False,
    )


def _rules(tmp_path):
    (tmp_path / ZYGOSITY_FILE).write_text(RULES_YML)
    return load_zygosity_rules(tmp_path)


def test_a_missing_file_means_no_rules(tmp_path):
    assert load_zygosity_rules(tmp_path) == []


def test_the_real_database_lists_candida_albicans():
    from pathlib import Path
    db = Path(__file__).resolve().parents[2] / "db"
    rules = load_zygosity_rules(db)
    assert [r.taxid for r in rules] == [5476]
    assert rules[0].evidence.endswith("NOTE.md")


def test_one_idiomorph_in_a_listed_species_is_unknown(tmp_path):
    z = genome_zygosity([_result("A")], 5476, _rules(tmp_path), lambda t: [])
    assert z["status"] == "unknown"
    assert z["taxid"] == 5476
    assert "collapse" in z["reason"]
    assert z["evidence"] == "results/2026-09-26_calbicans_mtl_reads/NOTE.md"


def test_a_strain_taxid_matches_through_its_lineage(tmp_path):
    z = genome_zygosity([_result("alpha")], 237561, _rules(tmp_path),
                        lambda t: [131567, 4892, 5476])
    assert z["status"] == "unknown"


def test_two_loci_of_one_idiomorph_are_still_one_idiomorph(tmp_path):
    z = genome_zygosity([_result("A"), _result("A")], 5476, _rules(tmp_path), lambda t: [])
    assert z["status"] == "unknown"


def test_both_idiomorphs_are_left_alone(tmp_path):
    assert genome_zygosity([_result("A"), _result("alpha")], 5476, _rules(tmp_path),
                           lambda t: []) is None


def test_an_undetermined_call_does_not_count_as_an_idiomorph(tmp_path):
    assert genome_zygosity([_result("undetermined")], 5476, _rules(tmp_path),
                           lambda t: []) is None


def test_an_unlisted_species_is_left_alone(tmp_path):
    assert genome_zygosity([_result("A")], 498019, _rules(tmp_path),
                           lambda t: [131567, 4892, 498019]) is None


def test_no_call_means_no_zygosity_claim(tmp_path):
    assert genome_zygosity([], 5476, _rules(tmp_path), lambda t: []) is None


def test_a_failed_lineage_lookup_is_recorded_not_guessed(tmp_path):
    def broken(taxid):
        raise RuntimeError("efetch 400")
    z = genome_zygosity([_result("A")], 237561, _rules(tmp_path), broken)
    assert z["status"] == "unchecked"
    assert "efetch 400" in z["reason"]


def test_the_pipeline_and_report_carry_it(tmp_path):
    from MATPredict.detect.pipeline import run_pipeline
    from MATPredict.detect.report import write_detection_report
    from MATPredict.detect.search import SearchHit
    from tests.detect.test_pipeline import _no_polish, _write_order, _write_record

    (tmp_path / ZYGOSITY_FILE).write_text(RULES_YML.replace("5476", "1"))
    order = (
        "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: enum\n"
        "    idiomorph_values: [a1, a2]\n    taxonomic_scope: [1]\n"
        "    genes:\n      - {name: mfa1, role: core_MAT, present_in_idiomorphs: [a1]}\n"
        "      - {name: flk1, role: flanking_conserved}\n"
    )
    _write_order(tmp_path, order)
    _write_record(tmp_path)
    key = FamilyKey("P", "aLocus")
    hits = [SearchHit(key, "mfa1", "core_MAT", "c1", 1000, 1300, "+", 95.0, "rec1",
                      "diamond_proteome", coverage=90.0),
            SearchHit(key, "flk1", "flanking_conserved", "c1", 3000, 3300, "+", 95.0,
                      "rec1", "diamond_proteome", coverage=90.0)]
    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=1, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: hits, search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        genetic_code=1, zygosity_lineage_resolver=lambda t: [],
    )
    [r] = outcome.results
    assert r.idiomorph == "a1"
    assert outcome.zygosity["status"] == "unknown"

    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["zygosity"]["status"] == "unknown"
    assert doc["detected"][0]["idiomorph"] == "a1"


def test_the_report_writes_null_when_no_rule_applies(tmp_path):
    from MATPredict.detect.pipeline import DetectionOutcome
    from MATPredict.detect.report import write_detection_report
    path = tmp_path / "report.yaml"
    write_detection_report(DetectionOutcome(results=[]), path)
    doc = yaml.safe_load(path.read_text())
    assert "zygosity" in doc and doc["zygosity"] is None


def test_the_rollout_lists_zygosity_unknown_genomes(tmp_path):
    base = {"routing_mode": "lineage", "families_attempted": ["Ascomycota:MTL"],
            "detected": [], "not_detected": []}
    docs = {"5476_a": dict(base, zygosity={"status": "unknown"}),
            "5476_b": dict(base, zygosity=None)}
    paths = []
    for name, doc in docs.items():
        (tmp_path / name).mkdir()
        p = tmp_path / name / "detection_report.yaml"
        p.write_text(yaml.safe_dump(doc))
        paths.append(p)
    summary = aggregate_reports(paths, lineage_resolver=lambda taxid: None)
    assert summary.zygosity_unknown == ["5476_a"]
    assert summary.to_doc()["zygosity_unknown"] == ["5476_a"]
