"""A PR call admitted only through a CAAX-scan precursor, in a family where the
CAAX-gained calls showed low enrichment for the curated mating-receptor clades,
is unverified.

Curator's ruling 2026-09-27 (results/2026-09-27_caax_precursor/per_family.tsv,
gained.tsv): gained calls whose receptor's nearest neighbour fell in the
curated mating-receptor subclades were Agrocybaceae 8/16, Mycenaceae 2/6,
Physalacriaceae 2/5, Galerinaceae 0/2 -- against 62% overall. The taxa are
curated data (`db/caax_unverified_taxa.yml`), not code. The label never changes
confidence; the call is still reported, and the rollout summary counts it.
"""
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult
from MATPredict.detect.rollout_aggregate import aggregate_reports
from MATPredict.detect.verification import (
    CAAX_UNVERIFIED_FILE, label_caax_unverified, load_caax_unverified_rules,
)

PR = FamilyKey("Basidiomycota", "PR")
AGROCYBACEAE = 3710399

RULES_YML = """\
caax_unverified_taxa:
  - taxid: 3710399
    name: Agrocybaceae
    reason: low enrichment of CAAX-gained receptors in mating-receptor clades (8/16)
    evidence: results/2026-09-27_caax_precursor/per_family.tsv
"""


def _r(caax_dependent=True, verification=None, conf="medium"):
    return DetectionResult(
        family_key=PR, contig="c1", start=1, end=100, confidence=conf,
        idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pheromone_receptor", "caax_precursor"], genes_missing=[],
        fragmented=False, caax_dependent=caax_dependent, verification=verification,
    )


def _rules(tmp_path):
    (tmp_path / CAAX_UNVERIFIED_FILE).write_text(RULES_YML)
    return load_caax_unverified_rules(tmp_path)


def test_a_missing_file_means_no_rules(tmp_path):
    assert load_caax_unverified_rules(tmp_path) == []


def test_the_real_database_lists_the_four_families():
    from pathlib import Path
    db = Path(__file__).resolve().parents[2] / "db"
    rules = load_caax_unverified_rules(db)
    assert {r.name for r in rules} == {"Agrocybaceae", "Mycenaceae", "Physalacriaceae", "Galerinaceae"}
    assert {r.taxid for r in rules} == {3710399, 2024004, 862241, 3710460}
    assert all(r.evidence.startswith("results/2026-09-27_caax_precursor/") for r in rules)


def test_a_caax_dependent_call_in_a_listed_family_is_unverified(tmp_path):
    out = label_caax_unverified([_r()], 999, _rules(tmp_path), lambda t: [1, AGROCYBACEAE, 5338])
    v = out[0].verification
    assert v["status"] == "unverified"
    assert v["taxid"] == AGROCYBACEAE and v["taxon"] == "Agrocybaceae"
    assert "enrichment" in v["reason"]
    assert out[0].confidence == "medium"


def test_the_genome_taxid_itself_may_be_listed(tmp_path):
    calls = []
    out = label_caax_unverified([_r()], AGROCYBACEAE, _rules(tmp_path),
                                lambda t: calls.append(t) or [])
    assert out[0].verification["status"] == "unverified"
    assert calls == []  # no lineage lookup needed


def test_a_call_not_dependent_on_the_scan_is_untouched(tmp_path):
    r = _r(caax_dependent=False)
    assert label_caax_unverified([r], 999, _rules(tmp_path), lambda t: [AGROCYBACEAE]) == [r]


def test_an_unlisted_family_is_untouched(tmp_path):
    r = _r()
    assert label_caax_unverified([r], 999, _rules(tmp_path), lambda t: [1, 5338]) == [r]


def test_no_taxid_or_no_rules_is_untouched(tmp_path):
    r = _r()
    assert label_caax_unverified([r], None, _rules(tmp_path), lambda t: [AGROCYBACEAE]) == [r]
    assert label_caax_unverified([r], 999, [], lambda t: [AGROCYBACEAE]) == [r]


def test_an_existing_label_is_kept(tmp_path):
    prior = {"status": "unverified", "reason": "override route"}
    r = _r(verification=prior)
    assert label_caax_unverified([r], 999, _rules(tmp_path), lambda t: [AGROCYBACEAE])[0].verification == prior


def test_a_failed_lineage_lookup_leaves_the_call_unlabelled(tmp_path):
    def boom(t):
        raise RuntimeError("efetch down")
    r = _r()
    assert label_caax_unverified([r], 999, _rules(tmp_path), boom) == [r]


def test_the_rollout_summary_counts_the_label(tmp_path):
    rep = {
        "genome_id": "g1", "taxid": 999,
        "detected": [{"family": "Basidiomycota:PR", "confidence": "medium",
                      "locus_class": "idiomorph_gene_only",
                      "verification": {"status": "unverified", "reason": "x"}}],
        "not_detected": [], "families_attempted": ["Basidiomycota:PR"],
    }
    p = tmp_path / "g1" / "detection_report.yaml"
    p.parent.mkdir()
    p.write_text(yaml.safe_dump(rep))
    summary = aggregate_reports([p], lineage_resolver=lambda t: "Basidiomycota")
    assert summary.unverified_calls == {"g1": 1}


def test_admission_only_through_the_scan():
    from MATPredict.detect.caax import admitted_only_through_scan
    scan = {"caax_precursor"}
    # receptor + precursor: 1 distinct homology gene -> dependent
    assert admitted_only_through_scan({"pheromone_receptor", "caax_precursor"}, scan, 1, 2)
    # two homology genes, both modelled -> not dependent
    assert not admitted_only_through_scan(
        {"pheromone_receptor", "fungal_mating_type_pheromone", "caax_precursor"}, scan, 2, 2)
    # two homology genes but only one modelled -> the modelled bar needs the scan
    assert admitted_only_through_scan(
        {"pheromone_receptor", "fungal_mating_type_pheromone", "caax_precursor"}, scan, 1, 2)
    # no scan precursor -> never dependent
    assert not admitted_only_through_scan({"pheromone_receptor"}, set(), 1, 2)
