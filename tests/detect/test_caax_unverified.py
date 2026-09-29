"""Every PR call admitted only through a CAAX-scan precursor is unverified.

Curator's ruling 2026-09-28 (review finding F4), replacing the 2026-09-27
four-family list. The negative control (results/2026-09-28_validation_f3_f4/
NOTE.md) could not bound the finder's false-positive rate: 6/9 curated mating
receptors flagged, 0/25 non-mating STE3 copies (95% upper limit 13.7%),
random windows 2.3%, so ~13 of the 118 gained Agaricales calls are expected
by chance alone and the true rate may be higher. The label never changes
confidence; the call is still reported. REVIEW LATER, once a labelled set of
>= 100 non-mating STE3 loci exists.
"""
import yaml

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult
from MATPredict.detect.rollout_aggregate import aggregate_reports
from MATPredict.detect.verification import CAAX_UNVERIFIED_EVIDENCE, label_caax_unverified

PR = FamilyKey("Basidiomycota", "PR")


def _r(caax_dependent=True, verification=None, conf="medium"):
    return DetectionResult(
        family_key=PR, contig="c1", start=1, end=100, confidence=conf,
        idiomorph="undetermined", ambiguous_with=[],
        genes_found=["pheromone_receptor", "caax_precursor"], genes_missing=[],
        fragmented=False, caax_dependent=caax_dependent, verification=verification,
    )


def test_every_caax_dependent_call_is_unverified():
    (r,) = label_caax_unverified([_r()])
    assert r.verification["status"] == "unverified"
    assert "strict-CAAX" in r.verification["reason"]
    assert r.verification["evidence"] == CAAX_UNVERIFIED_EVIDENCE
    assert r.confidence == "medium"


def test_the_evidence_is_the_validation_note():
    assert CAAX_UNVERIFIED_EVIDENCE == "results/2026-09-28_validation_f3_f4/NOTE.md"


def test_a_call_not_dependent_on_the_scan_is_untouched():
    r = _r(caax_dependent=False)
    assert label_caax_unverified([r]) == [r]


def test_an_existing_label_is_kept():
    existing = {"status": "unverified", "reason": "out-of-phylum override"}
    (r,) = label_caax_unverified([_r(verification=existing)])
    assert r.verification == existing


def test_the_four_family_list_is_gone():
    from pathlib import Path
    db = Path(__file__).resolve().parents[2] / "db"
    assert not (db / "caax_unverified_taxa.yml").exists()


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
