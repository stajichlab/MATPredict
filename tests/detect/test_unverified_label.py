"""A call made by searching a family outside its curated phylum is unverified.

Curator's ruling 2026-09-26: the Mortierellomycota and Kickxellomycota calls
(made with `--phylum Mucoromycota`) are "unverified -- not clear what they are
yet". The flank-ortholog synteny check (results/2026-09-26_flank_synteny_ED/
NOTE.md) found the Mucorales tptA-HMG-rnhA arrangement in 0/100 and 0/190 of
those genomes, against 117/293 in the Mucoromycota control, so none of their
13 calls is supported. The label is `verification: {status: unverified, ...}`
on the call; its confidence is not changed.

When is the family "outside its curated phylum"? Only an override route can
put it there -- `explicit_phylum` (`--phylum`) or `exhaustive`. Taxid routing
(`direct`, `lineage`, `phylum_fallback`) searches the genome's own phylum by
construction.

* override route, genome phylum known and different: unverified;
* override route, genome phylum the same: no label;
* `exhaustive`, genome phylum unknown: unverified -- nothing ties the genome
  to any curated phylum;
* `explicit_phylum`, genome phylum unknown: no label -- the operator stated
  the phylum, and there is nothing to contradict it.
"""
import yaml

from MATPredict.detect.family_registry import FamilyKey, RoutingDecision, load_all_families
from MATPredict.detect.pipeline import DetectionResult, run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.rollout_aggregate import aggregate_reports
from MATPredict.detect.search import SearchHit
from MATPredict.detect.verification import label_verification

from tests.detect.test_pipeline import _no_polish, _write_order, _write_record

MUC = FamilyKey("Mucoromycota", "MAT")


def _result(key=MUC):
    return DetectionResult(key, "c1", 1, 1000, "medium", "Minus", [], ["sexM"], [], False)


def test_an_override_on_another_phylum_is_unverified():
    [r] = label_verification([_result()], "explicit_phylum", "Kickxellomycota")
    assert r.verification["status"] == "unverified"
    assert r.verification["family_phylum"] == "Mucoromycota"
    assert r.verification["genome_phylum"] == "Kickxellomycota"
    assert "flank_synteny_ED" in r.verification["evidence"]
    assert r.confidence == "medium"


def test_an_override_on_the_same_phylum_is_not_labelled():
    [r] = label_verification([_result()], "explicit_phylum", "Mucoromycota")
    assert r.verification is None


def test_an_explicit_phylum_with_unknown_genome_phylum_is_not_labelled():
    [r] = label_verification([_result()], "explicit_phylum", None)
    assert r.verification is None


def test_an_exhaustive_search_with_unknown_genome_phylum_is_unverified():
    [r] = label_verification([_result()], "exhaustive", None)
    assert r.verification["status"] == "unverified"
    assert r.verification["genome_phylum"] is None


def test_taxid_routing_is_never_labelled():
    for mode in ("direct", "lineage", "phylum_fallback"):
        [r] = label_verification([_result()], mode, "Kickxellomycota")
        assert r.verification is None, mode


ORDER = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: enum\n"
    "    idiomorph_values: [a1, a2]\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT, present_in_idiomorphs: [a1]}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
)


def _run(tmp_path, genome_phylum=None, resolver=None):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    key = FamilyKey("P", "aLocus")
    hits = [SearchHit(key, "mfa1", "core_MAT", "c1", 1000, 1300, "+", 95.0, "rec1",
                      "diamond_proteome", coverage=90.0),
            SearchHit(key, "flk1", "flanking_conserved", "c1", 3000, 3300, "+", 95.0,
                      "rec1", "diamond_proteome", coverage=90.0)]
    routing = RoutingDecision(families=load_all_families(tmp_path),
                              routing_mode="explicit_phylum", phylum="P")
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=7, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: hits, search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        genetic_code=1, zygosity_lineage_resolver=lambda t: [], routing=routing,
        phylum_name_resolver=resolver or (lambda t: genome_phylum),
    )


def test_the_pipeline_labels_and_the_report_writes_it(tmp_path):
    outcome = _run(tmp_path, "Q")
    [r] = outcome.results
    assert r.verification["status"] == "unverified"
    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["detected"][0]["verification"]["status"] == "unverified"


def test_the_pipeline_does_not_label_a_same_phylum_call(tmp_path):
    [r] = _run(tmp_path, "P").results
    assert r.verification is None


def test_a_failed_phylum_lookup_leaves_an_explicit_call_unlabelled(tmp_path):
    """A resolver that raises is handled like an unknown phylum."""
    def broken(taxid):
        raise RuntimeError("efetch 400")
    [r] = _run(tmp_path, resolver=broken).results
    assert r.verification is None


def test_the_rollout_counts_unverified_calls(tmp_path):
    base = {"routing_mode": "explicit_phylum", "families_attempted": ["Mucoromycota:MAT"],
            "not_detected": []}
    call = {"family": "Mucoromycota:MAT", "confidence": "medium", "locus_class": "mat_locus"}
    docs = {
        "kick_a": dict(base, detected=[dict(call, verification={"status": "unverified"})]),
        "muc_b": dict(base, detected=[dict(call, verification=None)]),
    }
    paths = []
    for name, doc in docs.items():
        (tmp_path / name).mkdir()
        p = tmp_path / name / "detection_report.yaml"
        p.write_text(yaml.safe_dump(doc))
        paths.append(p)
    summary = aggregate_reports(paths, lineage_resolver=lambda taxid: None)
    assert summary.unverified_calls == {"kick_a": 1}
    assert summary.to_doc()["unverified_calls"] == {"kick_a": 1}
