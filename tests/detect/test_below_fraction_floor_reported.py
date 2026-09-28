"""A modelled cluster dropped at the fraction floor is reported, not lost.

Review finding F2 (2026-09-28, results/2026-09-28_fable_review/): the curated
Syncephalastrum racemosum NRRL 2496 record's own locus was admitted, polished
at 100%, then dropped at `ambiguity_floor` (2 of 5 live genes = 0.4 < 0.5) by
a bare `continue`, so it appeared in neither `detected`, `suppressed_loci` nor
`not_detected`. Such a cluster is now listed in `suppressed_loci` with
`withheld_reason: below_fraction_floor` and its fraction and best identity.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import BELOW_FRACTION_FLOOR, run_pipeline
import yaml

from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

KEY = FamilyKey("P", "MAT")


def _order(n_flanks):
    return (
        "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
        "    idiomorph_values: [Plus, Minus]\n    taxonomic_scope: [1]\n    genes:\n"
        + "".join(f"      - {{name: f{i}, role: flanking_variable}}\n" for i in range(1, n_flanks + 1))
        + "      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
    )


def _run(tmp_path, *, model_genes, relaxed=True):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(_order(4))
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text(
        "record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n")

    def localize(*a, **k):
        return [SearchHit(KEY, "sexP", "core_MAT", "c1", 1000, 1600, "+", 100.0, "rec1", "tblastn_genome"),
                SearchHit(KEY, "f1", "flanking_variable", "c1", 3000, 4000, "+", 99.0, "rec1", "tblastn_genome")]

    spans = {"sexP": (1000, 1600), "f1": (3000, 4000)}

    def model(gene_name, method):
        if gene_name not in model_genes:
            return None
        return _model(gene_name, "c1", *spans[gene_name], family_key=KEY,
                      role="core_MAT" if gene_name == "sexP" else "flanking_variable", method=method)

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
        relaxed_second_pass=relaxed,
    )


def test_a_modelled_cluster_below_the_floor_is_listed_as_suppressed(tmp_path):
    # 2 of 5 genes = 0.4 < 0.5; relaxed pass off, so nothing else reports it.
    outcome = _run(tmp_path, model_genes={"sexP"}, relaxed=False)
    assert outcome.results == []
    floor = [s for s in outcome.suppressed_loci if s.withheld_reason == BELOW_FRACTION_FLOOR]
    assert len(floor) == 1
    assert floor[0].contig == "c1"
    assert floor[0].withheld_detail["fraction_found"] == 0.4
    assert floor[0].withheld_detail["best_identity"] == 100.0


def test_the_report_carries_the_floor_entry(tmp_path):
    outcome = _run(tmp_path, model_genes={"sexP"}, relaxed=False)
    out = tmp_path / "report.yaml"
    write_detection_report(outcome, out)
    report = yaml.safe_load(out.read_text())
    (row,) = [s for s in report["suppressed_loci"] if s["withheld_reason"] == BELOW_FRACTION_FLOOR]
    assert row["fraction_found"] == 0.4


def test_an_unmodelled_cluster_below_the_floor_is_not_listed(tmp_path):
    outcome = _run(tmp_path, model_genes=set(), relaxed=False)
    assert not [s for s in outcome.suppressed_loci if s.withheld_reason == BELOW_FRACTION_FLOOR]


def test_a_floor_cluster_later_reported_by_the_relaxed_pass_is_not_also_suppressed(tmp_path):
    outcome = _run(tmp_path, model_genes={"sexP", "f1"}, relaxed=True)
    assert [r.detection_pass for r in outcome.results] == ["relaxed"]
    assert not [s for s in outcome.suppressed_loci if s.withheld_reason == BELOW_FRACTION_FLOOR]
