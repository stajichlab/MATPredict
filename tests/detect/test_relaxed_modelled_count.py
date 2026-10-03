"""A relaxed-pass call carries its real modelled-gene count.

Defect found 2026-09-26 (results/2026-09-26_phycomyces_trace/NOTE.md): since
6bc985d the modelled-gene bar (`min_polished_genes`) is applied after the
relaxed pass, but `_relaxed_results` never set `polished_genes`, so it stayed
at its default 0 and EVERY relaxed call was withheld. Phycomyces blakesleeanus
NRRL_1554 on the contig input -- sexP + rnhA on one contig, both modelled --
lost its correct Plus call; 0 relaxed calls exist across 294 Mucoromycota and
2,368 Serinales reports since.

The fix counts modelled genes with the SAME rule as the strict path
(`_modelled_gene_names`), through one shared helper.
"""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.pipeline import EvidenceFloor, _relaxed_results, run_pipeline
from MATPredict.detect.polish import STATUS_AGREE, STATUS_UNPOLISHED, PolishOutcome
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

RICH = Family(
    FamilyKey("Ascomycota", "MAT"), "enum", ["MAT1-1", "MAT1-2"], None,
    [{"name": f"g{i}", "role": "core_MAT" if i == 1 else "flanking_variable"}
     for i in range(1, 9)],
    [4890],
)
SEARCHABLE = {RICH.key: {f"g{i}" for i in range(1, 9)}}


def _tb(gene, role, start, end):
    return SearchHit(RICH.key, gene, role, "c1", start, end, "+", 80.0, "rec1",
                     "tblastn_genome")


def _cluster():
    return GeneCluster("c1", 1, 1000, [_tb("g1", "core_MAT", 1, 100),
                                        _tb("g2", "flanking_variable", 400, 500)])


def _outcome(status, gene):
    m = _model(gene, "c1", 1, 100, family_key=RICH.key) if status != STATUS_UNPOLISHED else None
    return PolishOutcome(status, m, m, m)


def test_two_modelled_genes_are_counted():
    c = _cluster()
    polish_by = {(id(c), RICH.key, "g1"): _outcome(STATUS_AGREE, "g1"),
                 (id(c), RICH.key, "g2"): _outcome(STATUS_AGREE, "g2")}
    (r,) = _relaxed_results([c], [RICH], searchable_genes=SEARCHABLE,
                            evidence_floor=EvidenceFloor(), polish_by=polish_by)
    assert r.polished_genes == 2


def test_an_unpolished_gene_is_not_counted():
    c = _cluster()
    polish_by = {(id(c), RICH.key, "g1"): _outcome(STATUS_AGREE, "g1"),
                 (id(c), RICH.key, "g2"): _outcome(STATUS_UNPOLISHED, "g2")}
    (r,) = _relaxed_results([c], [RICH], searchable_genes=SEARCHABLE,
                            evidence_floor=EvidenceFloor(), polish_by=polish_by)
    assert r.polished_genes == 1


def test_fast_path_genes_count_as_modelled_as_in_the_strict_path():
    c = GeneCluster("c1", 1, 1000, [
        SearchHit(RICH.key, "g1", "core_MAT", "c1", 1, 100, "+", 90.0, "rec1", "diamond_proteome"),
        SearchHit(RICH.key, "g2", "flanking_variable", "c1", 400, 500, "+", 90.0, "rec1", "diamond_proteome"),
    ])
    (r,) = _relaxed_results([c], [RICH], searchable_genes=SEARCHABLE,
                            evidence_floor=EvidenceFloor())
    assert r.polished_genes == 2


ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [Plus, Minus]\n    taxonomic_scope: [1]\n    genes:\n"
    + "".join(f"      - {{name: f{i}, role: flanking_variable}}\n" for i in range(1, 7))
    + "      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
)


def _run(tmp_path, *, model_genes, diagnostics=None):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text(
        "record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n")
    key = FamilyKey("P", "MAT")

    def localize(*a, **k):
        return [SearchHit(key, "sexP", "core_MAT", "c1", 1000, 1600, "+", 90.0, "rec1", "tblastn_genome"),
                SearchHit(key, "f1", "flanking_variable", "c1", 3000, 4000, "+", 90.0, "rec1", "tblastn_genome")]

    def model(gene_name, method):
        if gene_name not in model_genes:
            return None
        span = {"sexP": (1000, 1600), "f1": (3000, 4000)}[gene_name]
        return _model(gene_name, "c1", *span, family_key=key,
                      role="core_MAT" if gene_name == "sexP" else "flanking_variable", method=method)

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
        evidence_diagnostics_path=diagnostics,
    )


def test_a_relaxed_call_with_two_modelled_genes_is_reported_end_to_end(tmp_path):
    """2 of 7 = 0.29: below the strict floor, so only the relaxed pass can call it."""
    outcome = _run(tmp_path, model_genes={"sexP", "f1"})
    (r,) = outcome.results
    assert r.detection_pass == "relaxed"
    assert r.polished_genes == 2
    assert r.idiomorph == "Plus"
    assert r.confidence == "medium"


def test_a_relaxed_call_with_one_modelled_gene_is_still_withheld(tmp_path):
    outcome = _run(tmp_path, model_genes={"sexP"})
    assert outcome.results == []
    assert outcome.suppressed_unpolished == 1


def test_each_relaxed_call_logs_its_uncapped_tier(tmp_path):
    """Exploration data for the curator's medium-cap question."""
    import json
    diag = tmp_path / "diag.jsonl"
    _run(tmp_path, model_genes={"sexP", "f1"}, diagnostics=diag)
    rows = [json.loads(l) for l in diag.read_text().splitlines()]
    (row,) = [r for r in rows if r["kind"] == "relaxed_call"]
    assert row["confidence_reported"] == "medium"
    assert row["tier_uncapped"] in ("high", "medium", "low")
    assert row["polished_genes"] == 2
