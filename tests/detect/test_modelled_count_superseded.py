"""A gene superseded AFTER polishing must not count toward the modelled bar.

The path is RESCUE polishing. A cluster that evidences no idiomorph yet gets
both idiomorphs' missing core genes rescued in its window. sexM and sexP share
an HMG box, so both rescues can model the SAME real gene. Each rescued model
is recorded as modelled, and only afterwards does the post-polish idiomorph
pass supersede the loser -- which was then still counted. The loser is the same physical gene as the winner, so counting both
lets ONE gene clear a bar of two.

Measured before the fix: 25 of 3,101 Saccharomyces loci carried a
`polished_genes` larger than their `genes_found` (MATALPHA2 superseded by
MATA1 after polishing, still counted). None crossed the bar only because of
it, but the count is wrong and nothing prevented the crossing.
"""
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

ORDER = (
    "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
    "    idiomorph_values: [Plus, Minus]\n    taxonomic_scope: [1]\n"
    "    genes:\n"
    "      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
    "      - {name: sexM, role: core_MAT, present_in_idiomorphs: [Minus]}\n"
    "      - {name: sexC, role: core_MAT, present_in_idiomorphs: [Plus, Minus]}\n"
    "      - {name: tptA, role: flanking_conserved}\n"
)


def _setup(tmp_path):
    (tmp_path / "P").mkdir()
    (tmp_path / "P" / "order.yml").write_text(ORDER)
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text(
        "record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n"
    )


def _tblastn(gene, start, end, identity, role="core_MAT"):
    from MATPredict.detect.family_registry import FamilyKey
    return SearchHit(FamilyKey("P", "MAT"), gene, role, "c1", start, end, "+",
                     identity, "rec1", "tblastn_genome", coverage=60.0)


def _run(tmp_path, *, flank_models, localized=None):
    """sexC (both idiomorphs, so no idiomorph is evidenced) and tptA localize;
    sexP and sexM are rescued and BOTH model the one gene at 1000-1600, where
    sexM loses. sexC is never modelled; tptA only when `flank_models`."""
    from MATPredict.detect.family_registry import FamilyKey
    key = FamilyKey("P", "MAT")
    _setup(tmp_path)
    spans = {"sexP": (1000, 1600), "sexM": (1000, 1600), "tptA": (3000, 4000)}
    ident = {"sexP": 60.0, "sexM": 35.0, "tptA": 80.0}
    roles = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved"}

    def localize(*a, **k):
        return localized or [_tblastn("sexC", 2000, 2600, 50.0),
                             _tblastn("tptA", 3000, 4000, 80.0, role="flanking_conserved")]

    def model(gene_name, method):
        if gene_name == "sexC" or (gene_name == "tptA" and not flank_models):
            return None
        return _model(gene_name, "c1", *spans[gene_name], identity=ident[gene_name],
                      family_key=key, role=roles[gene_name], method=method)

    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )


def test_two_real_genes_clear_the_bar(tmp_path):
    """Control: sexP and tptA are both modelled -- two genes, reported."""
    outcome = _run(tmp_path, flank_models=True)
    assert len(outcome.results) == 1
    assert outcome.results[0].polished_genes == 2


def test_a_superseded_gene_does_not_count_toward_the_bar(tmp_path):
    """Only sexP is a real modelled gene; sexM's model is the same gene."""
    outcome = _run(tmp_path, flank_models=False)
    assert outcome.results == []
    assert outcome.suppressed_unpolished == 1


def test_each_resolution_is_reported_once(tmp_path):
    """The pre- and post-polish passes must not both record the same event.

    sexP and sexM already overlap as raw hits, so the pre-polish pass resolves
    them, and the post-polish pass sees the same pair again.
    """
    outcome = _run(tmp_path, flank_models=True, localized=[
        _tblastn("sexP", 1000, 1600, 60.0), _tblastn("sexM", 1000, 1600, 35.0),
        _tblastn("tptA", 3000, 4000, 80.0, role="flanking_conserved"),
    ])
    events = [(e.contig, e.winner, e.loser) for e in outcome.results[0].idiomorph_resolutions]
    assert len(events) == len(set(events))
