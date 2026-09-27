"""Model BOTH idiomorph-alternative genes, and decide on the two models.

Curator ruling, J. Stajich, 2026-09-26: sexM and sexP "are hard to tell apart",
so model both.

WHY. Before this, the pre-polish overlap resolution picked a winner from raw
tblastn identity, the loser was marked superseded, and a superseded hit was
never polished (`not_polish_candidate`) and never voted. In 13 Mucoromycota
Minus calls whose HMG box falls in the supported sexP clade of the
2026-09-26 HMG-box tree (Backusella x6, Apophysomyces x4, Blakeslea
trispora, Umbelopsis x2), that first-pass margin was 0.1-2 identity points
at ~34-36% -- the whole call rested on it, and the sexP gene was never given
a model to compete with.

With `model_idiomorph_alternatives: true` on a family, both genes of an
overlapping mutually exclusive pair are polished, and the pair is decided on
the tools' own alignment scores for the two models in the same window
(miniprot's score when both have a miniprot model, else exonerate's when
both have an exonerate model; the two tools' scores are never compared with
each other). The loser stays superseded, so it neither counts toward the
modelled-gene bar nor votes. No usable pair of scores: the first-pass
verdict stands.
"""
import dataclasses

from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _model

KEY = FamilyKey("P", "MAT")


def _order(flag):
    extra = "    model_idiomorph_alternatives: true\n" if flag else ""
    return (
        "phylum: P\nloci:\n  - locus_name: MAT\n    vocabulary_type: enum\n"
        "    idiomorph_values: [Plus, Minus]\n    taxonomic_scope: [1]\n" + extra +
        "    genes:\n"
        "      - {name: sexP, role: core_MAT, present_in_idiomorphs: [Plus]}\n"
        "      - {name: sexM, role: core_MAT, present_in_idiomorphs: [Minus]}\n"
        "      - {name: tptA, role: flanking_conserved}\n"
        "      - {name: rnhA, role: flanking_conserved}\n"
    )


def _setup(tmp_path, flag):
    (tmp_path / "P").mkdir(parents=True)
    (tmp_path / "P" / "order.yml").write_text(_order(flag))
    rec = tmp_path / "P" / "Fam" / "rec1"
    rec.mkdir(parents=True)
    rec.joinpath("metadata.yaml").write_text(
        "record_id: rec1\nmating_type: {locus_name: MAT, idiomorphs: [Plus]}\n"
    )


def _tblastn(gene, start, end, identity, bits, role="core_MAT"):
    return SearchHit(KEY, gene, role, "c1", start, end, "+", identity, "rec1",
                     "tblastn_genome", bitscore=bits)


ROLES = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved",
         "rnhA": "flanking_conserved"}
SPANS = {"sexP": (1000, 1600), "sexM": (1000, 1600), "tptA": (3000, 4000),
         "rnhA": (5000, 6000)}


def _run(tmp_path, *, flag, scores, model_fails=()):
    """sexM wins the first pass on raw identity (36.1 vs 34.0) over one HMG
    gene; `scores` gives each core gene's model score (same for both tools)."""
    _setup(tmp_path, flag)
    polished = []

    def localize(*a, **k):
        return [_tblastn("sexM", 1000, 1600, 36.1, 60.0),
                _tblastn("sexP", 1000, 1600, 34.0, 55.0),
                _tblastn("tptA", 3000, 4000, 90.0, 400.0, role="flanking_conserved"),
                _tblastn("rnhA", 5000, 6000, 90.0, 400.0, role="flanking_conserved")]

    def model(gene_name, method):
        polished.append((gene_name, method))
        if gene_name in model_fails:
            return None
        m = _model(gene_name, "c1", *SPANS[gene_name], identity=35.0, family_key=KEY,
                   role=ROLES[gene_name], method=method)
        return dataclasses.replace(m, score=scores.get(gene_name, 500.0))

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )
    return outcome, polished


def _evidence(result, gene):
    return next(e for e in result.gene_evidence if e.gene == gene)


def test_without_the_setting_the_first_pass_decides_and_the_loser_is_not_modelled(tmp_path):
    outcome, polished = _run(tmp_path, flag=False, scores={"sexM": 150.0, "sexP": 300.0})
    (r,) = outcome.results
    assert r.idiomorph == "Minus"
    assert "sexP" not in {g for g, _ in polished}


def test_with_the_setting_both_are_modelled_and_the_models_decide(tmp_path):
    outcome, polished = _run(tmp_path, flag=True, scores={"sexM": 150.0, "sexP": 300.0})
    (r,) = outcome.results
    assert {g for g, _ in polished} >= {"sexM", "sexP"}
    assert r.idiomorph == "Plus"
    assert "sexP" in r.genes_found and "sexM" not in r.genes_found
    (event,) = [e for e in r.idiomorph_resolutions if {e.winner, e.loser} == {"sexM", "sexP"}]
    assert (event.winner, event.loser) == ("sexP", "sexM")
    assert event.basis == "model_score_miniprot"
    assert (event.winner_model_score, event.loser_model_score) == (300.0, 150.0)


def test_the_models_can_confirm_the_first_pass(tmp_path):
    outcome, _ = _run(tmp_path, flag=True, scores={"sexM": 300.0, "sexP": 150.0})
    (r,) = outcome.results
    assert r.idiomorph == "Minus"
    (event,) = [e for e in r.idiomorph_resolutions if {e.winner, e.loser} == {"sexM", "sexP"}]
    assert event.winner == "sexM" and event.basis == "model_score_miniprot"


def test_the_loser_does_not_count_toward_the_modelled_bar(tmp_path):
    outcome, _ = _run(tmp_path, flag=True, scores={"sexM": 150.0, "sexP": 300.0})
    (r,) = outcome.results
    assert r.polished_genes == 3          # sexP, tptA, rnhA -- not sexM


def test_a_loser_with_no_model_keeps_the_first_pass_and_does_not_cap_the_tier(tmp_path):
    base, _ = _run(tmp_path / "a", flag=False, scores={})
    outcome, _ = _run(tmp_path / "b", flag=True, scores={}, model_fails=("sexP",))
    (r,) = outcome.results
    (b,) = base.results
    assert r.idiomorph == "Minus"
    assert r.confidence == b.confidence
    event = next(e for e in r.idiomorph_resolutions if {e.winner, e.loser} == {"sexM", "sexP"})
    assert event.basis == "first_pass_identity"


def test_the_setting_is_read_from_the_curated_mucoromycota_roster():
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families
    db = Path(__file__).resolve().parents[2] / "db"
    fams = {(f.key.phylum, f.key.locus_name): f for f in load_all_families(db)}
    assert fams[("Mucoromycota", "MAT")].model_idiomorph_alternatives is True
    assert not any(f.model_idiomorph_alternatives for k, f in fams.items()
                   if k != ("Mucoromycota", "MAT"))


def test_a_losing_model_never_counts_even_when_its_gene_has_another_live_hit(tmp_path):
    """Found in the 2026-09-26 validation run (Umbelopsis vinacea
    GCA_016758895.1): sexP lost the model decision at the HMG gene, but a
    separate weak sexP hit elsewhere in the cluster kept the NAME live, so the
    losing sexP model was counted as a modelled gene. One physical gene then
    cleared a bar of two. Here only sexM and nothing else is modelled, so the
    call must be withheld."""
    _setup(tmp_path, True)

    def localize(*a, **k):
        return [_tblastn("sexM", 1000, 1600, 36.1, 60.0),
                _tblastn("sexP", 1000, 1600, 34.0, 55.0),
                _tblastn("sexP", 2400, 2500, 30.0, 30.0),     # another HMG, sexP-labelled
                _tblastn("tptA", 3000, 4000, 90.0, 400.0, role="flanking_conserved")]

    def model(gene_name, method):
        if gene_name not in ("sexM", "sexP"):
            return None                                        # flank not modelled
        m = _model(gene_name, "c1", *SPANS[gene_name], identity=35.0, family_key=KEY,
                   role="core_MAT", method=method)
        return dataclasses.replace(m, score={"sexM": 300.0, "sexP": 150.0}[gene_name])

    outcome = run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=None, taxid=None,
        db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **kw: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **kw: model(gene_name, "miniprot_refine"),
    )
    assert outcome.results == []
    assert outcome.suppressed_unpolished == 1
