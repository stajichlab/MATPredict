"""Which loser of a resolved sexM/sexP pair the locus verdict sets aside.

0d216c4 (#64) made the profile-HMM locus verdict score "the winners of resolved
pairs only", by setting aside every gene that lost a pair ANYWHERE in the
cluster. The clean-run regression panel (2026-10-22) found three failures:

(a) R. pusillus JY54 / Radiomyces spectabilis: sexM and sexP each lose at a
    different HMG gene, so both polished models were set aside, the verdict fell
    back to HSP fragments and the MAT-gene gate withheld the Plus call;
(b) M. irregularis B50: a 48-aa first-pass sexP tblastn fragment "won" inside
    the sexM gene and set aside the polished sexM model, so Minus/high became
    undetermined;
(c) Zygo23 scaffold (proteome) inputs: the winner is an annotated protein the
    classifier never scores, the loser's model was set aside, and the classifier
    no longer ran (Cunninghamella Minus high -> medium).

Curator's ruling 2026-10-22: set a loser aside ONLY where it overlaps its winner
AND the winner has a scored protein there; never discard both genes of a pair;
a first-pass fragment must not displace a polished model; proteome inputs keep a
scored protein. The original intent (a loser model overlapping the winner's
model is not scored) still holds.
"""
import dataclasses

from MATPredict.detect.idiomorph import IdiomorphResolution
from MATPredict.detect.pipeline import _pair_losers_set_aside, run_pipeline
from MATPredict.detect.search import SearchHit

from tests.detect.test_idiomorph_classifier import CODON, KEY, SEXM, SEXP, _setup
from tests.detect.test_pipeline import _model

ROLES = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved",
         "rnhA": "flanking_conserved"}
FLANKS = {"tptA": (12000, 12500), "rnhA": (13000, 13500)}


def _cds(protein, stop=True):
    return "".join(CODON[a] for a in protein) + ("TAA" if stop else "")


def _hit(gene, span, identity, method="tblastn_genome", **kw):
    return SearchHit(KEY, gene, ROLES[gene], "c1", *span, "+", identity, "rec1", method,
                     bitscore=identity, **kw)


def _run(tmp_path, genome, *, localize_hits, models, min_margin=10.0, flank_identity=90.0,
         proteome_hits=None):
    """Run the pipeline on `genome` (contig c1) with stubbed searches and polishing.

    `models[gene]` is a (start, end) span the polishing tools return for that gene,
    or None for a gene neither tool can model. Flanks are always found and modelled.
    """
    _setup(tmp_path, min_margin=min_margin)
    (tmp_path / "genome.fa").write_text(f">c1\n{genome.ljust(14000, 'A')}\n")
    flank_hits = [_hit(g, s, 90.0) for g, s in FLANKS.items()]
    spans = {**FLANKS, **models}

    def model(gene_name, method):
        span = spans.get(gene_name)
        if span is None:
            return None
        identity = flank_identity if ROLES[gene_name] != "core_MAT" else 60.0
        return _model(gene_name, "c1", *span, identity=identity, family_key=KEY,
                      role=ROLES[gene_name], method=method)

    kw = {}
    if proteome_hits is not None:
        (tmp_path / "proteome.faa").write_text(">p1\nM\n")
        kw["search_fast_path"] = lambda *a, **k: list(proteome_hits) + [
            _hit(g, s, 90.0, method="diamond_proteome") for g, s in FLANKS.items()]
        localize = lambda *a, **k: list(localize_hits)  # noqa: E731
    else:
        localize = lambda *a, **k: list(localize_hits) + flank_hits  # noqa: E731
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa",
        proteome_fasta=(tmp_path / "proteome.faa") if proteome_hits is not None else None,
        taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_localize=localize,
        polish_with_exonerate=lambda *, gene_name, **k: model(gene_name, "exonerate_refine"),
        polish_with_miniprot=lambda *, gene_name, **k: model(gene_name, "miniprot_refine"),
        **kw,
    )


def test_a_pair_that_loses_at_different_positions_keeps_both_models(tmp_path):
    """(a) JY54 / Radspe1 shape. Two HMG genes: a partial sexM-like gene and a
    full sexP gene. At each, both names hit and the right one wins the first
    pass, so sexM and sexP each lose somewhere. Flanks modelled below 40%, so a
    fragment-typed verdict is withheld by the MAT-gene gate: only scoring the
    polished models gives a call."""
    m_cds, p_cds = _cds(SEXM[: len(SEXM) // 3]), _cds(SEXP)
    a0 = 1000
    b0 = a0 + len(m_cds) + 500
    genome = "A" * (a0 - 1) + m_cds + "A" * 500 + p_cds
    m_span, p_span = (a0, a0 + len(m_cds) - 4), (b0, b0 + len(p_cds) - 4)
    out = _run(tmp_path, genome, flank_identity=35.0,
               localize_hits=[_hit("sexM", m_span, 60.0), _hit("sexP", m_span, 30.0),
                              _hit("sexP", p_span, 60.0), _hit("sexM", p_span, 30.0)],
               models={"sexM": m_span, "sexP": p_span})
    assert [s.withheld_reason for s in out.suppressed_loci] == []
    (r,) = out.results
    c = r.idiomorph_classifier
    assert c["classifier_input"] == "model"
    assert r.idiomorph == "Plus"
    assert {"sexM", "sexP"} <= set(r.genes_found)


def test_a_first_pass_fragment_does_not_displace_the_polished_model(tmp_path):
    """(b) B50 shape. One gene encoding the sexM protein, modelled full-length.
    A short sexP tblastn fragment inside it has the higher identity and wins the
    first pass; sexP itself cannot be modelled. The polished sexM model must
    still be scored (Minus, decisive); the fragment's HSP verdict alone is under
    the margin (on main: undetermined, and the spurious sexP joined the genes)."""
    filler = "DEKTLRVNPQ" * 4      # 40 aa inside the gene that neither profile scores
    m_cds = _cds(SEXM[:100] + filler + SEXM[100:])
    a0 = 1000
    genome = "A" * (a0 - 1) + m_cds
    m_span = (a0, a0 + len(m_cds) - 4)
    frag = (a0 + 3 * 100, a0 + 3 * (100 + len(filler)) - 1)
    out = _run(tmp_path, genome, min_margin=25.0,
               localize_hits=[_hit("sexM", m_span, 30.0), _hit("sexP", frag, 37.5)],
               models={"sexM": m_span, "sexP": None})
    (r,) = out.results
    c = r.idiomorph_classifier
    assert c["classifier_input"] == "model"
    assert r.idiomorph == "Minus"
    assert c["margin"] >= 25.0
    assert "sexP" not in r.genes_found
    assert r.confidence == "high"


def test_a_proteome_winner_keeps_the_classifier_running(tmp_path):
    """(c) Zygo23 scaffold shape. The annotated protein (diamond_proteome) wins
    sexM; the opposite name is found by tblastn rescue inside it and modelled.
    The annotated protein is never classifier input, so the loser's model is
    the only scored protein and must stay: the classifier still runs and the
    call keeps its classifier verdict and confidence."""
    m_cds = _cds(SEXM)
    a0 = 1000
    genome = "A" * (a0 - 1) + m_cds
    m_span = (a0, a0 + len(m_cds) - 4)
    out = _run(tmp_path, genome,
               proteome_hits=[_hit("sexM", m_span, 90.0, method="diamond_proteome", coverage=95.0)],
               localize_hits=[_hit("sexP", m_span, 30.0)],
               models={"sexP": m_span})
    (r,) = out.results
    c = r.idiomorph_classifier
    assert c is not None and c["idiomorph"] == "Minus"
    assert r.idiomorph == "Minus"
    assert r.confidence == "high"


def test_a_loser_model_overlapping_the_winner_model_is_still_set_aside(tmp_path):
    """The case 0d216c4 was written for, kept: the sexM model overlaps the
    winning sexP model (the pair is decided on profile bitscores, sexP wins) and
    reads on into a sexM-like protein. Its protein must not enter the verdict:
    the Minus score is the winner's own protein against the sexM profile."""
    p_cds = _cds(SEXP, stop=False)
    tail = _cds(SEXM[: len(SEXM) // 4])
    a0 = 1000
    genome = "A" * (a0 - 1) + p_cds + tail
    p_span = (a0, a0 + len(p_cds) - 1)
    m_span = (a0 + len(p_cds) - 3 * 60, a0 + len(p_cds) + len(tail) - 4)
    out = _run(tmp_path, genome,
               localize_hits=[_hit("sexP", p_span, 60.0), _hit("sexM", m_span, 40.0)],
               models={"sexP": p_span, "sexM": m_span})
    (r,) = out.results
    (e,) = [e for e in r.idiomorph_resolutions if {e.winner, e.loser} == {"sexM", "sexP"}]
    assert (e.winner, e.loser) == ("sexP", "sexM")
    c = r.idiomorph_classifier
    assert c["classifier_input"] == "model" and r.idiomorph == "Plus"
    # Pooled, the loser's sexM-like read-through would score well against sexM.
    loser_only = _run(tmp_path / "x", "A" * (a0 - 1) + tail,
                      localize_hits=[_hit("sexM", (a0, a0 + len(tail) - 4), 40.0)],
                      models={"sexM": (a0, a0 + len(tail) - 4)})
    (lr,) = loser_only.results
    assert c["scores"]["Minus"] < lr.idiomorph_classifier["scores"]["Minus"] - 10


# --- the rule itself, on synthetic entries ---------------------------------

def _event(winner, loser, wspan, lspan):
    return IdiomorphResolution(contig="c1", winner=winner, loser=loser, winner_identity=60.0,
                               loser_identity=30.0, overlap_fraction=1.0,
                               winner_span=wspan, loser_span=lspan)


def test_rule_sets_aside_only_where_the_winner_has_a_scored_protein():
    entries = [("sexP", ("c1", 100, 400), "P"), ("sexM", ("c1", 150, 450), "M")]
    ev = [_event("sexP", "sexM", (100, 400), (150, 450))]
    assert _pair_losers_set_aside(entries, ev) == {1}
    # Winner unmodelled (fragment or proteome protein): the loser stays.
    assert _pair_losers_set_aside(entries[1:], ev) == set()


def test_rule_is_positional():
    # sexM lost at 2000-2300, not where its model is.
    entries = [("sexP", ("c1", 100, 400), "P"), ("sexM", ("c1", 150, 450), "M"),
               ("sexP", ("c1", 2000, 2300), "P2")]
    ev = [_event("sexP", "sexM", (2000, 2300), (2000, 2300))]
    assert _pair_losers_set_aside(entries, ev) == set()


def test_rule_never_discards_both_genes_of_a_pair():
    entries = [("sexP", ("c1", 100, 400), "P"), ("sexM", ("c1", 150, 450), "M")]
    ev = [_event("sexP", "sexM", (100, 400), (150, 450)),
          _event("sexM", "sexP", (150, 450), (100, 400))]
    assert _pair_losers_set_aside(entries, ev) == set()


def test_rule_does_nothing_without_events():
    entries = [("sexP", ("c1", 100, 400), "P"), ("sexM", ("c1", 150, 450), "M")]
    assert _pair_losers_set_aside(entries, []) == set()


# --- Fix A: a first-pass-identity winner must be (nearly) as long as the loser ----

def _based(winner, loser, wspan, lspan, basis):
    return dataclasses.replace(_event(winner, loser, wspan, lspan), basis=basis)


def test_a_short_first_pass_winner_keeps_the_full_length_loser():
    """Fable review of #140, finding 3: no tool modelled both genes, so the HMM pass
    recorded nothing and the first-pass event stands. A 60-aa sexP model must not set
    aside the full-length sexM model it sits inside."""
    entries = [("sexP", ("c1", 1000, 1180), "P" * 60), ("sexM", ("c1", 400, 1900), "M" * 480)]
    ev = [_based("sexP", "sexM", (1000, 1180), (400, 1900), "first_pass_identity")]
    assert _pair_losers_set_aside(entries, ev) == set()


def test_a_model_level_decision_sets_aside_even_for_a_short_winner():
    """The profile bitscore is authoritative: the same short winner still sets the loser aside."""
    entries = [("sexP", ("c1", 1000, 1180), "P" * 60), ("sexM", ("c1", 400, 1900), "M" * 480)]
    for basis in ("hmm_bitscore", "model_score_miniprot", "model_score_exonerate"):
        ev = [_based("sexP", "sexM", (1000, 1180), (400, 1900), basis)]
        assert _pair_losers_set_aside(entries, ev) == {1}, basis


def test_a_long_first_pass_winner_still_sets_aside_a_short_loser():
    entries = [("sexP", ("c1", 400, 1900), "P" * 480), ("sexM", ("c1", 1000, 1180), "M" * 60)]
    ev = [_based("sexP", "sexM", (400, 1900), (1000, 1180), "first_pass_identity")]
    assert _pair_losers_set_aside(entries, ev) == {1}
    # The boundary: a winner of exactly 0.8 x the loser's length still displaces it.
    entries = [("sexP", ("c1", 100, 400), "P" * 80), ("sexM", ("c1", 150, 450), "M" * 100)]
    ev = [_based("sexP", "sexM", (100, 400), (150, 450), "first_pass_identity")]
    assert _pair_losers_set_aside(entries, ev) == {1}
    entries[0] = ("sexP", ("c1", 100, 400), "P" * 79)
    assert _pair_losers_set_aside(entries, ev) == set()


def test_the_guard_never_adds_a_set_aside_to_a_mutual_pair():
    """Fable review of the guard, finding 1: a long sexP and a short sexM that each lost to the
    other at the same place (two first-pass events). The unguarded rule keeps both; the guard
    must not turn that into 'the long sexP displaces the short sexM'."""
    entries = [("sexP", ("c1", 400, 1900), "P" * 480), ("sexM", ("c1", 900, 1100), "M" * 60)]
    ev = [_based("sexP", "sexM", (400, 1900), (900, 1100), "first_pass_identity"),
          _based("sexM", "sexP", (900, 1100), (400, 1900), "first_pass_identity")]
    assert _pair_losers_set_aside(entries, ev) == set()


def test_an_unlisted_basis_is_held_to_the_length_guard():
    """An idiomorph-resolution basis added later and not in `_MODEL_LEVEL_BASES` fails toward
    keeping the loser's model: a short winner under it must not set a full-length loser aside."""
    entries = [("sexP", ("c1", 1000, 1180), "P" * 60), ("sexM", ("c1", 400, 1900), "M" * 480)]
    ev = [_based("sexP", "sexM", (1000, 1180), (400, 1900), "some_future_basis")]
    assert _pair_losers_set_aside(entries, ev) == set()
    # ... and the classifier's own decision is listed as model-level.
    ev = [_based("sexP", "sexM", (1000, 1180), (400, 1900), "hmm_classifier")]
    assert _pair_losers_set_aside(entries, ev) == {1}
