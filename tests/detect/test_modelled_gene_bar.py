"""A locus must rest on real gene models, not on bare alignments.

Curator's rulings, 2026-09-22, made against the Pezizomycotina panels
(46,647 lineage-routed loci):

1. "without polishing it is low" -- a locus whose every gene is a raw tblastn
   HSP no tool could model is `low`, not `medium`.
2. "same as 1, must have polished" -- `mat_locus` is the strongest claim the
   pipeline makes; 3,388 loci were making it on entirely unpolished evidence.
3. "yes require more than 1" -- two modelled genes to be reported at all.

Why two and not one: loci with exactly ONE modelled gene numbered 1,161 and
produced ZERO high-confidence calls, so the stricter bar is free. All 2,957
high-confidence calls carry two or more.

The bar counts ANNOTATED genes too (`diamond_proteome`), not just polished
ones. Gating on polish alone would withhold the real MAT locus of every fully
annotated genome, because an annotated gene is never put to the tools.
"""
from pathlib import Path

from MATPredict.detect.idiomorph import LOCUS_CLASS_MAT, LOCUS_CLASS_PARTIAL
from MATPredict.detect.pipeline import MIN_POLISHED_GENES, run_pipeline
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import (  # the established fixtures
    FAMILY, _no_polish, _write_order, _write_record,
)


def _raw(gene, start, end, ident=31.0):
    """A bare tblastn HSP: localized, and no tool could model it."""
    return SearchHit(FAMILY.key, gene, "core_MAT", "c1", start, end, "+", ident,
                     "rec1", "tblastn_genome", coverage=12.0)


def _annotated(gene, start, end):
    """A gene from the annotated fast path -- a model somebody already called."""
    return SearchHit(FAMILY.key, gene, "core_MAT", "c1", start, end, "+", 91.5,
                     "rec1", "diamond_proteome", coverage=88.0)


def _run(tmp_path, hits, **kw):
    _write_order(tmp_path)
    _write_record(tmp_path)
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: hits,
        # Stubbed so a cluster missing a core gene does not reach for a real
        # makeblastdb/tblastn in the genome-wide rescue.
        search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
        **kw,
    )


def test_the_default_bar_is_two():
    assert MIN_POLISHED_GENES == 2


def test_an_all_unpolished_locus_is_withheld(tmp_path):
    outcome = _run(tmp_path, [_raw("mfa1", 100, 200), _raw("pra1", 300, 400)])
    assert outcome.results == []
    assert outcome.suppressed_unpolished == 1


def test_a_withheld_family_still_says_why(tmp_path):
    outcome = _run(tmp_path, [_raw("mfa1", 100, 200), _raw("pra1", 300, 400)])
    reasons = [n.reason for n in outcome.not_detected if n.family_key == FAMILY.key]
    assert reasons, "a withheld family must still appear in not_detected"
    assert "modelled" in reasons[0]
    assert str(MIN_POLISHED_GENES) in reasons[0]


def test_an_annotated_locus_is_reported(tmp_path):
    """The regression that matters: annotated genomes must still work.

    Both genes come from diamond and are never polished. Gating on polish
    alone would withhold this, and with it every ZygoLife/BFD-proteome call.
    """
    outcome = _run(tmp_path, [_annotated("mfa1", 100, 200), _annotated("pra1", 300, 400)])
    assert len(outcome.results) == 1
    assert outcome.results[0].polished_genes == 2
    assert outcome.suppressed_unpolished == 0


def test_one_modelled_gene_is_not_enough(tmp_path):
    outcome = _run(tmp_path, [_annotated("mfa1", 100, 200), _raw("pra1", 300, 400)])
    assert outcome.results == []
    assert outcome.suppressed_unpolished == 1


def test_the_bar_can_be_lowered_and_then_the_locus_is_low_and_partial(tmp_path):
    """Rules 1 and 2, observed on a locus the bar would otherwise withhold."""
    outcome = _run(tmp_path, [_raw("mfa1", 100, 200), _raw("pra1", 300, 400)],
                   min_polished_genes=0)
    assert len(outcome.results) == 1
    r = outcome.results[0]
    assert r.polished_genes == 0
    assert r.confidence == "low"
    assert r.locus_class != LOCUS_CLASS_MAT


def test_repeated_hits_of_one_gene_do_not_add_up_to_the_bar(tmp_path):
    """Three HSP fragments of ONE alpha-box are one gene, not three."""
    outcome = _run(tmp_path, [
        _annotated("mfa1", 100, 200), _annotated("mfa1", 210, 300),
        _annotated("mfa1", 310, 400),
    ])
    assert outcome.results == [] or outcome.results[0].polished_genes == 1
