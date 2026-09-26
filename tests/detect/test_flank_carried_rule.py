"""A call whose core genes were never modelled rests on its flanks alone.

Curator's ruling 2026-09-26 (docs/notes/2026-09-26_cauris-flank-carried-and-
calbicans-zygosity.md; measured in docs/notes/2026-09-26_polish-cap-measured-
and-serinales-scan.md), revised the same day after the Ascomycota audit
(results/2026-09-26_flank_rule_ascomycota/NOTE.md). The rule judges the
STRONGEST core hit (lowest e-value):

1. E <= 1e-5 and within the family's `flank_carried_window_bp` of the flank
   span: keep, cap at `low`, class `partial_locus`, flag
   `idiomorph_unmodelled: true`;
2. otherwise: withhold, like any bar failure.

Generic over the family's own flanking genes (role `flanking_*`), not
hardcoded to PAP1/OBP1/PIK1.
"""
import yaml

from MATPredict.detect.family_registry import FamilyKey, load_all_families
from MATPredict.detect.flank_carried import (
    DEFAULT_FLANK_CARRIED_WINDOW_BP, FLANK_CARRIED_MAX_EVALUE, WITHHELD_FLANK_CARRIED,
    apply_flank_carried_rule,
)
from MATPredict.detect.idiomorph import LOCUS_CLASS_MAT, LOCUS_CLASS_PARTIAL
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence, run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _no_polish, _write_order, _write_record

KEY = FamilyKey("P", "aLocus")
STRONG, WEAK, NOISE = 1e-13, 1e-3, 2.0


def _ev(gene, role, start, end, status, method="tblastn_genome", contig="c1", evalue=None):
    return GeneEvidence(gene, role, contig, start, end, "+", 40.0, 50.0, "rec1",
                        method, status=status, evalue=evalue)


def _core(start, end, evalue=STRONG, contig="c1", gene="mfa1"):
    return _ev(gene, "core_MAT", start, end, "unpolished", contig=contig, evalue=evalue)


def _result(evidence, confidence="medium", locus_class=LOCUS_CLASS_MAT, key=KEY):
    return DetectionResult(
        family_key=key, contig="c1", start=min(e.start for e in evidence),
        end=max(e.end for e in evidence), confidence=confidence, idiomorph="a1",
        ambiguous_with=[], genes_found=sorted({e.gene_name for e in evidence}),
        genes_missing=[], fragmented=False, gene_evidence=evidence,
        locus_class=locus_class, polished_genes=2,
    )


FLANKS = [
    _ev("flk1", "flanking_conserved", 1_000, 2_000, "polished_single"),
    _ev("flk2", "flanking_variable", 10_000, 11_000, "polished_agree"),
]


def test_the_defaults_are_the_rulings_numbers():
    assert DEFAULT_FLANK_CARRIED_WINDOW_BP == 3_000
    assert FLANK_CARRIED_MAX_EVALUE == 1e-5


def test_a_call_with_a_modelled_core_gene_is_untouched():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 50_000, 51_000, "polished_single")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [r] and withheld == []


def test_an_annotated_core_gene_counts_as_modelled():
    r = _result(FLANKS + [_ev("mfa1", "core_MAT", 50_000, 51_000, "not_polish_candidate",
                              method="diamond_proteome")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [r] and withheld == []


def test_a_strong_core_hit_inside_the_flank_span_is_kept_at_low():
    r = _result(FLANKS + [_core(5_000, 5_500)], confidence="high")
    [kept], withheld = apply_flank_carried_rule([r])
    assert withheld == []
    assert kept.confidence == "low"
    assert kept.locus_class == LOCUS_CLASS_PARTIAL
    assert kept.idiomorph_unmodelled is True


def test_the_default_window_extends_the_span_on_both_sides():
    just_inside = _result(FLANKS + [_core(13_500, 14_000)])
    just_outside = _result(FLANKS + [_core(14_001, 14_300)])
    kept, withheld = apply_flank_carried_rule([just_inside, just_outside])
    assert [k.end for k in kept] == [14_000]
    assert [w.end for w in withheld] == [14_300]


def test_a_family_window_widens_the_span():
    """SLA2/APN2/COX13 sit outside the idiomorph: the core can be ~10 kb off."""
    r = _result(FLANKS + [_core(20_000, 20_500)])
    assert apply_flank_carried_rule([r])[0] == []
    [kept], withheld = apply_flank_carried_rule([r], {KEY: 20_000})
    assert withheld == [] and kept.confidence == "low"


def test_a_family_window_can_also_narrow_it():
    r = _result(FLANKS + [_core(12_000, 12_500)])
    kept, withheld = apply_flank_carried_rule([r], {KEY: 500})
    assert kept == [] and len(withheld) == 1


def test_a_strong_core_hit_outside_the_window_is_withheld():
    """The Debaryomyces case: the core fragment sits 8 kb beyond the block."""
    r = _result(FLANKS + [_core(19_000, 19_400)])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == []
    [w] = withheld
    assert w.withheld_reason == WITHHELD_FLANK_CARRIED


def test_a_weak_or_noise_strongest_hit_is_withheld_even_inside():
    for evalue in (WEAK, NOISE, None):
        r = _result(FLANKS + [_core(5_000, 5_500, evalue=evalue)])
        kept, withheld = apply_flank_carried_rule([r])
        assert kept == [] and len(withheld) == 1, evalue


def test_the_evalue_floor_is_inclusive():
    r = _result(FLANKS + [_core(5_000, 5_500, evalue=1e-5)])
    [kept], _ = apply_flank_carried_rule([r])
    assert kept.idiomorph_unmodelled is True


def test_a_stray_weak_hit_far_away_no_longer_withholds_a_good_call():
    """Didymobotryum rigidum: MAT genes between SLA2 and APN2, a stray MAT1-2-4
    hit 60 kb off. Only the strongest core hit is judged."""
    r = _result(FLANKS + [_core(5_000, 5_500, evalue=STRONG),
                          _core(70_000, 70_500, evalue=NOISE, gene="pra1")])
    [kept], withheld = apply_flank_carried_rule([r])
    assert withheld == [] and kept.confidence == "low"


def test_when_the_strongest_hit_is_far_the_call_is_withheld():
    """A nearer but weaker hit does not rescue it."""
    r = _result(FLANKS + [_core(5_000, 5_500, evalue=WEAK),
                          _core(70_000, 70_500, evalue=STRONG, gene="pra1")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [] and len(withheld) == 1


def test_a_strongest_hit_on_another_contig_is_outside():
    r = _result(FLANKS + [_core(5_000, 5_500, contig="c2")])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [] and len(withheld) == 1


def test_a_call_with_no_core_hit_at_all_is_withheld():
    """Nothing places a MAT gene between the flanks; there is no call to make."""
    kept, withheld = apply_flank_carried_rule([_result(FLANKS)])
    assert kept == [] and len(withheld) == 1


def test_a_call_with_no_modelled_flank_is_left_to_the_other_bars():
    """Not flank-carried: the modelled-gene bar and tiering already rule on it."""
    r = _result([_ev("flk1", "flanking_conserved", 1_000, 2_000, "unpolished"),
                 _core(5_000, 5_500)])
    kept, withheld = apply_flank_carried_rule([r])
    assert kept == [r] and withheld == []


def test_the_curated_windows(tmp_path):
    """Serinales flanks sit inside the idiomorph (3 kb); SLA2/APN2/COX13 and the
    Mucorales tptA/rnhA flanks sit outside it (20 kb)."""
    from pathlib import Path

    db = Path(__file__).resolve().parents[2] / "db"
    windows = {f.key: f.flank_carried_window_bp for f in load_all_families(db)}
    assert windows[FamilyKey("Ascomycota", "MTL")] == 3_000
    for locus in ("MAT", "MATtub", "MATyl", "MATsc"):
        assert windows[FamilyKey("Ascomycota", locus)] == 20_000, locus
    assert windows[FamilyKey("Mucoromycota", "MAT")] == 20_000


ORDER = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
    "      - {name: flk2, role: flanking_conserved}\n"
)


def _annotated(gene, role, start, end):
    return SearchHit(KEY, gene, role, "c1", start, end, "+", 91.5, "rec1",
                     "diamond_proteome", coverage=88.0)


def _raw(gene, start, end, evalue):
    return SearchHit(KEY, gene, "core_MAT", "c1", start, end, "+", 31.0, "rec1",
                     "tblastn_genome", coverage=12.0, evalue=evalue)


def _run(tmp_path, core_start, core_end, evalue=STRONG, extra=()):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    hits = [_annotated("flk1", "flanking_conserved", 100, 200),
            _annotated("flk2", "flanking_conserved", 5_000, 5_100),
            _raw("mfa1", core_start, core_end, evalue), *extra]
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: hits, search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )


def test_the_pipeline_keeps_a_strong_inside_call_at_low(tmp_path):
    outcome = _run(tmp_path, 2_000, 2_300)
    [r] = outcome.results
    assert r.confidence == "low"
    assert r.locus_class == LOCUS_CLASS_PARTIAL
    assert r.idiomorph_unmodelled is True


def test_the_evidence_carries_the_genes_best_evalue(tmp_path):
    """The report's hit for a gene is chosen by method and identity; its
    e-value is the gene's best in the cluster, as the audit measured it."""
    outcome = _run(tmp_path, 2_000, 2_300, evalue=0.5,
                   extra=[_raw("mfa1", 2_400, 2_500, 1e-9)])
    [r] = outcome.results
    [core] = [e for e in r.gene_evidence if e.role == "core_MAT"]
    assert core.evalue == 1e-9


def test_the_pipeline_withholds_a_noise_call_and_says_why(tmp_path):
    outcome = _run(tmp_path, 2_000, 2_300, evalue=NOISE)
    assert outcome.results == []
    [w] = outcome.suppressed_loci
    assert w.withheld_reason == WITHHELD_FLANK_CARRIED
    assert outcome.suppressed_flank_carried == 1
    [reason] = [n.reason for n in outcome.not_detected if n.family_key == KEY]
    assert "flank" in reason

    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["suppressed_flank_carried"] == 1
    assert doc["suppressed_loci"][0]["withheld_reason"] == WITHHELD_FLANK_CARRIED


def test_the_pipeline_withholds_an_outside_call(tmp_path):
    outcome = _run(tmp_path, 9_000, 9_300)
    assert outcome.results == [] and outcome.suppressed_flank_carried == 1


def test_the_report_writes_the_flag_and_the_evalue(tmp_path):
    outcome = _run(tmp_path, 2_000, 2_300)
    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["detected"][0]["idiomorph_unmodelled"] is True
    [core] = [e for e in doc["detected"][0]["gene_evidence"] if e["role"] == "core_MAT"]
    assert core["evalue"] == STRONG
