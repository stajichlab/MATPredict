"""A MAT locus split across contigs by the assembly is reported, not withheld.

Curator's ruling 2026-09-27 (results/2026-09-27_rarrhizus_uncalled/NOTE.md):
a single MODELLED core gene at >= 95% identity near a contig end, whose
family's flanking genes are found on OTHER contigs at similar identity, is
reported as `partial_locus`, confidence `low`, flagged `split_locus`, instead
of being withheld by the two-modelled-gene bar.

The motivating case: 12 Rhizopus arrhizus/delemar GL-series assemblies carry
sexP at 98.4% alone on a contig 1-199 bp from its end, tptA/btbA on a second
contig and rnhA on a third.
"""
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult, GeneEvidence
from MATPredict.detect.search import SearchHit
from MATPredict.detect.split_locus import (
    DEFAULT_SPLIT_LOCUS, SplitLocusParams, evaluate_split_locus, split_locus_params,
)

KEY = FamilyKey("Mucoromycota", "MAT")
ROLES = {"sexP": "core_MAT", "sexM": "core_MAT", "tptA": "flanking_conserved",
         "rnhA": "flanking_conserved", "btbA": "flanking_variable",
         "glrA": "flanking_variable"}
LENGTHS = {"cA": 3_932, "cB": 10_802, "cC": 8_177, "cD": 50_000}


def _hit(gene, contig, start, end, identity, bitscore, key=KEY):
    return SearchHit(key, gene, ROLES[gene], contig, start, end, "+", identity,
                     "rec", "tblastn_genome", bitscore=bitscore)


def _ev(gene, contig, start, end, identity, status="polished_agree"):
    return GeneEvidence(gene, ROLES[gene], contig, start, end, "+", identity, 100.0,
                        "rec", "miniprot_refine", status=status)


def _candidate(core_contig="cA", core_start=2_803, core_end=3_741, identity=98.4,
               status="polished_agree", gene="sexP"):
    ev = [_ev(gene, core_contig, core_start, core_end, identity, status)]
    return DetectionResult(
        family_key=KEY, contig=core_contig, start=core_start, end=core_end,
        confidence="low", idiomorph="Plus", ambiguous_with=[], genes_found=[gene],
        genes_missing=["tptA", "rnhA"], fragmented=False, gene_evidence=ev,
        polished_genes=1,
    )


def _genome_hits(**over):
    hits = [
        _hit("sexP", "cA", 2_803, 3_741, 98.4, 520.0),
        _hit("tptA", "cB", 5_320, 6_649, 99.4, 900.0),
        _hit("btbA", "cB", 7_540, 10_680, 98.6, 2_000.0),
        _hit("rnhA", "cC", 333, 4_219, 100.0, 1_500.0),
        # an HMG paralog elsewhere: weaker, must not matter
        _hit("sexP", "cD", 14_929, 15_300, 40.5, 60.0),
    ]
    return hits


def _eval(candidate=None, hits=None, reported=(), params=DEFAULT_SPLIT_LOCUS):
    return evaluate_split_locus(
        candidate if candidate is not None else _candidate(),
        hits if hits is not None else _genome_hits(),
        family_roles=ROLES, contig_lengths=LENGTHS, reported_spans=list(reported),
        params=params,
    )


def test_defaults_are_the_rulings_numbers():
    assert DEFAULT_SPLIT_LOCUS.enabled is True
    assert DEFAULT_SPLIT_LOCUS.min_core_identity == 95.0
    assert DEFAULT_SPLIT_LOCUS.max_edge_bp == 500
    assert DEFAULT_SPLIT_LOCUS.min_flank_genes == 2
    assert DEFAULT_SPLIT_LOCUS.flank_identity_gap == 10.0


def test_the_rhizopus_gl_case_is_accepted():
    info = _eval()
    assert info is not None
    assert info["core_gene"] == "sexP" and info["core_contig"] == "cA"
    assert info["edge_distance"] == 191
    assert {f["gene"] for f in info["flanks"]} == {"tptA", "btbA", "rnhA"}
    assert set(info["contigs"]) == {"cA", "cB", "cC"}


def test_an_unmodelled_core_gene_is_not_enough():
    assert _eval(_candidate(status="unpolished")) is None


def test_an_annotated_fast_path_core_gene_counts_as_modelled():
    c = _candidate(status="not_polish_candidate")
    c = DetectionResult(**{**c.__dict__, "gene_evidence": [
        GeneEvidence("sexP", "core_MAT", "cA", 2_803, 3_741, "+", 98.4, 100.0, "rec",
                     "diamond_proteome", status="not_polish_candidate")]})
    assert _eval(c) is not None


def test_a_core_gene_below_95_percent_is_rejected():
    assert _eval(_candidate(identity=94.9)) is None


def test_a_core_gene_far_from_a_contig_end_is_rejected():
    # 20 kb into a 50 kb contig: the locus is not broken by the assembly here
    assert _eval(_candidate(core_contig="cD", core_start=20_000, core_end=21_000)) is None


def test_the_left_contig_end_counts_too():
    c = _candidate(core_start=101, core_end=1_000)
    hits = [h for h in _genome_hits() if not (h.gene_name == "sexP" and h.contig == "cA")]
    hits.append(_hit("sexP", "cA", 101, 1_000, 98.4, 520.0))
    info = _eval(c, hits)
    assert info is not None and info["edge_distance"] == 100


def test_flanks_must_reach_similar_identity():
    # flanks 30 points below the core: they may be paralogs, not this locus
    hits = [h if h.role == "core_MAT" else
            SearchHit(**{**h.__dict__, "identity": h.identity - 30.0})
            for h in _genome_hits()]
    assert _eval(hits=hits) is None


def test_at_least_two_distinct_flank_genes_are_required():
    hits = [h for h in _genome_hits() if h.gene_name not in ("btbA", "rnhA")]
    assert _eval(hits=hits) is None


def test_a_flank_must_lie_on_another_contig():
    # every flank on the core's own contig is an ordinary locus, not a split
    hits = [_hit("sexP", "cA", 2_803, 3_741, 98.4, 520.0),
            _hit("tptA", "cA", 100, 1_400, 99.4, 900.0),
            _hit("rnhA", "cA", 1_500, 2_700, 100.0, 1_500.0)]
    assert _eval(hits=hits) is None


def test_the_core_hit_must_be_the_genomes_best_for_its_gene():
    # a stronger sexP elsewhere means this one is the paralog
    hits = _genome_hits() + [_hit("sexP", "cD", 30_000, 31_000, 99.9, 900.0)]
    assert _eval(hits=hits) is None


def test_a_flank_inside_a_reported_locus_does_not_count():
    # tptA and btbA belong to a locus already reported elsewhere
    reported = [("cB", 1, 10_802)]
    assert _eval(reported=reported) is None


def test_the_flanks_best_hit_is_used_not_any_hit():
    # rnhA's best hit is on cC; a weak rnhA hit elsewhere changes nothing
    hits = _genome_hits() + [_hit("rnhA", "cD", 40_000, 41_000, 35.0, 50.0)]
    info = _eval(hits=hits)
    assert info is not None
    assert next(f for f in info["flanks"] if f["gene"] == "rnhA")["contig"] == "cC"


def test_a_family_can_turn_it_off():
    assert _eval(params=SplitLocusParams(enabled=False)) is None


def test_roster_overrides_are_read():
    p = split_locus_params({"split_locus": {"max_edge_bp": 1_000, "min_flank_genes": 1}})
    assert p.max_edge_bp == 1_000 and p.min_flank_genes == 1
    assert p.min_core_identity == 95.0
    assert split_locus_params({}) == DEFAULT_SPLIT_LOCUS
    assert split_locus_params({"split_locus": False}).enabled is False


def test_modelled_statuses_match_polish():
    from MATPredict.detect.polish import STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE
    from MATPredict.detect.split_locus import MODELLED_STATUSES
    assert MODELLED_STATUSES == {STATUS_AGREE, STATUS_DISAGREE, STATUS_SINGLE}


def test_every_curated_family_gets_the_default_rule():
    from pathlib import Path
    from MATPredict.detect.family_registry import load_all_families
    db = Path(__file__).resolve().parents[2] / "db"
    assert {f.split_locus for f in load_all_families(db)} == {DEFAULT_SPLIT_LOCUS}


# ---- pipeline level -------------------------------------------------------
import yaml as _yaml

from MATPredict.detect.idiomorph import LOCUS_CLASS_PARTIAL
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.report import write_detection_report

from tests.detect.test_pipeline import _no_polish, _write_order, _write_record

PKEY = FamilyKey("P", "aLocus")
PORDER = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
    "      - {name: flk2, role: flanking_conserved}\n"
    "      - {name: flk3, role: flanking_variable}\n"
)


def _ann(gene, role, contig, start, end, identity=98.5):
    return SearchHit(PKEY, gene, role, contig, start, end, "+", identity, "rec1",
                     "diamond_proteome", coverage=99.0, bitscore=500.0)


def _genome(tmp_path, lengths):
    with open(tmp_path / "genome.fa", "w") as fh:
        for name, n in lengths.items():
            fh.write(f">{name}\n{'A' * n}\n")


def _prun(tmp_path, hits, lengths, order=PORDER):
    _write_order(tmp_path, order)
    _write_record(tmp_path)
    _genome(tmp_path, lengths)
    return run_pipeline(
        genome_fasta=tmp_path / "genome.fa", proteome_fasta=tmp_path / "proteome.faa",
        taxid=None, db_root=tmp_path, reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: hits, search_localize=lambda *a, **k: [],
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )


SPLIT_LENGTHS = {"cA": 4_000, "cB": 12_000, "cC": 9_000}
SPLIT_HITS = [
    _ann("mfa1", "core_MAT", "cA", 2_800, 3_800),        # 200 bp from the end
    _ann("flk1", "flanking_conserved", "cB", 5_000, 6_000),
    _ann("flk3", "flanking_variable", "cB", 7_000, 9_000),
    _ann("flk2", "flanking_conserved", "cC", 300, 4_000),
]


def test_the_pipeline_reports_a_split_locus_at_low(tmp_path):
    outcome = _prun(tmp_path, SPLIT_HITS, SPLIT_LENGTHS)
    [r] = outcome.results
    assert r.family_key == PKEY
    assert r.confidence == "low"
    assert r.locus_class == LOCUS_CLASS_PARTIAL
    assert r.split_locus["core_gene"] == "mfa1"
    assert r.split_locus["contigs"] == ["cA", "cB", "cC"]
    assert not [n for n in outcome.not_detected if n.family_key == PKEY]

    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = _yaml.safe_load(path.read_text())
    assert doc["detected"][0]["split_locus"]["edge_distance"] == 200


def test_the_pipeline_leaves_a_core_gene_deep_in_a_contig_withheld(tmp_path):
    lengths = {"cA": 40_000, "cB": 12_000, "cC": 9_000}
    hits = [_ann("mfa1", "core_MAT", "cA", 20_000, 21_000)] + SPLIT_HITS[1:]
    outcome = _prun(tmp_path, hits, lengths)
    assert outcome.results == []


def test_the_roster_can_turn_it_off(tmp_path):
    order = PORDER.replace("taxonomic_scope: [1]\n", "taxonomic_scope: [1]\n    split_locus: false\n")
    outcome = _prun(tmp_path, SPLIT_HITS, SPLIT_LENGTHS, order=order)
    assert outcome.results == []


def test_an_ordinary_called_locus_is_not_touched(tmp_path):
    # core and two flanks on one contig: the normal path calls it; no split flag
    hits = [_ann("mfa1", "core_MAT", "cA", 1_000, 2_000),
            _ann("flk1", "flanking_conserved", "cA", 3_000, 4_000),
            _ann("flk2", "flanking_conserved", "cA", 5_000, 6_000)]
    outcome = _prun(tmp_path, hits, {"cA": 50_000})
    [r] = outcome.results
    assert r.split_locus is None
