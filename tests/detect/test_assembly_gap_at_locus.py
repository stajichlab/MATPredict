"""An assembly gap where the locus should be is reported, not silently missed.

Curator's ruling 2026-09-26 (results/2026-09-26_cauris_uncalled/NOTE.md). Ten
uncalled C. auris genomes are reference-consensus assemblies with 8.1-8.9 kb
of N runs exactly over the reference idiomorph: the locus is not absent, the
assembly never resolved it. When a family's flank gene anchors a position
and a long N block sits there, the report says `assembly_gap_at_locus`.

The real gaps are not one run but a block of runs split by short islands of
sequence (16-143 bp), so runs are merged before the length is judged.
"""
import yaml

from MATPredict.detect.assembly_gap import (
    ANCHOR_MIN_IDENTITY, ANCHOR_WINDOW_BP, GAP_MERGE_BP, MIN_GAP_N_BP,
    find_gaps_at_locus, n_blocks,
)
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import run_pipeline
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.rollout_aggregate import aggregate_reports
from MATPredict.detect.search import SearchHit

from tests.detect.test_pipeline import _no_polish, _write_order, _write_record

KEY = FamilyKey("P", "aLocus")


def _hit(gene, role, start, end, contig="c1", identity=91.5):
    return SearchHit(KEY, gene, role, contig, start, end, "+", identity, "rec1",
                     "tblastn_genome", coverage=10.0)


def test_the_thresholds_are_stated():
    assert MIN_GAP_N_BP == 100
    assert GAP_MERGE_BP == 500
    assert ANCHOR_WINDOW_BP == 5_000
    assert ANCHOR_MIN_IDENTITY == 50.0


def test_n_blocks_are_one_based_and_inclusive():
    seq = "A" * 10 + "N" * 200 + "A" * 10
    assert n_blocks(seq) == [(11, 210, 200)]


def test_short_islands_between_runs_merge_into_one_block():
    """The C. auris SCO shape: runs split by 16-56 bp of sequence."""
    seq = "A" * 50 + "N" * 1270 + "A" * 42 + "N" * 1318 + "a" * 56 + "N" * 4521 + "A" * 50
    [(start, end, n)] = n_blocks(seq)
    assert (start, n) == (51, 1270 + 1318 + 4521)
    assert end == 50 + 1270 + 42 + 1318 + 56 + 4521


def test_runs_further_apart_than_the_merge_distance_stay_separate():
    seq = "N" * 150 + "A" * (GAP_MERGE_BP + 1) + "N" * 150
    assert len(n_blocks(seq)) == 2


def test_a_block_below_the_minimum_is_ignored():
    assert n_blocks("A" * 10 + "N" * (MIN_GAP_N_BP - 1) + "A" * 10) == []
    assert n_blocks("acgt" + "n" * MIN_GAP_N_BP + "acgt") == [(5, 104, 100)]


def _genome(tmp_path, seq, name="c1"):
    path = tmp_path / "genome.fa"
    path.write_text(f">{name} desc\n" + "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60)) + "\n")
    return path


def test_a_flank_hit_inside_the_gap_block_anchors_it(tmp_path):
    """The England C. auris shape: a 141 bp PIK1 island inside the N block."""
    seq = "A" * 1000 + "N" * 1569 + "A" * 21 + "N" * 3848 + "A" * 141 + "N" * 3357 + "A" * 1000
    genome = _genome(tmp_path, seq)
    island = 1000 + 1569 + 21 + 3848 + 1
    [gap] = find_gaps_at_locus({KEY: [_hit("PIK1", "flanking_variable", island, island + 140)]},
                               genome)
    assert gap.family_key == KEY and gap.contig == "c1"
    assert gap.start == 1001 and gap.n_bases == 1569 + 3848 + 3357
    assert gap.anchors == ["PIK1"]


def test_a_flank_hit_beside_the_block_anchors_it_within_the_window(tmp_path):
    seq = "A" * 10_000 + "N" * 500 + "A" * 20_000
    genome = _genome(tmp_path, seq)
    near = _hit("SLA2", "flanking_conserved", 10_500 + ANCHOR_WINDOW_BP, 10_600 + ANCHOR_WINDOW_BP)
    far = _hit("SLA2", "flanking_conserved", 10_502 + ANCHOR_WINDOW_BP, 10_700 + ANCHOR_WINDOW_BP)
    assert len(find_gaps_at_locus({KEY: [near]}, genome)) == 1
    assert find_gaps_at_locus({KEY: [far]}, genome) == []


def test_a_paralog_level_flank_hit_does_not_anchor(tmp_path):
    """Measured on 16 C. auris genomes: every false anchor was a flank paralog
    at 22.9-35.7% beside a 117-353 bp N block; the two real anchors were
    islands inside the MTL gap at 91.5% and 100%."""
    genome = _genome(tmp_path, "A" * 1000 + "N" * 133 + "A" * 1000)
    weak = _hit("PAP1", "flanking_variable", 900, 1000, identity=ANCHOR_MIN_IDENTITY - 0.1)
    strong = _hit("PAP1", "flanking_variable", 900, 1000, identity=ANCHOR_MIN_IDENTITY)
    assert find_gaps_at_locus({KEY: [weak]}, genome) == []
    assert len(find_gaps_at_locus({KEY: [strong]}, genome)) == 1


def test_a_core_hit_does_not_anchor(tmp_path):
    """The ruling anchors on flanks; a lone core fragment is not a position."""
    seq = "A" * 1000 + "N" * 500 + "A" * 1000
    genome = _genome(tmp_path, seq)
    assert find_gaps_at_locus({KEY: [_hit("mfa1", "core_MAT", 900, 1000)]}, genome) == []


def test_a_flank_hit_on_another_contig_does_not_anchor(tmp_path):
    genome = _genome(tmp_path, "A" * 1000 + "N" * 500 + "A" * 1000)
    assert find_gaps_at_locus(
        {KEY: [_hit("PIK1", "flanking_variable", 900, 1000, contig="c2")]}, genome) == []


def test_an_unreadable_genome_reports_no_gap(tmp_path):
    assert find_gaps_at_locus({KEY: [_hit("PIK1", "flanking_variable", 1, 90)]},
                              tmp_path / "missing.fa") == []


ORDER = (
    "phylum: P\nloci:\n  - locus_name: aLocus\n    vocabulary_type: pattern\n"
    "    idiomorph_pattern: \"^a[0-9]+$\"\n    taxonomic_scope: [1]\n"
    "    genes:\n      - {name: mfa1, role: core_MAT}\n"
    "      - {name: flk1, role: flanking_conserved}\n"
)


def _run(tmp_path, hits, seq):
    _write_order(tmp_path, ORDER)
    _write_record(tmp_path)
    genome = _genome(tmp_path, seq)
    return run_pipeline(
        genome_fasta=genome, proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: [], search_localize=lambda *a, **k: hits,
        polish_with_exonerate=_no_polish, polish_with_miniprot=_no_polish,
    )


def test_the_pipeline_reports_a_gap_for_a_family_with_no_call(tmp_path):
    seq = "A" * 2000 + "N" * 3000 + "A" * 2000
    outcome = _run(tmp_path, [_hit("flk1", "flanking_conserved", 1500, 1700)], seq)
    assert outcome.results == []
    [gap] = outcome.assembly_gaps_at_locus
    assert (gap.start, gap.end, gap.n_bases) == (2001, 5000, 3000)

    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    assert doc["assembly_gap_at_locus"] == [{
        "family": "P:aLocus", "contig": "c1", "start": 2001, "end": 5000,
        "n_bases": 3000, "anchors": ["flk1"],
    }]


def test_the_report_writes_an_empty_list_when_there_is_no_gap(tmp_path):
    outcome = _run(tmp_path, [_hit("flk1", "flanking_conserved", 1500, 1700)], "A" * 7000)
    assert outcome.assembly_gaps_at_locus == []
    path = tmp_path / "report.yaml"
    write_detection_report(outcome, path)
    assert yaml.safe_load(path.read_text())["assembly_gap_at_locus"] == []


def test_the_rollout_counts_gap_genomes_separately(tmp_path):
    base = {"routing_mode": "lineage", "families_attempted": ["P:aLocus"], "detected": [],
            "not_detected": [{"family": "P:aLocus", "reason": "x"}]}
    gap = dict(base, assembly_gap_at_locus=[{"family": "P:aLocus", "contig": "c1", "start": 1,
                                             "end": 200, "n_bases": 200, "anchors": ["flk1"]}])
    paths = []
    for name, doc in (("1_gapped", gap), ("1_plain", base)):
        (tmp_path / name).mkdir()
        p = tmp_path / name / "detection_report.yaml"
        p.write_text(yaml.safe_dump(doc))
        paths.append(p)
    summary = aggregate_reports(paths, lineage_resolver=lambda taxid: None)
    assert summary.assembly_gap_at_locus == ["1_gapped"]
    assert summary.to_doc()["assembly_gap_at_locus"] == ["1_gapped"]
