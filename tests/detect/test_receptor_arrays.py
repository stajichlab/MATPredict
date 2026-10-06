"""Report-only pheromone-receptor arrays and the `array_support` flag.

Study: analysis/2026-10-06_agaricomycetes-pr-arrays.md, options 1 and 2. STE3-like
loci on one contig within 50 kb are one array; an array is `supported` when it
holds >= 2 loci, or a precursor-homology hit, or >= 2 distinct strict-CAAX ORFs.
The flag never changes a call, its confidence, its verification or the counts,
and array membership does not make a locus a mating receptor.
"""
from dataclasses import replace

import yaml

from MATPredict.detect import pipeline
from MATPredict.detect.caax import CAAX_METHOD
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionOutcome, DetectionResult, run_pipeline
from MATPredict.detect.receptor_arrays import (
    ARRAY_GAP_BP, LOCI_ARRAY_COLUMNS, RECEPTOR_ARRAYS_NOTE, SUPPORTED, UNSUPPORTED,
    attach_array_support, build_receptor_arrays, group_arrays, loci_columns, merge_receptor_loci,
)
from MATPredict.detect.report import write_detection_report
from MATPredict.detect.search import SearchHit

from tests.detect.test_caax_precursor import _genome, _order, _orf, _receptor
from tests.detect.test_pipeline import _model

KEY = FamilyKey("P", "PR")
OTHER = FamilyKey("P", "HD")


class _Fam:
    key = KEY
    pheromone_precursor_scan = {
        "gene": "caax_precursor", "receptor_genes": ["pheromone_receptor"],
        "motif": "C[VI][IV][AVMG]", "window_bp": 10_000, "min_codons": 20, "max_codons": 130,
    }
    genes = [
        {"name": "pheromone_receptor", "role": "core_MAT", "gene_class": "pheromone_receptor"},
        {"name": "pheromone_B43", "role": "core_MAT", "gene_class": "pheromone_precursor"},
        {"name": "caax_precursor", "role": "core_MAT", "gene_class": "pheromone_precursor",
         "optional": True},
    ]


class _NoScan:
    key = OTHER
    pheromone_precursor_scan = None
    genes = []


def rec(start, end, strand="+", contig="c1", ref="rec1"):
    return SearchHit(KEY, "pheromone_receptor", "core_MAT", contig, start, end, strand, 55.0,
                     ref, "tblastn_genome", coverage=60.0)


def caax(start, end, strand="+", contig="c1"):
    return SearchHit(KEY, "caax_precursor", "core_MAT", contig, start, end, strand, 0.0,
                     "caax_scan:CVIA", CAAX_METHOD, align_length_aa=30)


def prec(start, end, contig="c1"):
    return SearchHit(KEY, "pheromone_B43", "core_MAT", contig, start, end, "+", 40.0,
                     "rec1", "tblastn_genome")


def arrays(hits):
    return build_receptor_arrays([_Fam()], hits)


# ---- grouping -----------------------------------------------------------

def test_the_gap_is_the_studys_50_kb_between_locus_ends():
    assert ARRAY_GAP_BP == 50_000
    a = rec(1_000, 2_000)
    joined = rec(2_000 + ARRAY_GAP_BP, 2_000 + ARRAY_GAP_BP + 900)
    apart = rec(2_000 + ARRAY_GAP_BP + 1, 2_000 + ARRAY_GAP_BP + 900)
    assert [x.size for x in arrays([a, joined])] == [2]
    assert sorted(x.size for x in arrays([a, apart])) == [1, 1]


def test_single_linkage_chains_through_the_furthest_end():
    # a long locus keeps the array open: the gap is measured from the furthest end
    hits = [rec(1_000, 80_000), rec(20_000, 21_000, "-"), rec(120_000, 121_000)]
    [arr] = arrays(hits)
    assert arr.size == 3 and (arr.start, arr.end) == (1_000, 121_000)


def test_fragments_of_one_gene_are_one_locus_not_an_array():
    # an intron-split HSP pair and a second reference over the same gene
    hits = [rec(1_000, 1_500), rec(1_650, 2_400), rec(1_020, 2_390, ref="rec2")]
    [arr] = arrays(hits)
    assert arr.size == 1 and arr.members == ((1_000, 2_400, "+"),)


def test_an_array_holds_loci_of_both_strands_and_opposite_strands_are_not_merged():
    hits = [rec(1_000, 2_500, "+"), rec(1_200, 2_600, "-"), rec(9_000, 10_500, "-")]
    [arr] = arrays(hits)
    assert arr.size == 3
    assert [m[2] for m in arr.members] == ["+", "-", "-"]
    assert arr.member_strings()[0] == "1000-2500:+"


def test_arrays_never_span_a_contig_boundary():
    hits = [rec(1_000, 2_000, contig="c1"), rec(1_500, 2_500, contig="c2"),
            rec(3_000, 4_000, contig="c1")]
    out = arrays(hits)
    assert sorted((a.contig, a.size) for a in out) == [("c1", 2), ("c2", 1)]
    # an array near a contig end is still bounded by that contig's own loci
    [c2] = [a for a in out if a.contig == "c2"]
    assert (c2.start, c2.end) == (1_500, 2_500)


def test_helpers_group_and_merge_directly():
    loci = merge_receptor_loci([rec(10, 20), rec(30, 40, "-")])
    assert loci == [("c1", 10, 20, "+"), ("c1", 30, 40, "-")]
    assert len(group_arrays(loci)) == 1


def test_array_id_is_stable_and_names_family_contig_and_span():
    [arr] = arrays([rec(5_000, 6_000), rec(9_000, 9_900)])
    assert arr.array_id == "P:PR:c1:5000-9900"
    assert arrays([rec(9_000, 9_900), rec(5_000, 6_000)])[0].array_id == arr.array_id


# ---- support ------------------------------------------------------------

def test_supported_by_array_size():
    [arr] = arrays([rec(1_000, 2_000), rec(20_000, 21_000)])
    assert arr.support == SUPPORTED and arr.reasons == ("array_size>=2",)


def test_supported_by_precursor_homology():
    [arr] = arrays([rec(1_000, 2_000), prec(5_000, 5_150)])
    assert arr.size == 1 and arr.support == SUPPORTED
    assert arr.reasons == ("precursor_homology",) and arr.n_precursor_hits == 1


def test_precursor_homology_beyond_the_window_does_not_count():
    [arr] = arrays([rec(1_000, 2_000), prec(2_000 + 10_001, 12_300)])
    assert arr.support == UNSUPPORTED and arr.n_precursor_hits == 0


def test_supported_by_two_distinct_strict_caax_orfs():
    [arr] = arrays([rec(1_000, 2_000), caax(3_000, 3_090), caax(6_000, 6_090, "-")])
    assert arr.size == 1 and arr.support == SUPPORTED
    assert arr.reasons == ("caax_orfs>=2",) and arr.n_caax_orfs == 2


def test_one_caax_orf_is_not_enough():
    [arr] = arrays([rec(1_000, 2_000), caax(3_000, 3_090)])
    assert arr.support == UNSUPPORTED and arr.n_caax_orfs == 1
    assert "strict_caax_orfs=1(<2)" in arr.reasons


def test_the_same_orf_seen_twice_counts_once():
    [arr] = arrays([rec(1_000, 2_000), caax(3_000, 3_090), caax(3_090, 3_000)])
    assert arr.n_caax_orfs == 1 and arr.support == UNSUPPORTED


def test_unsupported_singleton_lists_every_reason():
    [arr] = arrays([rec(1_000, 2_000)])
    assert arr.support == UNSUPPORTED
    assert arr.reasons == ("single_locus", "no_precursor_homology", "strict_caax_orfs=0(<2)")


def test_caax_and_precursor_hits_of_another_contig_or_a_superseded_receptor_are_ignored():
    sup = replace(rec(50_000, 51_000), superseded_by="other")
    out = arrays([rec(1_000, 2_000), caax(3_000, 3_090, contig="c9"), caax(4_000, 4_090, contig="c9"),
                  prec(3_000, 3_100, contig="c9"), sup])
    [arr] = out
    assert arr.size == 1 and arr.support == UNSUPPORTED


def test_no_arrays_for_a_family_without_the_scan():
    assert build_receptor_arrays([_NoScan()], [rec(1, 100)]) == []


# ---- attaching to calls -------------------------------------------------

def _call(family=KEY, start=900, end=21_100, **kw):
    return DetectionResult(
        family_key=family, contig="c1", start=start, end=end, confidence="medium",
        idiomorph="undetermined", ambiguous_with=[], genes_found=["pheromone_receptor"],
        genes_missing=[], fragmented=False, **kw)


def test_a_call_gets_its_array_and_only_that_field_changes():
    out = arrays([rec(1_000, 2_000), rec(20_000, 21_000)])
    call = _call(verification={"status": "unverified", "reason": "x"}, caax_dependent=True)
    [new], [arr] = attach_array_support([call], [_Fam()], out)
    assert new.receptor_array["array_id"] == out[0].array_id
    assert new.receptor_array["array_size"] == 2
    assert new.receptor_array["array_members"] == ["1000-2000:+", "20000-21000:+"]
    assert new.receptor_array["array_support"] == SUPPORTED
    assert replace(new, receptor_array=None) == call
    assert arr.n_calls == 1


def test_two_calls_in_one_array_list_the_array_once_with_both_calls():
    out = arrays([rec(1_000, 2_000), rec(20_000, 21_000)])
    calls = [_call(start=900, end=2_100), _call(start=19_900, end=21_100)]
    new, arrs = attach_array_support(calls, [_Fam()], out)
    assert len(arrs) == 1 and arrs[0].n_calls == 2
    assert new[0].receptor_array["array_id"] == new[1].receptor_array["array_id"]


def test_a_call_that_is_not_pr_is_untouched_and_a_merged_pr_call_is_found():
    out = arrays([rec(1_000, 2_000)])
    hd = _call(family=OTHER)
    merged = _call(family=FamilyKey("P", "Bbeta"), merged_from=[{"family": "P:PR"}])
    new, _ = attach_array_support([hd, merged], [_Fam()], out)
    assert new[0].receptor_array is None
    assert new[1].receptor_array["array_size"] == 1


def test_a_pr_call_beside_no_array_reports_nulls_not_a_guess():
    other_contig = replace(_call(), contig="c9")
    [new], _ = attach_array_support([other_contig], [_Fam()], arrays([rec(1_000, 2_000)]))
    assert new.receptor_array["array_id"] is None and new.receptor_array["array_support"] is None


# ---- report and loci.tsv columns ---------------------------------------

def test_the_report_carries_call_fields_arrays_and_the_note(tmp_path):
    out = arrays([rec(1_000, 2_000), rec(20_000, 21_000)])
    new, arrs = attach_array_support([_call(), _call(family=OTHER)], [_Fam()], out)
    outcome = DetectionOutcome(results=new, receptor_arrays=arrs)
    path = tmp_path / "r.yaml"
    write_detection_report(outcome, path)
    doc = yaml.safe_load(path.read_text())
    pr, hd = doc["detected"]
    assert pr["array_id"] == "P:PR:c1:1000-21000" and pr["array_size"] == 2
    assert pr["array_support"] == SUPPORTED
    assert pr["array_support_reasons"] == ["array_size>=2"]
    assert not any(k.startswith("array_") for k in hd)
    [a] = doc["receptor_arrays"]
    assert a["calls"] == 1 and a["array_members"] == ["1000-2000:+", "20000-21000:+"]
    assert "does not establish" in doc["receptor_arrays_note"]
    assert "both mating and non-mating" in doc["receptor_arrays_note"]
    assert doc["receptor_arrays_note"] == RECEPTOR_ARRAYS_NOTE


def test_a_report_without_arrays_writes_an_empty_list(tmp_path):
    path = tmp_path / "r.yaml"
    write_detection_report(DetectionOutcome(results=[]), path)
    assert yaml.safe_load(path.read_text())["receptor_arrays"] == []


def test_loci_columns_for_a_pr_call_and_a_non_pr_call():
    assert LOCI_ARRAY_COLUMNS == ("array_id", "array_size", "array_members", "array_support",
                                  "array_support_reasons")
    pr = {"array_id": "P:PR:c1:1-9", "array_size": 2, "array_members": ["1-2:+", "5-9:-"],
          "array_support": SUPPORTED, "array_support_reasons": ["array_size>=2"]}
    assert loci_columns(pr) == {
        "array_id": "P:PR:c1:1-9", "array_size": 2, "array_members": "1-2:+|5-9:-",
        "array_support": SUPPORTED, "array_support_reasons": "array_size>=2"}
    assert set(loci_columns({}).values()) == {""}


# ---- the pipeline: calls do not change ---------------------------------

def _run(tmp_path):
    _order(tmp_path)
    seq = "C" * 1_000 + "A" * 900 + "C" * 2_000 + _orf() + "C" * 5_000
    genome = _genome(tmp_path, seq)

    def polish(**kw):
        if kw["gene_name"] == "pheromone_receptor":
            return _model("pheromone_receptor", "c1", 1_001, 1_900, family_key=KEY)
        return None

    return run_pipeline(
        genome_fasta=genome, proteome_fasta=None, taxid=None, db_root=tmp_path,
        reference_fasta=tmp_path / "reference.faa",
        search_fast_path=lambda *a, **k: [],
        search_localize=lambda *a, **k: [_receptor(1_001, 1_900)],
        polish_with_exonerate=polish, polish_with_miniprot=polish,
    )


def test_the_pipeline_flags_the_call_and_changes_nothing_else(tmp_path, monkeypatch):
    with_arrays = _run(tmp_path)
    monkeypatch.setattr(pipeline, "build_receptor_arrays", lambda families, hits: [])
    without = _run(tmp_path)
    assert len(with_arrays.results) == len(without.results) == 1
    [a], [b] = with_arrays.results, without.results
    # unchanged: tier, confidence, verification label, counts, every other field
    assert replace(a, receptor_array=None) == replace(b, receptor_array=None)
    assert (a.confidence, a.verification) == (b.confidence, b.verification)
    assert with_arrays.suppressed_loci == without.suppressed_loci
    # the flag: one receptor locus, one CAAX ORF -> an unsupported singleton
    assert a.receptor_array["array_size"] == 1
    assert a.receptor_array["array_support"] == UNSUPPORTED
    assert a.verification["status"] == "unverified"
    assert [x.n_calls for x in with_arrays.receptor_arrays] == [1]
    assert without.receptor_arrays == []
