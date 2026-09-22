"""`distinct_intervals` -- how many separate places on the genome a call's
evidence actually occupies.

Reported as a FIELD, deliberately not used as a gate. Measured on 334
Tremellales genomes after the pheromone dedup:

    high    n=  134   single-interval=    0   median intervals 4
    medium  n=  435   single-interval=    0   median intervals 2
    low     n=2,157   single-interval=2,157   median intervals 1

Every low call is one ~150 bp HSP; every high and medium call occupies two or
more places. Gene COUNT does not separate these -- one ORF used to carry three
gene names -- but interval count does.

It is a field first so the curator has calibration data before any threshold
is set, and because an interval gate on REPORTING would lose the case the
`EvidenceFloor` docstring says must not be lost: a genuine lone gene on a tiny
contig in a fragmented assembly.
"""
from MATPredict.detect.clustering import GeneCluster
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.search import SearchHit
from MATPredict.detect.scoring import count_distinct_intervals

KEY = FamilyKey("P", "MAT")


def _hit(gene, start, end, contig="c1", superseded=None):
    return SearchHit(
        family_key=KEY, gene_name=gene, role="core_MAT", contig=contig,
        start=start, end=end, strand="+", identity=50.0,
        reference_record_id="r1", method="tblastn_genome",
        superseded_by=superseded,
    )


def test_three_gene_names_on_one_interval_count_as_one():
    """The measured failure: MFa1/MFa2/MFa3 on one 95 bp ORF."""
    hits = [_hit("MFa1", 100, 195), _hit("MFa2", 100, 195), _hit("MFa3", 100, 195)]
    assert count_distinct_intervals(hits) == 1


def test_genes_far_apart_count_separately():
    hits = [_hit("SXI1", 100, 1400), _hit("STE3", 70_000, 71_300)]
    assert count_distinct_intervals(hits) == 2


def test_tandem_duplicates_are_two_intervals_not_one():
    """Tandem pheromone copies are real and ADJACENT, not overlapping, so an
    interval rule must count them separately. This is why interval-counting is
    safe for duplication where name-collapsing would not be."""
    hits = [_hit("MFa", 100, 195), _hit("MFa", 400, 495)]
    assert count_distinct_intervals(hits) == 2


def test_partly_overlapping_hits_of_one_gene_are_one_interval():
    """Transitive grouping at >=50% of the shorter hit, the same test
    `idiomorph._overlap_groups` already uses."""
    hits = [_hit("SXI1", 100, 200), _hit("SXI1", 150, 260)]
    assert count_distinct_intervals(hits) == 1


def test_different_contigs_never_share_an_interval():
    hits = [_hit("SXI1", 100, 200), _hit("SXI2", 100, 200, contig="c2")]
    assert count_distinct_intervals(hits) == 2


def test_superseded_hits_do_not_count():
    """A resolved idiomorph cross-match is one gene seen twice, and the loser
    is kept in the report as evidence. It must not inflate the interval count."""
    hits = [_hit("sexP", 100, 300), _hit("sexM", 100, 300, superseded="sexP")]
    assert count_distinct_intervals(hits) == 1


def test_no_hits_is_zero():
    assert count_distinct_intervals([]) == 0


def test_the_field_reaches_the_report():
    from MATPredict.detect.report import _result_doc  # noqa: F401
