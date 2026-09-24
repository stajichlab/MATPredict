"""tblastn and diamond now retain the fields needed to judge an HSP.

Before this, a `SearchHit` carried `identity` and nothing else comparable:
no bitscore, no e-value, no alignment length, and `coverage=None` on every
tblastn hit. That is why a 31-codon, 40.6%-identity HSP was indistinguishable
from a real gene -- percent identity over an unrecorded length says almost
nothing.

`coverage` was NOT unavailable, it was never requested. tblastn 2.17.0 offers
`qcovhsp`, and because tblastn swaps the roles relative to the diamond fast
path -- here the CURATED PROTEIN is the query and the genome is the database
(`search.py`'s `_TBLASTN_OUTFMT` comment) -- `qcovhsp` is the same quantity
diamond's `scovhsp` gives: percent of the curated reference protein covered
by this HSP.

These fields are captured for calibration. Nothing ranks on them yet: the
existing `identity`-based `_rank` stays until a replacement is checked against
the 23/23 Mucoromycota ground truth.
"""
from MATPredict.detect.search import (
    _DIAMOND_OUTFMT,
    _TBLASTN_OUTFMT,
    SearchHit,
    search_localize,
)


def test_tblastn_asks_for_the_scoring_columns():
    cols = _TBLASTN_OUTFMT.split()
    for needed in ("qcovhsp", "evalue", "bitscore", "qlen", "length"):
        assert needed in cols, f"tblastn outfmt must request {needed}"


def test_diamond_asks_for_the_scoring_columns():
    for needed in ("evalue", "bitscore"):
        assert needed in _DIAMOND_OUTFMT, f"diamond outfmt must request {needed}"


def test_searchhit_carries_them_and_they_are_optional():
    """Optional because a polished model has no BLAST statistics of its own;
    `polish` constructs SearchHits directly."""
    hit = SearchHit(
        family_key=None, gene_name="g", role="core_MAT", contig="c1",
        start=1, end=9, strand="+", identity=50.0, reference_record_id="r",
        method="tblastn_genome",
    )
    assert hit.bitscore is None and hit.evalue is None
    assert hit.align_length_aa is None and hit.reference_length_aa is None


def _fake_tblastn(stdout):
    class R:
        returncode = 0
    def runner(cmd, **kw):
        r = R()
        r.stdout = stdout if cmd[0] == "tblastn" else ""
        r.stderr = ""
        return r
    return runner


def test_tblastn_hits_carry_coverage_from_qcovhsp(tmp_path):
    """The field that was previously hardcoded None on this path."""
    from MATPredict.detect.family_registry import Family, FamilyKey

    key = FamilyKey("P", "MAT")
    fam = Family(key=key, vocabulary_type="enum", idiomorph_values=["a"],
                 idiomorph_pattern=None,
                 genes=[{"name": "geneA", "role": "core_MAT"}],
                 taxonomic_scope=[1])
    # makeblastdb goes through the same runner and just needs to succeed.
    # qseqid sseqid pident length qlen qcovhsp evalue bitscore sstart send sframe
    row = "rec1|gene0|geneA\tcontig_1\t88.5\t120\t150\t80\t1e-40\t155.2\t500\t860\t1"
    hits = search_localize(
        tmp_path / "g.fa", [fam], tmp_path / "ref.faa", {"rec1": key},
        runner=_fake_tblastn(row + "\n"),
    )
    assert len(hits) == 1
    h = hits[0]
    assert h.coverage == 80.0, "coverage must come from qcovhsp, not be None"
    assert h.bitscore == 155.2
    assert h.evalue == 1e-40
    assert h.align_length_aa == 120
    assert h.reference_length_aa == 150
    assert h.identity == 88.5 and h.contig == "contig_1"
    assert (h.start, h.end, h.strand) == (500, 860, "+")
