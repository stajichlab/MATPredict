"""The flank-only helixer rule: flanking genes take helixer's ends; everything else keeps homology."""
from dataclasses import dataclass

from MATPredict.detect.helixer_rule import KEPT, REFINED, refine_gene, refine_genes
from MATPredict.detect.tracks import Track, TrackGene


@dataclass
class G:
    gene_name: str
    role: str
    contig: str
    start: int
    end: int
    strand: str


def hx(start, end, strand="+", contig="c1", complete=True, name="h1"):
    return TrackGene(name=name, contig=contig, start=start, end=end, strand=strand,
                     exons=((start, end),), source="helixer", complete_orf=complete)


def test_flank_gene_takes_helixer_ends_plus_strand():
    r = refine_gene(G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+"), Track("helixer", (hx(950, 2100),)))
    assert (r.action, r.reason, r.start, r.end) == (REFINED, "helixer_extension", 950, 2100)
    assert (r.original_start, r.original_end) == (1000, 2000)
    assert (r.delta_5p, r.delta_3p) == (50, 100)   # both ends grew
    assert r.helixer_gene == "h1"


def test_deltas_follow_the_gene_strand():
    r = refine_gene(G("APN2", "flanking_variable", "c1", 1000, 2000, "-"), Track("helixer", (hx(900, 2010, "-"),)))
    assert (r.delta_5p, r.delta_3p) == (10, 100)   # on - the 5' end is the right-hand end


def test_helixer_never_shortens_a_flank_model():
    r = refine_gene(G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+"), Track("helixer", (hx(1030, 1990),)))
    assert (r.action, r.reason, r.start, r.end) == (KEPT, "helixer_not_longer", 1000, 2000)


def test_each_end_moves_on_its_own():
    # helixer reaches 80 bp further at the 5' end but stops 40 bp short at the 3' end: only the 5' end moves
    r = refine_gene(G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+"), Track("helixer", (hx(920, 1960),)))
    assert (r.action, r.start, r.end, r.delta_5p, r.delta_3p) == (REFINED, 920, 2000, 80, 0)
    r = refine_gene(G("SLA2", "flanking_conserved", "c1", 1000, 2000, "-"), Track("helixer", (hx(920, 1960, "-"),)))
    assert (r.start, r.end, r.delta_5p, r.delta_3p) == (920, 2000, 0, 80)   # on - the low end is the 3' end


def test_core_mat_keeps_homology_even_when_helixer_disagrees():
    r = refine_gene(G("sexM", "core_MAT", "c1", 1000, 2000, "+"), Track("helixer", (hx(500, 3000),)))
    assert (r.action, r.reason, r.start, r.end) == (KEPT, "homology_stratum", 1000, 2000)


def test_missing_helixer_model_keeps_homology():
    r = refine_gene(G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+"), Track("helixer", ()))
    assert (r.action, r.reason) == (KEPT, "no_helixer_model")


def test_other_strand_or_contig_or_small_overlap_is_not_a_match():
    g = G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+")
    for model in (hx(1000, 2000, "-"), hx(1000, 2000, contig="c2"), hx(1900, 5000)):
        assert refine_gene(g, Track("helixer", (model,))).reason == "no_helixer_model"


def test_incomplete_orf_blocks_but_unknown_does_not():
    g = G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+")
    assert refine_gene(g, Track("helixer", (hx(950, 2100, complete=False),))).reason == "helixer_orf_incomplete"
    assert refine_gene(g, Track("helixer", (hx(950, 2100, complete=None),))).action == REFINED


def test_identical_ends_are_recorded_as_agreeing():
    r = refine_gene(G("SLA2", "flanking_conserved", "c1", 1000, 2000, "+"), Track("helixer", (hx(1000, 2000),)))
    assert (r.action, r.reason) == (KEPT, "ends_agree")


def test_largest_overlap_wins_and_one_row_per_gene():
    genes = [G("a", "flanking_conserved", "c1", 1000, 2000, "+"), G("b", "core_MAT", "c1", 5000, 5600, "+")]
    track = Track("helixer", (hx(1500, 2400, name="small"), hx(990, 2020, name="big")))
    rows = refine_genes(genes, track)
    assert [r.gene for r in rows] == ["a", "b"] and rows[0].helixer_gene == "big"
