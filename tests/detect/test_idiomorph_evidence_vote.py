"""The idiomorph is decided by the STRONGEST evidence, not by presence/absence.

Curator's rulings, 2026-09-21, from the Saccharomyces cassette work:

1. A gene can be marked `idiomorph_informative: false`. It still marks the
   cassette; it just stops voting on which allele the cassette carries --
   "just knowing where the cassettes are is sufficient".
2. Rank the idiomorphs by their best evidence and take the winner.
3. On an exact tie, REPORT BOTH rather than try to break it.

WHY, measured on S. cerevisiae S288C (tblastn against the three chromosome-III
cassettes, bitscore shown):

                    HML(alpha)   MAT(alpha)   HMR(a)
      MATA2             246          246        246
      MATALPHA2         422          422        374    <- wins everywhere
      MATALPHA1         365          365         82    (38 aa, 22% coverage)
      MATA1               -            -        236

`a2` and `alpha2` are homologous over the X region shared by all three
cassettes, so `MATA2` hits at 100% identity and 100% coverage everywhere and
`MATALPHA2`, being the longer protein, outscores it at EVERY cassette --
including HMR, which is actually `a`. Any rule that lets that pair vote calls
HMR alpha. Marking the pair non-informative leaves `MATA1` vs `MATALPHA1`,
which is right at all three.

Presence/absence could not work either: every cassette carries genes indicated
for both idiomorphs, so the old rule returned `undetermined` for all three.
"""
from MATPredict.detect.family_registry import Family, FamilyKey
from MATPredict.detect.idiomorph import assign_idiomorph, idiomorph_candidates
from MATPredict.detect.search import SearchHit

KEY = FamilyKey("Ascomycota", "MATsc")
FAM = Family(
    key=KEY, vocabulary_type="enum", idiomorph_values=["MATa", "MATalpha"],
    idiomorph_pattern=None,
    genes=[
        {"name": "MATALPHA1", "role": "core_MAT", "present_in_idiomorphs": ["MATalpha"]},
        {"name": "MATA1", "role": "core_MAT", "present_in_idiomorphs": ["MATa"]},
        {"name": "MATALPHA2", "role": "core_MAT", "present_in_idiomorphs": ["MATalpha"],
         "idiomorph_informative": False},
        {"name": "MATA2", "role": "core_MAT", "present_in_idiomorphs": ["MATa"],
         "idiomorph_informative": False},
    ],
    taxonomic_scope=[4893],
)


def _h(gene, bits, start=1000):
    return SearchHit(KEY, gene, "core_MAT", "c1", start, start + 300, "+",
                     100.0, "r", "tblastn_genome", bitscore=bits)


def test_an_alpha_cassette_is_called_alpha():
    hits = [_h("MATA2", 246), _h("MATALPHA2", 422), _h("MATALPHA1", 365)]
    names = [h.gene_name for h in hits]
    assert assign_idiomorph(FAM, names, hits) == "MATalpha"


def test_the_a_cassette_is_called_a_despite_alpha2_outscoring_everything():
    """The whole point. MATALPHA2 has the highest bitscore at HMR (374) but is
    non-informative, so the call comes from MATA1 (236) vs MATALPHA1 (82)."""
    hits = [_h("MATA2", 246), _h("MATALPHA2", 374),
            _h("MATA1", 236), _h("MATALPHA1", 82)]
    names = [h.gene_name for h in hits]
    assert assign_idiomorph(FAM, names, hits) == "MATa"


def test_presence_alone_would_have_said_undetermined():
    """Guard against regressing to the old rule: both idiomorphs are indicated
    by gene NAME at every cassette."""
    names = ["MATA2", "MATALPHA2", "MATA1", "MATALPHA1"]
    assert assign_idiomorph(FAM, names) == "undetermined"   # no hits -> old rule


def test_an_exact_tie_reports_both():
    hits = [_h("MATA1", 200), _h("MATALPHA1", 200)]
    names = [h.gene_name for h in hits]
    assert assign_idiomorph(FAM, names, hits) == "undetermined"
    assert idiomorph_candidates(FAM, names, hits) == [
        {"idiomorph": "MATa", "score": 200.0},
        {"idiomorph": "MATalpha", "score": 200.0},
    ]


def test_candidates_are_ranked_and_reported_even_when_there_is_a_winner():
    hits = [_h("MATA1", 236), _h("MATALPHA1", 82)]
    names = [h.gene_name for h in hits]
    assert idiomorph_candidates(FAM, names, hits) == [
        {"idiomorph": "MATa", "score": 236.0},
        {"idiomorph": "MATalpha", "score": 82.0},
    ]


def test_non_informative_genes_never_appear_as_candidates():
    hits = [_h("MATA2", 246), _h("MATALPHA2", 422)]
    names = [h.gene_name for h in hits]
    assert idiomorph_candidates(FAM, names, hits) == []
    assert assign_idiomorph(FAM, names, hits) == "undetermined"


def test_superseded_hits_do_not_vote():
    hits = [_h("MATA1", 500), _h("MATALPHA1", 82)]
    hits[0] = SearchHit(KEY, "MATA1", "core_MAT", "c1", 1000, 1300, "+", 100.0,
                        "r", "tblastn_genome", bitscore=500.0,
                        superseded_by="MATALPHA1")
    names = [h.gene_name for h in hits]
    assert assign_idiomorph(FAM, names, hits) == "MATalpha"


def test_identity_is_the_fallback_when_no_bitscore_is_recorded():
    """Polished models carry no BLAST statistics."""
    a = SearchHit(KEY, "MATA1", "core_MAT", "c1", 1, 300, "+", 90.0, "r", "exonerate_refine")
    b = SearchHit(KEY, "MATALPHA1", "core_MAT", "c1", 1, 300, "+", 40.0, "r", "exonerate_refine")
    assert assign_idiomorph(FAM, ["MATA1", "MATALPHA1"], [a, b]) == "MATa"


def test_a_pattern_family_is_still_undetermined():
    fam = Family(key=FamilyKey("B", "HD"), vocabulary_type="pattern",
                 idiomorph_values=None, idiomorph_pattern="^A[0-9]+$",
                 genes=[{"name": "HD1", "role": "core_MAT"}], taxonomic_scope=[1])
    assert assign_idiomorph(fam, ["HD1"], [_h("HD1", 100)]) == "undetermined"
