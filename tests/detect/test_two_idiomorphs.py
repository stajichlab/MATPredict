"""Genome-level two-idiomorph statement (curator's ruling 2026-09-29).

A genome with a determined Plus call and a determined Minus call in one
family is NOT asserted to be homothallic: the statement lists the arrangement,
the evidence, and every possible cause with the support each one has.
"""
from MATPredict.detect.two_idiomorphs import (
    CAUSES, two_idiomorph_statements, translate_exons,
)

FAM = "Mucoromycota:MAT"

# A 300-aa-ish protein-coding stretch: random-looking but fixed.
_PROT_A = ("MSTKRPLNAFMLFRSEMQKQHPELSNPEISKLLGERWRALSEEEKAPYYEEEARLKAEHAEKYPD"
           "YKYQPRRKTGSKRKKSSDSEQQHQQQPLLTQESSSSPQSSTPSPLQHNQYDSIASSSSSSAPTFQ")
_PROT_B = ("MSTKRPLNAFMLFRSEMQKQHPELSNPEISKLLGERWRALSEEEKAPYYEEEARLKAEHAEKYPD"
           "YKYQPRRKTGSKRKKSSDSEQQHQQQPLLTQESSSSPQSSTPSPLQHNQYDSIASSSSSSAPTFQ")
_CODON = {aa: c for c, aa in [
    ("GCT", "A"), ("CGT", "R"), ("AAT", "N"), ("GAT", "D"), ("TGT", "C"), ("CAA", "Q"),
    ("GAA", "E"), ("GGT", "G"), ("CAT", "H"), ("ATT", "I"), ("CTT", "L"), ("AAA", "K"),
    ("ATG", "M"), ("TTT", "F"), ("CCT", "P"), ("TCT", "S"), ("ACT", "T"), ("TGG", "W"),
    ("TAT", "Y"), ("GTT", "V")]}


def _dna(prot):
    return "".join(_CODON[a] for a in prot)


def _mutate(prot, every):
    """Change one residue in every `every` positions (A<->S), lowering identity."""
    out = list(prot)
    for i in range(0, len(out), every):
        out[i] = "S" if out[i] != "S" else "A"
    return "".join(out)


def _flank(name, contig, start, prot_len):
    return {"gene": name, "role": "flanking_conserved", "contig": contig,
            "start": start, "end": start + prot_len * 3 - 1, "strand": "+",
            "identity": 80.0, "status": "polished_agree",
            "exons": [(start, start + prot_len * 3 - 1)]}


def _call(idiomorph, contig, start, end, *, confidence="high", flanks=(),
          locus_class="mat_locus", classifier_input="model", margin=100.0):
    return {"family": FAM, "contig": contig, "start": start, "end": end,
            "idiomorph": idiomorph, "confidence": confidence,
            "locus_class": locus_class, "classifier_input": classifier_input,
            "margin": margin, "gene_evidence": list(flanks)}


def _one(calls, seqs=None):
    out = two_idiomorph_statements(calls, {FAM}, seqs or {}, genetic_code=1)
    assert len(out) == 1
    return out[0]


def test_no_statement_for_one_idiomorph_or_undetermined():
    calls = [_call("Plus", "c1", 1, 100), _call("undetermined", "c2", 1, 100)]
    assert two_idiomorph_statements(calls, {FAM}, {}, genetic_code=1) == []


def test_family_not_enabled_gives_no_statement():
    calls = [_call("Plus", "c1", 1, 100), _call("Minus", "c2", 1, 100)]
    assert two_idiomorph_statements(calls, {"Other:MAT"}, {}, genetic_code=1) == []


def test_arrangements():
    assert _one([_call("Plus", "c1", 1, 5000),
                 _call("Minus", "c2", 1, 5000)])["arrangement"] == "unlinked"
    assert _one([_call("Plus", "c1", 1, 5000),
                 _call("Minus", "c1", 150_000, 155_000)])["arrangement"] == "same_contig_distant"
    assert _one([_call("Plus", "c1", 1, 5000),
                 _call("Minus", "c1", 10_000, 15_000)])["arrangement"] == "same_locus"


def test_every_cause_is_listed_and_homothallism_is_never_asserted():
    st = _one([_call("Plus", "c1", 1, 5000), _call("Minus", "c2", 1, 5000)])
    assert [c["cause"] for c in st["possible_causes"]] == list(CAUSES)
    assert "homothallic" not in st  # never a verdict field
    assert st["supported_causes"] == []


def test_two_calls_of_one_idiomorph_support_duplication():
    st = _one([_call("Plus", "c1", 1, 5000), _call("Plus", "c3", 1, 5000),
               _call("Minus", "c2", 1, 5000)])
    assert "duplication" in st["supported_causes"]


def test_gc_difference_supports_mixed_culture():
    seqs = {"c1": "GC" * 5000, "c2": "AT" * 4000 + "GC" * 1000}
    st = _one([_call("Plus", "c1", 1, 5000), _call("Minus", "c2", 1, 5000)], seqs)
    assert "mixed_culture_or_heterokaryon" in st["supported_causes"]
    assert st["evidence"]["gc_difference_pct"] > 5


def test_near_identical_shared_flank_supports_fusion_not_homothallism():
    da = _dna(_PROT_A)
    seqs = {"c1": da + "A" * 6000, "c2": da + "A" * 6000}
    st = _one([_call("Plus", "c1", 1, 8000, flanks=[_flank("rnhA", "c1", 1, len(_PROT_A))]),
               _call("Minus", "c2", 1, 8000, flanks=[_flank("rnhA", "c2", 1, len(_PROT_A))])],
              seqs)
    pair = st["evidence"]["shared_flanks"][0]
    assert pair["gene"] == "rnhA" and pair["protein_identity"] >= 99.0
    assert "hybrid_or_fusion" in st["supported_causes"]
    assert "homothallism" not in st["supported_causes"]


def test_divergent_shared_flank_supports_homothallism():
    da, db = _dna(_PROT_A), _dna(_mutate(_PROT_B, 4))
    seqs = {"c1": da + "A" * 6000, "c2": db + "A" * 6000}
    st = _one([_call("Plus", "c1", 1, 8000, flanks=[_flank("rnhA", "c1", 1, len(_PROT_A))]),
               _call("Minus", "c2", 1, 8000, flanks=[_flank("rnhA", "c2", 1, len(_PROT_A))])],
              seqs)
    ident = st["evidence"]["shared_flanks"][0]["protein_identity"]
    assert 55.0 <= ident < 95.0
    assert "homothallism" in st["supported_causes"]


def test_same_locus_supports_homothallism():
    st = _one([_call("Plus", "c1", 1, 5000, locus_class="homothallic_candidate"),
               _call("Minus", "c1", 1, 5000, locus_class="homothallic_candidate")])
    assert st["arrangement"] == "same_locus"
    assert "homothallism" in st["supported_causes"]


def test_translate_exons_reverse_strand():
    dna = _dna("MKW")
    rc = dna[::-1].translate(str.maketrans("ACGT", "TGCA"))
    assert translate_exons(rc, [(1, 9)], "-", 1) == "MKW"


def test_pipeline_helper_only_runs_for_enabled_families_and_changes_no_call():
    from types import SimpleNamespace
    from MATPredict.detect.family_registry import FamilyKey
    from MATPredict.detect.pipeline import DetectionResult, _two_idiomorph_statements

    key = FamilyKey("Mucoromycota", "MAT")

    def res(idio, contig):
        return DetectionResult(family_key=key, contig=contig, start=1, end=5000,
                               confidence="high", idiomorph=idio, ambiguous_with=[],
                               genes_found=["sexP"], genes_missing=[], fragmented=False)

    results = [res("Plus", "c1"), res("Minus", "c2")]
    before = [(r.idiomorph, r.confidence) for r in results]
    on = [SimpleNamespace(key=key, two_idiomorphs_report=True)]
    off = [SimpleNamespace(key=key, two_idiomorphs_report=False)]
    assert _two_idiomorph_statements(results, off, None, 1) == []
    out = _two_idiomorph_statements(results, on, None, 1)
    assert len(out) == 1 and out[0]["arrangement"] == "unlinked"
    assert [(r.idiomorph, r.confidence) for r in results] == before
