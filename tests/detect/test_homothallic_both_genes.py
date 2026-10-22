"""`--homothallic-policy both_genes`: a complete, profile-confirmed alpha-box gene and a complete, profile-confirmed
HMG-box gene at one locus is a `homothallic_candidate`; fragments, HMG paralogs and off-locus genes are not; the
policy is off by default. Pure rule tests plus a synthetic-genome test on both strands with the real profiles."""
from argparse import Namespace

import pytest

from MATPredict.detect import homothallic as hm
from MATPredict.detect.homothallic import CoreModel, evaluate


def _m(gene, start, end, alpha=0.0, hmg=0.0, aa=300, cds=300, contig="c1", lost_to=None):
    return CoreModel(gene, contig, start, end, aa, cds, {"alpha_box": alpha, "HMG_box": hmg}, lost_to)


ALPHA = dict(alpha=180.0, hmg=3.0)
HMG = dict(alpha=4.0, hmg=170.0)


def test_policy_defaults_to_off():
    from MATPredict.detect.cli import _homothallic_policy_from_args

    assert hm.DEFAULT_HOMOTHALLIC_POLICY == hm.POLICY_OFF == "off"
    assert _homothallic_policy_from_args(Namespace()) == "off"
    assert _homothallic_policy_from_args(Namespace(homothallic_policy="both_genes")) == "both_genes"


def test_cli_flag_parses_and_defaults_off():
    import argparse

    from MATPredict.detect.cli import register_subcommands

    parser = argparse.ArgumentParser()
    register_subcommands(parser.add_subparsers(dest="cmd"))
    base = ["detect", "--genome", "g.fa", "--taxid", "1"]
    assert parser.parse_args(base).homothallic_policy == "off"
    assert parser.parse_args(base + ["--homothallic-policy", "both_genes"]).homothallic_policy == "both_genes"


def test_two_complete_profile_confirmed_genes_make_a_candidate():
    note = evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), _m("MAT1-2-1", 3000, 4100, **HMG)], 20_000)
    assert note["decision"] == "homothallic_candidate" and note["separation_bp"] == 800
    assert note["alpha_box_gene"]["gene"] == "MAT1-1-1" and note["hmg_box_gene"]["gene"] == "MAT1-2-1"


def test_an_hmg_gene_that_lost_a_resolution_still_counts_when_the_profile_confirms_it():
    # the Pyrenophora tritici-repentis pattern: MAT1-2-1 set aside for MAT1-1-3, the protein is the MAT1-2-1 HMG box
    note = evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), _m("MAT1-2-1", 3000, 4100, lost_to="MAT1-1-3", **HMG)], 20_000)
    assert note is not None and note["hmg_box_gene"]["lost_to"] == "MAT1-1-3"


def test_an_hmg_paralog_below_the_mat1_2_1_floor_does_not_trigger():
    # curated MAT1-1-3 scores 98-101 bits on the HMG profile; the floor is 110
    assert evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), _m("MAT1-2-1", 3000, 4100, alpha=2.0, hmg=99.0)], 20_000) is None


def test_an_hmg_gene_off_the_locus_does_not_trigger():
    assert evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), _m("MAT1-2-1", 400_000, 401_100, **HMG)], 20_000) is None
    assert evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), _m("MAT1-2-1", 3000, 4100, contig="c2", **HMG)], 20_000) is None


def test_overlapping_models_are_one_region_not_two_genes():
    assert evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), _m("MAT1-2-1", 2100, 3100, **HMG)], 20_000) is None


def test_a_fragment_does_not_trigger():
    short = _m("MAT1-2-1", 3000, 3300, aa=90, cds=100, **HMG)          # 90 aa ORF
    truncated = _m("MAT1-2-1", 3000, 4100, aa=150, cds=366, **HMG)    # stop inside a 366-codon model
    for frag in (short, truncated):
        assert evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA), frag], 20_000) is None


def test_one_idiomorph_alone_and_unconfirmed_proteins_do_not_trigger():
    assert evaluate([_m("MAT1-1-1", 1000, 2200, **ALPHA)], 20_000) is None
    # an alpha-box protein that the HMG profile also scores closely is not confirmed as either
    assert evaluate([_m("MAT1-1-1", 1000, 2200, alpha=150.0, hmg=140.0), _m("MAT1-2-1", 3000, 4100, **HMG)], 20_000) is None


# ---- synthetic genome, both strands, real shipped profiles -----------------------------------------------------
_CODON = {"A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "I": "ATT", "K": "AAA",
          "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA", "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT",
          "W": "TGG", "Y": "TAT"}
_COMP = str.maketrans("ACGT", "TGCA")


def _protein(db_root, record, gene):
    faa = db_root / "Ascomycota" / "Mycosphaerellales" / record / "proteins.faa"
    seq, keep = [], False
    for line in faa.read_text().splitlines():
        if line.startswith(">"):
            keep = f"name={gene}" in line
        elif keep:
            seq.append(line.strip())
    return "".join(seq)


def _dna(protein):
    return "".join(_CODON[a] for a in protein) + "TAA"


def _family_and_profiles(real_db_root):
    pytest.importorskip("pyhmmer")
    from MATPredict.detect import partial_core
    from MATPredict.detect.family_registry import FamilyKey, load_all_families

    family = next(f for f in load_all_families(real_db_root) if f.key == FamilyKey("Ascomycota", "MAT"))
    clf = partial_core.load_core_profiles(real_db_root, family)
    if clf is None:
        pytest.skip("alpha/HMG profiles not available")
    return family, clf


def _evidence(gene, contig, start, end, strand, lost_to=None):
    from MATPredict.detect.pipeline import GeneEvidence

    return GeneEvidence(gene_name=gene, role="core_MAT", contig=contig, start=start, end=end, strand=strand,
                        identity=50.0, coverage=100.0, reference_record_id="rec", method="exonerate_refine",
                        status="polished_agree", exons=((start, end),), lost_to=lost_to)


@pytest.mark.parametrize("alpha_strand,hmg_strand", [("+", "-"), ("-", "+"), ("+", "+"), ("-", "-")])
def test_synthetic_locus_on_both_strands_is_a_candidate(real_db_root, tmp_path, alpha_strand, hmg_strand):
    from MATPredict.detect.pipeline import _translate_model

    family, clf = _family_and_profiles(real_db_root)
    alpha_dna = _dna(_protein(real_db_root, "1047171_ipo323_MAT_MAT1-1", "MAT1-1-1"))
    hmg_dna = _dna(_protein(real_db_root, "1047171_ipo94269_MAT_MAT1-2", "MAT1-2-1"))

    def placed(dna, strand):
        return dna if strand == "+" else dna.translate(_COMP)[::-1]

    spacer = "ACGTTGCA" * 100
    a_seq, h_seq = placed(alpha_dna, alpha_strand), placed(hmg_dna, hmg_strand)
    genome = "GATTACA" * 50 + a_seq + spacer + h_seq + "GATTACA" * 50
    a_start = 7 * 50 + 1
    a_end = a_start + len(a_seq) - 1
    h_start = a_end + len(spacer) + 1
    h_end = h_start + len(h_seq) - 1
    fasta = tmp_path / "g.fa"
    fasta.write_text(">c1\n" + genome + "\n")
    evidence = [_evidence("MAT1-1-1", "c1", a_start, a_end, alpha_strand),
                _evidence("MAT1-2-1", "c1", h_start, h_end, hmg_strand, lost_to="MAT1-1-3")]
    cache: dict = {}
    note = hm.note_for_evidence(evidence, family, clf, lambda span: _translate_model(fasta, span, 1, cache),
                                family.max_homothallic_separation_bp)
    assert note is not None, "both complete genes, on either strand, must make a candidate"
    assert note["alpha_box_gene"]["alpha_bits"] >= hm.MIN_PROFILE_BITS["alpha_box"]
    assert note["hmg_box_gene"]["hmg_bits"] >= hm.MIN_PROFILE_BITS["HMG_box"]
    assert note["separation_bp"] == len(spacer) + 1  # start of the later gene minus end of the earlier, as classify_locus


def test_synthetic_hmg_fragment_and_off_locus_hmg_do_not_trigger(real_db_root, tmp_path):
    from MATPredict.detect.pipeline import _translate_model

    family, clf = _family_and_profiles(real_db_root)
    alpha_dna = _dna(_protein(real_db_root, "1047171_ipo323_MAT_MAT1-1", "MAT1-1-1"))
    hmg_protein = _protein(real_db_root, "1047171_ipo94269_MAT_MAT1-2", "MAT1-2-1")
    far = "ACGT" * 50_000                      # 200 kb between the genes
    genome = alpha_dna + "ACGTTGCA" * 100 + _dna(hmg_protein[:90]) + far + _dna(hmg_protein)
    fasta = tmp_path / "g.fa"
    fasta.write_text(">c1\n" + genome + "\n")
    a_end = len(alpha_dna)
    frag_start = a_end + 800 + 1
    frag_end = frag_start + 90 * 3 + 2
    far_start = frag_end + len(far) + 1
    far_end = far_start + len(hmg_protein) * 3 + 2
    cache: dict = {}
    translate = lambda span: _translate_model(fasta, span, 1, cache)  # noqa: E731
    fragment = [_evidence("MAT1-1-1", "c1", 1, a_end, "+"), _evidence("MAT1-2-1", "c1", frag_start, frag_end, "+")]
    assert hm.note_for_evidence(fragment, family, clf, translate, family.max_homothallic_separation_bp) is None
    off_locus = [_evidence("MAT1-1-1", "c1", 1, a_end, "+"), _evidence("MAT1-2-1", "c1", far_start, far_end, "+")]
    assert hm.note_for_evidence(off_locus, family, clf, translate, family.max_homothallic_separation_bp) is None


def test_length_floor_is_190_aa():
    assert hm.MIN_PROTEIN_AA == 190
    ok = _m("MAT1-1-1", 1000, 2200, aa=196, cds=196, **ALPHA)      # Pyrenophora-sized MAT1-1-1
    short = _m("MAT1-1-1", 1000, 2200, aa=189, cds=189, **ALPHA)
    assert evaluate([ok, _m("MAT1-2-1", 3000, 4100, **HMG)], 20_000) is not None
    assert evaluate([short, _m("MAT1-2-1", 3000, 4100, **HMG)], 20_000) is None


def _alpha(identity, coverage):
    import dataclasses

    return dataclasses.replace(_m("MAT1-1-1", 1000, 2200, **ALPHA), identity=identity, coverage=coverage)


def test_truncated_copy_of_a_known_gene_is_a_remnant_not_a_candidate():
    hmg = _m("MAT1-2-1", 3000, 4100, **HMG)
    assert evaluate([_alpha(92.2, 62.0), hmg], 20_000) is None      # Botrytis: close to the reference, 62% of it


def test_divergent_gene_with_low_coverage_is_kept():
    hmg = _m("MAT1-2-1", 3000, 4100, **HMG)
    assert evaluate([_alpha(62.4, 58.0), hmg], 20_000) is not None  # Didymella: no close reference, coverage says nothing


def test_full_length_close_gene_and_unknown_evidence_are_kept():
    hmg = _m("MAT1-2-1", 3000, 4100, **HMG)
    assert evaluate([_alpha(96.4, 98.5), hmg], 20_000) is not None  # S. minor
    assert evaluate([_alpha(None, None), hmg], 20_000) is not None


def test_remnant_thresholds_are_85_identity_and_80_coverage():
    from MATPredict.detect.homothallic import is_remnant

    assert is_remnant(_alpha(85.0, 79.9)) and not is_remnant(_alpha(84.9, 50.0)) and not is_remnant(_alpha(95.0, 80.0))


def test_note_for_evidence_passes_identity_and_coverage():
    import types

    from MATPredict.detect.homothallic import note_for_evidence

    class FakeFamily:
        genes = [{"name": "MAT1-1-1", "gene_class": "alpha_box"}, {"name": "MAT1-2-1", "gene_class": "HMG_box"}]

    clf = types.SimpleNamespace(models={0: ("MAT1-1-1", None), 1: ("MAT1-2-1", None)})
    import MATPredict.detect.classifier as c
    p_a, p_h = "A" * 300, "H" * 300
    scores = {p_a: {0: 180.0, 1: 3.0}, p_h: {0: 4.0, 1: 170.0}}
    orig = c.score_proteins
    c.score_proteins = lambda _clf, prots: scores[prots[0]]
    try:
        def ev(gene, start, end, ident, cov):
            return types.SimpleNamespace(role="core_MAT", status="polished_agree", gene_name=gene, contig="c1", start=start,
                                         end=end, strand="+", exons=[(start, end)], lost_to=None, identity=ident, coverage=cov)
        a, h = ev("MAT1-1-1", 1000, 1000 + 3 * 300 - 1, 92.0, 62.0), ev("MAT1-2-1", 3000, 3000 + 3 * 300 - 1, 99.0, 100.0)
        translate = lambda span: p_a if span.start == 1000 else p_h
        assert note_for_evidence([a, h], FakeFamily, clf, translate, 20_000) is None
        a.identity, a.coverage = 60.0, 58.0
        assert note_for_evidence([a, h], FakeFamily, clf, translate, 20_000) is not None
    finally:
        c.score_proteins = orig


def _candidate_result(idiomorph="MAT1-1"):
    from MATPredict.detect.family_registry import FamilyKey
    from MATPredict.detect.pipeline import DetectionResult, GeneEvidence, LocusSegment

    key = FamilyKey("Ascomycota", "MAT")
    genes = [GeneEvidence("MAT1-1-1", "core_MAT", "c1", 1000, 2200, "+", 62.0, 50.0, "rec_MAT_1", "exonerate_refine"),
             GeneEvidence("MAT1-2-1", "core_MAT", "c1", 3000, 4100, "+", 47.0, 68.0, "rec_MAT_2", "exonerate_refine")]
    return DetectionResult(
        family_key=key, contig="c1", start=900, end=4200, confidence="medium", idiomorph=idiomorph, ambiguous_with=[],
        genes_found=[g.gene_name for g in genes], genes_missing=[], fragmented=False,
        segments=[LocusSegment("c1", 900, 4200, contig_edge_distance=99)], gene_evidence=genes,
        reference_records=["rec_MAT_1"], locus_class="mat_locus",
    )


def test_a_homothallic_candidate_reports_undetermined_and_keeps_the_vote_in_the_note():
    from MATPredict.detect.pipeline import _as_homothallic_candidate

    note = {"decision": "homothallic_candidate", "policy": "both_genes"}
    r = _as_homothallic_candidate(_candidate_result("MAT1-1"), note)
    assert r.locus_class == "homothallic_candidate"
    assert r.idiomorph == "undetermined"
    assert r.homothallic_note["vote_idiomorph"] == "MAT1-1"
    assert r.homothallic_note["previous_locus_class"] == "mat_locus"
    # nothing else about the call changes
    base = _candidate_result("MAT1-1")
    assert (r.confidence, r.genes_found, r.start, r.end) == (base.confidence, base.genes_found, base.start, base.end)


def test_the_vote_winner_is_kept_whichever_idiomorph_it_was():
    from MATPredict.detect.pipeline import _as_homothallic_candidate

    for vote in ("MAT1-2", "undetermined"):
        r = _as_homothallic_candidate(_candidate_result(vote), {"decision": "homothallic_candidate"})
        assert r.idiomorph == "undetermined" and r.homothallic_note["vote_idiomorph"] == vote


def test_idiomorph_class_is_none_for_an_undetermined_candidate():
    from MATPredict.detect.pipeline import idiomorph_class_for

    class Fam:
        idiomorph_classes = {"MAT1-1": "alpha", "MAT1-2": "HMG"}

    assert idiomorph_class_for(Fam(), "undetermined") is None
