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
