"""Frameshift-aware translation of exonerate models (curator ruling 2026-10-04).

Measured case (notable finding 026): Mucor griseocyanus CBS 116.08 (MinION-only
assembly) has a 2-nt frameshift in an A5 homopolymer inside its sexM; the exon
CDS translated to a ~112-aa fragment and the call fell to undetermined. The
exonerate `Align` blocks give the frame-corrected CDS; block coordinates were
measured on exonerate 2.4.0: plus [t, t+len-1], minus [t-len, t-1].
"""
import types

from MATPredict.detect.pipeline import _translate_model
from MATPredict.detect.search import _exonerate_cds_blocks, _exonerate_frameshifts

def _write(tmp_path, seq):
    p = tmp_path / "g.fa"
    p.write_text(">c1\n" + seq + "\n")
    return p


def test_block_coordinates_follow_the_measured_convention():
    sim = ["c1", "exonerate", "similarity", "1", "1", ".", "+", ".",
           "alignment_id 1 ; Query q ; Align 101 1 30 ; Align 133 11 60"]
    assert _exonerate_cds_blocks(sim, "+", 0) == ((101, 130), (133, 192))
    sim[8] = "alignment_id 1 ; Query q ; Align 500 1 30 ; Align 467 11 60"
    assert _exonerate_cds_blocks(sim, "-", 1000) == ((1470, 1499), (1407, 1466))


def test_frameshift_count_is_summed_from_exon_attributes():
    exons = [["c", "e", "exon", "1", "9", ".", "+", ".", "insertions 0 ; frameshifts 2"],
             ["c", "e", "exon", "20", "29", ".", "+", ".", "insertions 0 ; deletions 0"]]
    assert _exonerate_frameshifts(exons) == 2


def test_a_frameshifted_model_is_translated_from_its_blocks(tmp_path):
    prot = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQ"
    # one fixed codon per amino acid
    codon = {a: c for c, a in zip(
        ["GCT","CGT","AAT","GAT","TGT","CAA","GAA","GGT","CAT","ATT","CTG","AAA","ATG","TTT","CCT","TCT","ACT","TGG","TAT","GTT"],
        "ARNDCQEGHILKMFPSTWYV")}
    cds = "".join(codon[a] for a in prot) + "TAA"
    k = 60  # insert 2 bases (a frameshift) after codon 20
    genome = "GG" * 10 + cds[:k] + "AA" + cds[k:] + "GG" * 10
    fa = _write(tmp_path, genome)
    start = 21
    end = start + len(cds) + 2 - 1
    model = types.SimpleNamespace(contig="c1", start=start, end=end, strand="+",
                                  exons=[types.SimpleNamespace(start=start, end=end)],
                                  frameshifts=2,
                                  cds_blocks=((start, start + k - 1), (start + k + 2, end)))
    assert _translate_model(fa, model, 1, {}) == prot
    model.frameshifts = 0
    assert len(_translate_model(fa, model, 1, {})) < len(prot)
