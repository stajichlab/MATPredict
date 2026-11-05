"""`shared_flank_check`: cross-family cassette pairs that share flanking sequence are one heterozygous locus, not two cassettes (#189)."""
from __future__ import annotations

import random
import shutil
from dataclasses import replace
from types import SimpleNamespace

import pytest

from MATPredict.detect import family_status as fs
from MATPredict.detect import shared_flank as sf
from MATPredict.detect.family_registry import FamilyKey
from MATPredict.detect.pipeline import DetectionResult

SC, MTL = FamilyKey("Ascomycota", "MATsc"), FamilyKey("Ascomycota", "MTL")


def _call(key, idiomorph, contig, start, end):
    return DetectionResult(family_key=key, contig=contig, start=start, end=end, confidence="high", idiomorph=idiomorph, ambiguous_with=[],
                           genes_found=[], genes_missing=[], fragmented=False, gene_evidence=[])


def test_classify_separates_shared_flank_collinear_and_unrelated():
    assert sf.classify([], 20000, 3000)["status"] == "none"
    assert sf.classify([(97.2, 533, 100, 632, 1, 533), (97.3, 147, 900, 1046, 3361, 3507)], 40000, 3500)["status"] == "shared_flank"
    assert sf.classify([(60.0, 300, 1, 300, 1, 300)], 20000, 3000)["status"] == "none"          # too divergent
    assert sf.classify([(98.0, 80, 1, 80, 1, 80)], 20000, 3000)["status"] == "none"             # too short
    assert sf.classify([(95.0, 12000, 1, 12000, 1, 12000)], 40000, 20000)["status"] == "collinear"


def test_both_sides_needs_shared_sequence_left_and_right_of_the_idiomorph_in_both_loci():
    left, right = (97.0, 401, 100, 500, 1, 400), (97.0, 147, 900, 1046, 3361, 3507)
    core_a, core_b = (5601, 7600), (601, 2100)       # contig coordinates; windows start at contig 5001-... and 1
    # window A starts 5000 before core A's first residue: hits at window 100-500 and 900-1046 are contig 5100-5500 and 5900-6046: both LEFT of core A
    assert not sf.both_sides([left, right], 5000, 0, core_a, core_b)
    # a hit right of core A as well
    right_a = (97.0, 147, 3100, 3246, 3361, 3507)                                     # contig 8100-8246 (right of 7600); B side 3361-3507 is right of 2100
    assert sf.both_sides([left, right_a], 5000, 0, core_a, core_b)
    assert not sf.both_sides([left], 5000, 0, core_a, core_b)                      # one side only (the Williopsis-like single hit)
    assert not sf.both_sides([(97.0, 80, 100, 180, 1, 80), right_a], 5000, 0, core_a, core_b)   # a hit under 100 bp does not count


def test_the_role_stays_cassette_undetermined_unless_both_sides_are_shared_and_it_is_not_a_repeat():
    calls = [_call(SC, "MATalpha", "c1", 1000, 9000), _call(MTL, "A", "c2", 500, 2000)]
    cases = ((dict(status="none"), "cassette_undetermined"),
             (dict(status="shared_flank", both_sides=False, repeat=None), "cassette_undetermined"),     # one short hit
             (dict(status="collinear", both_sides=False, repeat=None), "cassette_undetermined"),        # measured only
             (dict(status="shared_flank", both_sides=True, repeat=True), "cassette_undetermined"),      # a repeat
             (dict(status="shared_flank", both_sides=True, repeat=None), "cassette_undetermined"),      # repeat check not possible
             (dict(status="shared_flank", both_sides=True, repeat=False), "allelic_pair"),
             (dict(status="not_run"), "cassette_undetermined"))
    for check, role in cases:
        out = fs.flag_cross_family_cassettes(calls, lambda a, b, c=check: dict(c))
        assert out[0].cassette["role"] == role and out[1].cassette["role"] == role
        assert out[0].cassette["shared_flank_check"]["status"] == check["status"] and out[0].cassette["shared_flank_check"]["partner_family"] == "Ascomycota:MTL"
        assert [r.idiomorph for r in out] == ["MATalpha", "A"]
    assert fs.flag_cross_family_cassettes(calls)[0].cassette["shared_flank_check"] == {"status": "not_run"}


def test_every_partner_is_checked_and_the_allelic_one_decides():
    third = FamilyKey("Ascomycota", "MATog")
    calls = [_call(SC, "MATalpha", "c1", 1000, 9000), _call(MTL, "A", "c2", 500, 2000), _call(third, "MATa", "c3", 10, 900)]
    def check(a, b):
        return dict(status="shared_flank", both_sides=True, repeat=False) if b.family_key == third else dict(status="none")
    out = fs.flag_cross_family_cassettes(calls, check)
    assert out[0].cassette["role"] == "allelic_pair" and out[0].cassette["shared_flank_check"]["partner_family"] == "Ascomycota:MATog"


@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not installed")
def test_blastn_reports_both_sides_for_an_allele_pair_and_a_repeat_for_a_multicopy_flank(tmp_path):
    rnd = random.Random(7)
    dna = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    flank_l, flank_r, alpha, a_idio = dna(600), dna(600), dna(2000), dna(1500)
    c20 = dna(5000) + flank_l + alpha + flank_r + dna(5000)                                   # alpha core 5601-7600
    mutate = lambda s: "".join(("A" if b != "A" else "C") if rnd.random() < 0.03 else b for b in s)      # about 3% divergence
    c274 = mutate(flank_l) + a_idio + mutate(flank_r)                                          # a core 601-2100
    c99 = dna(4000)
    a = SimpleNamespace(contig="c20", start=5601, end=7600)
    b = SimpleNamespace(contig="c274", start=601, end=2100)
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{n}\n{x}\n" for n, x in (("c20", c20), ("c274", c274), ("c99", c99))))
    got = sf.shared_flank_check(fa, a, b)
    assert got["status"] == "shared_flank" and got["both_sides"] is True and got["repeat"] is False
    fb = tmp_path / "g3.fa"                                                                     # the left flank also sits elsewhere: a repeat
    fb.write_text("".join(f">{n}\n{x}\n" for n, x in (("c20", c20), ("c274", c274), ("c99", c99 + flank_l + dna(500) + flank_l + dna(300)))))
    assert sf.shared_flank_check(fb, a, b)["repeat"] is True

    assert sf.shared_flank_check(fa, a, SimpleNamespace(contig="c274", start=2701, end=3000))["both_sides"] is False   # a core with nothing shared to its right
    assert sf.shared_flank_check(fa, a, SimpleNamespace(contig="c99", start=1, end=4000))["status"] == "none"
    assert sf.shared_flank_check(tmp_path / "missing.fa", a, b)["status"] == "not_run"
