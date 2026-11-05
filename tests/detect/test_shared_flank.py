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


def test_the_role_follows_the_check_and_the_call_is_untouched():
    calls = [_call(SC, "MATalpha", "c1", 1000, 9000), _call(MTL, "A", "c2", 500, 2000)]
    for status, role in (("none", "cassette_undetermined"), ("shared_flank", "allelic_pair"), ("collinear", "duplicate_haplotype"), ("not_run", "cassette_undetermined")):
        out = fs.flag_cross_family_cassettes(calls, lambda a, b, s=status: {"status": s})
        assert out[0].cassette["role"] == role and out[0].cassette["shared_flank_check"]["status"] == status and out[1].cassette["role"] == role
        assert [r.idiomorph for r in out] == ["MATalpha", "A"]
    assert fs.flag_cross_family_cassettes(calls)[0].cassette["shared_flank_check"] == {"status": "not_run"}


@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not installed")
def test_blastn_finds_a_shared_flank_and_not_an_unrelated_contig(tmp_path):
    rnd = random.Random(7)
    dna = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
    flank_l, flank_r = dna(600), dna(600)
    alpha = dna(2000); a_idio = dna(1500)
    c20 = dna(5000) + flank_l + alpha + flank_r + dna(5000)
    mutate = lambda s: "".join(("A" if b != "A" else "C") if rnd.random() < 0.03 else b for b in s)      # about 3% divergence
    c274 = mutate(flank_l) + a_idio + mutate(flank_r)
    c99 = dna(4000)
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{n}\n{s}\n" for n, s in (("c20", c20), ("c274", c274), ("c99", c99))))
    a = SimpleNamespace(contig="c20", start=5601, end=8600)
    assert sf.shared_flank_check(fa, a, SimpleNamespace(contig="c274", start=1, end=len(c274)))["status"] == "shared_flank"
    assert sf.shared_flank_check(fa, a, SimpleNamespace(contig="c99", start=1, end=len(c99)))["status"] == "none"
    assert sf.shared_flank_check(tmp_path / "missing.fa", a, SimpleNamespace(contig="c99", start=1, end=10))["status"] == "not_run"
