"""reads-type: idiomorph typing from raw reads by unique k-mers."""
from __future__ import annotations

import gzip
import random

import pytest

from MATPredict.reads.kmers import canonical_kmers, iter_fastq_reads, revcomp
from MATPredict.reads.panel import Panel
from MATPredict.reads.typing import call_idiomorphs, type_reads

K = 21
ABSENT = 0.2


def _rand(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


FLANK_L = _rand(300, 1)
FLANK_R = _rand(300, 2)
MAT1 = FLANK_L + _rand(2000, 3) + FLANK_R
MAT2 = FLANK_L + _rand(1800, 4) + FLANK_R


def _reads(seq: str, depth: float, length: int = 100, seed: int = 0, rc_half: bool = True):
    rng = random.Random(seed)
    n = int(depth * len(seq) / length)
    out = []
    for i in range(n):
        s = rng.randrange(0, len(seq) - length + 1)
        r = seq[s:s + length]
        if rc_half and i % 2:
            r = revcomp(r)
        out.append(r)
    return out


def _write_fastq(path, reads):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for i, r in enumerate(reads):
            fh.write(f"@r{i}\n{r}\n+\n{'I' * len(r)}\n")


def _panel() -> Panel:
    return Panel.from_sequences({"MAT1-1": MAT1, "MAT1-2": MAT2}, k=K)


def test_revcomp_and_canonical_are_strand_independent():
    assert revcomp("AACG") == "CGTT"
    seq = _rand(80, 9)
    assert canonical_kmers(seq, K) == canonical_kmers(revcomp(seq), K)


def test_panel_unique_kmers_exclude_shared_flanks():
    p = _panel()
    flank_kmers = canonical_kmers(FLANK_L, K) | canonical_kmers(FLANK_R, K)
    assert p.shared and p.shared >= flank_kmers
    for name in ("MAT1-1", "MAT1-2"):
        assert not (p.unique[name] & p.shared)
    assert not (p.unique["MAT1-1"] & p.unique["MAT1-2"])


def test_panel_needs_two_idiomorphs():
    with pytest.raises(ValueError):
        Panel.from_sequences({"MAT1-1": MAT1}, k=K)


def test_iter_fastq_reads_plain_and_gz(tmp_path):
    reads = _reads(MAT1, 1, seed=1)[:5]
    for name in ("a.fq", "a.fq.gz"):
        _write_fastq(tmp_path / name, reads)
        assert list(iter_fastq_reads(tmp_path / name)) == reads


@pytest.mark.parametrize("truth,expected", [("MAT1-1", "MAT1-1"), ("MAT1-2", "MAT1-2")])
def test_single_idiomorph_strain(tmp_path, truth, expected):
    seq = MAT1 if truth == "MAT1-1" else MAT2
    fq = tmp_path / "r.fq.gz"
    _write_fastq(fq, _reads(seq, 20, seed=5))
    res = type_reads(_panel(), [fq])
    assert res.call == expected
    assert res.breadth[expected] > 0.95
    other = "MAT1-2" if expected == "MAT1-1" else "MAT1-1"
    assert res.breadth[other] < 0.02


def test_both_idiomorphs_reports_depth_ratio(tmp_path):
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(MAT1, 10, seed=1) + _reads(MAT2, 30, seed=2))
    res = type_reads(_panel(), [fq])
    assert res.call == "both"
    ratio = res.depth["MAT1-2"] / res.depth["MAT1-1"]
    assert 2.0 < ratio < 4.5


def test_flanks_only_gives_none(tmp_path):
    """Reference gap: the shared flank control is covered, neither idiomorph is."""
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(FLANK_L, 40, seed=3) + _reads(FLANK_R, 40, seed=4))
    assert type_reads(_panel(), [fq]).call == "none"


def test_unrelated_reads_have_no_depth_control(tmp_path):
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(_rand(5000, 77), 20, seed=3))
    assert type_reads(_panel(), [fq]).call == "low_depth"


def test_diverged_reads_do_not_make_a_second_call(tmp_path):
    """A MAT1-2 strain plus a few percent-diverged MAT1-1-like reads (the 50a pattern)."""
    rng = random.Random(8)
    noisy = []
    for r in _reads(MAT1, 0.3, seed=4):
        r = list(r)
        for pos in rng.sample(range(len(r)), 8):
            r[pos] = rng.choice([b for b in "ACGT" if b != r[pos]])
        noisy.append("".join(r))
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(MAT2, 30, seed=2) + noisy)
    res = type_reads(_panel(), [fq])
    assert res.call == "MAT1-2"
    assert res.breadth["MAT1-1"] < ABSENT


def test_diverged_reference_allele_still_calls(tmp_path):
    """Strain allele differs from the panel reference by a SNP every 60 bp (about a third of the k-mers lost)."""
    allele = list(MAT2)
    for i in range(20, len(allele), 60):
        allele[i] = "A" if allele[i] != "A" else "C"
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads("".join(allele), 30, seed=2))
    res = type_reads(_panel(), [fq])
    assert res.call == "MAT1-2"
    assert res.breadth["MAT1-2"] < 0.9


def test_low_depth_still_calls(tmp_path):
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(MAT2, 4, seed=6))
    res = type_reads(_panel(), [fq])
    assert res.call == "MAT1-2"


def test_max_reads_caps_input(tmp_path):
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(MAT2, 20, seed=2))
    res = type_reads(_panel(), [fq], max_reads=10)
    assert res.reads_used == 10


def test_call_rule_is_explicit():
    assert call_idiomorphs({"A": 0.99, "B": 0.0}, {"A": 5.0, "B": 0.0}, 5.0) == ("A", [])
    assert call_idiomorphs({"A": 0.99, "B": 0.98}, {"A": 5.0, "B": 4.0}, 5.0)[0] == "both"
    assert call_idiomorphs({"A": 0.0, "B": 0.0}, {"A": 0.0, "B": 0.0}, 5.0)[0] == "none"
    assert call_idiomorphs({"A": 0.7, "B": 0.0}, {"A": 1.5, "B": 0.0}, 1.5)[0] == "A"   # breadth tracks depth
    assert call_idiomorphs({"A": 0.2, "B": 0.0}, {"A": 0.3, "B": 0.0}, 0.3)[0] == "low_depth"


def test_trace_signal_is_flagged_not_called():
    call, flags = call_idiomorphs({"A": 0.99, "B": 0.05}, {"A": 50.0, "B": 0.5}, 50.0)
    assert call == "A" and flags == ["trace_B"]


def test_low_coverage_sample_is_low_depth_not_none(tmp_path):
    fq = tmp_path / "r.fq"
    _write_fastq(fq, _reads(MAT2, 0.2, seed=6))
    assert type_reads(_panel(), [fq]).call == "low_depth"


def test_cli_writes_tsv(tmp_path):
    from MATPredict.__main__ import main
    for name, seq in (("m1.fa", MAT1), ("m2.fa", MAT2)):
        (tmp_path / name).write_text(f">{name}\n{seq}\n")
    fq = tmp_path / "s1.fq.gz"
    _write_fastq(fq, _reads(MAT1, 10, seed=1))
    out = tmp_path / "out.tsv"
    rc = main(["reads-type", "--idiomorph", f"MAT1-1={tmp_path / 'm1.fa'}", "--idiomorph", f"MAT1-2={tmp_path / 'm2.fa'}",
               "--k", "21", "--reads", str(fq), "--sample", "s1", "--out", str(out)])
    assert rc == 0
    rows = [l.split("\t") for l in out.read_text().splitlines()]
    assert rows[0][:2] == ["sample", "call"] and rows[1][:2] == ["s1", "MAT1-1"]


def test_cli_rejects_bad_idiomorph_spec(tmp_path):
    from MATPredict.__main__ import main
    assert main(["reads-type", "--idiomorph", "oops", "--reads", "x.fq", "--out", str(tmp_path / "o.tsv")]) == 1
