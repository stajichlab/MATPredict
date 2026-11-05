"""`shared_flank_check`: do the two loci of a cross-family cassette pair share sequence? (report only; calls are unchanged.)

Complementary idiomorphs called by two families at two positions are the two cassettes only when their surroundings are unrelated. If
the DNA around the two loci matches, they are the a and alpha alleles of one heterozygous locus (*Wickerhamomyces edaphicus*: both ends of the
2.4 kb MTL A contig match the MATsc alpha contig at 97%, either side of the idiomorph), or one locus seen on two haplotypes
(*Saccharomycopsis fibuligera*: collinear contigs at 92-98%). Curator ruling 2026-11-04 (#189), after an independent review of four exemplars.

The windows are each locus +/- WINDOW_BP, clipped to the contig. Statuses:
  shared_flank   a local alignment of at least MIN_ALIGNED_BP at MIN_IDENTITY % or more between the two windows
  collinear      alignments at COLLINEAR_IDENTITY % or more cover at least COLLINEAR_FRACTION of the shorter window (one locus on two haplotypes)
  none           no such alignment: the loci are unrelated, so a candidate second cassette
  not_run        the genome or blastn was not available
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path

from MATPredict.detect.caax import _read_contigs

WINDOW_BP = 20000
MIN_IDENTITY = 95.0
MIN_ALIGNED_BP = 100
COLLINEAR_IDENTITY = 90.0
COLLINEAR_FRACTION = 0.5
STATUSES = ("shared_flank", "collinear", "none", "not_run")


def classify(hsps: list[tuple], len_a: int, len_b: int) -> dict:
    """Status from blastn alignments, each (identity %, aligned bp, q_start, q_end, s_start, s_end) in window coordinates."""
    shared = [h for h in hsps if h[0] >= MIN_IDENTITY and h[1] >= MIN_ALIGNED_BP]
    covered = _union_bp([(min(h[2], h[3]), max(h[2], h[3])) for h in hsps if h[0] >= COLLINEAR_IDENTITY])
    if covered >= COLLINEAR_FRACTION * min(len_a, len_b):
        status = "collinear"
    elif shared:
        status = "shared_flank"
    else:
        status = "none"
    best = max(shared or hsps or [(None, 0)], key=lambda h: h[1])
    return {"status": status, "aligned_bp": sum(h[1] for h in shared), "best_identity": best[0], "window_bp": WINDOW_BP}


def _union_bp(spans: list[tuple[int, int]]) -> int:
    total, end = 0, None
    for a, b in sorted(spans):
        if end is None or a > end:
            total += b - a + 1
            end = b
        elif b > end:
            total += b - end
            end = b
    return total


def _window(seq: str, start: int, end: int) -> tuple[str, int]:
    lo = max(1, start - WINDOW_BP)
    hi = min(len(seq), end + WINDOW_BP)
    return seq[lo - 1:hi], lo


def shared_flank_check(genome_fasta: Path, a, b, *, blastn: str | None = None) -> dict:
    """The check for two calls `a` and `b` (contig, start, end); status `not_run` when the genome cannot be read or blastn is missing."""
    exe = blastn or shutil.which("blastn")
    if exe is None:
        return {"status": "not_run", "reason": "blastn not found"}
    seqs = _read_contigs(Path(genome_fasta), {a.contig, b.contig})
    if a.contig not in seqs or b.contig not in seqs:
        return {"status": "not_run", "reason": "genome contigs not readable"}
    wa, _ = _window(seqs[a.contig], a.start, a.end)
    wb, _ = _window(seqs[b.contig], b.start, b.end)
    with tempfile.TemporaryDirectory() as tmp:
        qa, qb = Path(tmp) / "a.fa", Path(tmp) / "b.fa"
        qa.write_text(f">a\n{wa}\n")
        qb.write_text(f">b\n{wb}\n")
        try:
            run = subprocess.run([exe, "-task", "megablast", "-query", str(qa), "-subject", str(qb), "-evalue", "1e-20",
                                  "-outfmt", "6 pident length qstart qend sstart send"], capture_output=True, text=True, timeout=300)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return {"status": "not_run", "reason": str(exc)[:80]}
    if run.returncode != 0:
        return {"status": "not_run", "reason": run.stderr.strip()[:80]}
    hsps = [(float(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5])) for f in (ln.split("\t") for ln in run.stdout.splitlines()) if len(f) >= 6]
    return classify(hsps, len(wa), len(wb))
