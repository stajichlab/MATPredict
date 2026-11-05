"""`shared_flank_check`: what do the two loci of a cross-family cassette pair share? (report only; calls are unchanged.)

Complementary idiomorphs called by two families at two positions are the two cassettes when their surroundings are unrelated. If the DNA around
the two loci matches they may be the a and alpha alleles of one heterozygous locus (*Wickerhamomyces edaphicus*: both ends of the 3.5 kb MTL A contig
match the MATsc alpha contig at 97%, one on each side of the idiomorph), but sharing alone proves nothing: real cassette systems share sequence
(the mirrored cassettes of *Saccharomycopsis schoenii* share 10-31 kb at 96-100%), a single short hit can be a repeat, and a collinear pair can be one
locus on two haplotypes. Curator ruling 2026-11-04 (#189, after an independent review): the check REPORTS the measurements and the role stays
`cassette_undetermined`; the role becomes `allelic_pair` only when shared sequence lies on BOTH sides of the idiomorph in BOTH loci and is not a
repeat (not found at more than two places in the genome).

The windows are each locus +/- WINDOW_BP, clipped to the contig. Measurements, in the returned dict:
  status         `shared_flank` (an alignment of at least MIN_ALIGNED_BP at MIN_IDENTITY % or more), `collinear` (alignments at COLLINEAR_IDENTITY % or more
                 cover at least COLLINEAR_FRACTION of the shorter window), `none`, or `not_run` (the genome or blastn was not available)
  aligned_bp, best_identity   the qualifying alignments
  both_sides     true when, in each locus, a qualifying alignment lies left of the idiomorph's genes and another right of them
  repeat         true when a flank segment has a copy at a third place in the genome; null when it could not be checked
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


def both_sides(hsps: list[tuple], off_a: int, off_b: int, core_a: tuple[int, int], core_b: tuple[int, int]) -> bool:
    """True when, in each locus, a qualifying alignment lies left of the locus's core genes (`core_*` = start, end in contig coordinates) and
    another right of them; an alignment counts as left or right by its midpoint, so a hit that runs a base or two into a gene still counts.
    `off_*` is the window start minus 1, turning window into contig coordinates."""
    shared = [h for h in hsps if h[0] >= MIN_IDENTITY and h[1] >= MIN_ALIGNED_BP]

    def sides(spans, core):
        lo, hi = core
        return any((b + e) / 2 < lo for b, e in spans), any((b + e) / 2 > hi for b, e in spans)

    spans_a = [(min(h[2], h[3]) + off_a, max(h[2], h[3]) + off_a) for h in shared]
    spans_b = [(min(h[4], h[5]) + off_b, max(h[4], h[5]) + off_b) for h in shared]
    return all(sides(spans_a, core_a)) and all(sides(spans_b, core_b))


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


def _core_span(call) -> tuple[int, int]:
    """Extent of the call's core_MAT genes on its contig; the call's own span when it lists none."""
    own = [g for g in (getattr(call, "gene_evidence", None) or []) if g.role == "core_MAT" and g.contig == call.contig and not getattr(g, "lost_to", None)]
    return (min(g.start for g in own), max(g.end for g in own)) if own else (call.start, call.end)


def _copies(exe: str, segments: list[str], genome_fasta: Path, tmp: Path) -> int | None:
    """Most places in the genome (alignments at 90% identity over 80% of the segment) that any flank segment has; None when not checkable.
    A true pair gives 2 (the two loci); more means a repeat."""
    if str(genome_fasta).endswith(".gz") or not Path(genome_fasta).exists():
        return None
    q = tmp / "seg.fa"
    q.write_text("".join(f">s{i}\n{x}\n" for i, x in enumerate(segments)))
    try:
        run = subprocess.run([exe, "-task", "megablast", "-query", str(q), "-subject", str(genome_fasta), "-evalue", "1e-20",
                              "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send"], capture_output=True, text=True, timeout=600)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if run.returncode != 0:
        return None
    qlen = {f"s{i}": len(x) for i, x in enumerate(segments)}
    places: dict[str, list[tuple[str, int]]] = {}
    for ln in run.stdout.splitlines():
        f = ln.split("\t")
        if len(f) < 8 or float(f[2]) < 90.0 or int(f[3]) < 0.8 * qlen[f[0]]:
            continue
        places.setdefault(f[0], []).append((f[1], (int(f[6]) + int(f[7])) // 2 // 1000))
    return max((len(set(v)) for v in places.values()), default=0)


def shared_flank_check(genome_fasta: Path, a, b, *, blastn: str | None = None) -> dict:
    """The measurements for two calls `a` and `b`; status `not_run` when the genome cannot be read or blastn is missing."""
    exe = blastn or shutil.which("blastn")
    if exe is None:
        return {"status": "not_run", "reason": "blastn not found"}
    seqs = _read_contigs(Path(genome_fasta), {a.contig, b.contig})
    if a.contig not in seqs or b.contig not in seqs:
        return {"status": "not_run", "reason": "genome contigs not readable"}
    wa, lo_a = _window(seqs[a.contig], a.start, a.end)
    wb, lo_b = _window(seqs[b.contig], b.start, b.end)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "a.fa").write_text(f">a\n{wa}\n")
        (tmp / "b.fa").write_text(f">b\n{wb}\n")
        try:
            run = subprocess.run([exe, "-task", "megablast", "-query", str(tmp / "a.fa"), "-subject", str(tmp / "b.fa"), "-evalue", "1e-20",
                                  "-outfmt", "6 pident length qstart qend sstart send"], capture_output=True, text=True, timeout=300)
        except (OSError, subprocess.TimeoutExpired) as exc:
            return {"status": "not_run", "reason": str(exc)[:80]}
        if run.returncode != 0:
            return {"status": "not_run", "reason": run.stderr.strip()[:80]}
        hsps = [(float(f[0]), int(f[1]), int(f[2]), int(f[3]), int(f[4]), int(f[5])) for f in (ln.split("\t") for ln in run.stdout.splitlines()) if len(f) >= 6]
        out = classify(hsps, len(wa), len(wb))
        out["both_sides"] = out["status"] in ("shared_flank", "collinear") and both_sides(hsps, lo_a - 1, lo_b - 1, _core_span(a), _core_span(b))
        out["repeat"] = None
        if out["both_sides"]:
            segs = [wa[min(h[2], h[3]) - 1:max(h[2], h[3])] for h in hsps if h[0] >= MIN_IDENTITY and h[1] >= MIN_ALIGNED_BP]
            copies = _copies(exe, segs, Path(genome_fasta), tmp)
            out["repeat"] = None if copies is None else copies > 2
    return out
