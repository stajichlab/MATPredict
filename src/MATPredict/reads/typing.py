"""Call the idiomorph(s) of a strain from reads by counting idiomorph-unique k-mers."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from MATPredict.reads.kmers import iter_fastq_reads
from MATPredict.reads.panel import Panel

BREADTH_FRACTION = 0.50   # observed breadth must reach this fraction of the breadth the depth predicts;
                          # exact k-mers lose ~20% of a reference to SNPs in a diverged allele (Fola MAT1-2: 0.8 of expected),
                          # scattered trace reads give <=0.15. Set from three strains, see analysis note.
MIN_REL_DEPTH = 0.10      # idiomorph depth relative to the shared single-copy k-mers; below = trace, not carried
MIN_SHARED_DEPTH = 1.0    # mean k-mer coverage of the shared control below which no call is made

@dataclass
class TypingResult:
    call: str
    breadth: dict[str, float]
    depth: dict[str, float]
    shared_depth: float
    reads_used: int
    flags: list[str] = field(default_factory=list)


def call_idiomorphs(breadth: dict[str, float], depth: dict[str, float], shared_depth: float) -> tuple[str, list[str]]:
    """Return (call, flags). An idiomorph is present when its depth is at least MIN_REL_DEPTH of the shared
    control and its breadth reaches BREADTH_FRACTION of the Poisson expectation 1 - exp(-depth), so low-depth
    data are not misread as partial coverage. Call: one idiomorph, 'both', 'none', or 'low_depth'."""
    if shared_depth < MIN_SHARED_DEPTH:
        return "low_depth", []
    present, flags = [], []
    for name, b in breadth.items():
        d = depth[name]
        expected = 1.0 - math.exp(-d)
        if d / shared_depth >= MIN_REL_DEPTH and expected > 0 and b >= BREADTH_FRACTION * expected:
            present.append(name)
        elif b > 0:
            flags.append(f"trace_{name}")
    call = present[0] if len(present) == 1 else "both" if present else "none"
    return call, flags


def type_reads(panel: Panel, fastqs: list[str | Path], max_reads: int | None = None) -> TypingResult:
    """Count panel k-mers in the reads of `fastqs` (all files pooled)."""
    table, classes = panel.lookup()
    counts = {name: [0] * len(kms) for name, kms in classes.items()}
    k = panel.k
    n_reads = 0
    for path in fastqs:
        for read in iter_fastq_reads(path):
            if max_reads is not None and n_reads >= max_reads:
                break
            n_reads += 1
            for i in range(len(read) - k + 1):
                hit = table.get(read[i:i + k])
                if hit is not None:
                    counts[hit[0]][hit[1]] += 1
        if max_reads is not None and n_reads >= max_reads:
            break
    breadth = {n: sum(1 for c in counts[n] if c) / len(counts[n]) for n in panel.unique}
    depth = {n: sum(counts[n]) / len(counts[n]) for n in panel.unique}
    shared = counts["shared"]
    shared_depth = sum(shared) / len(shared) if shared else 0.0
    call, flags = call_idiomorphs(breadth, depth, shared_depth)
    if call == "both" and min(depth.values()) / max(depth.values()) < 0.25:
        flags.append("idiomorph_depth_ratio_below_0.25")
    return TypingResult(call, breadth, depth, shared_depth, n_reads, flags)
