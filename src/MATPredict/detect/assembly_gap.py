"""`assembly_gap_at_locus`: the locus position is an assembly gap.

Curator's ruling 2026-09-26 (results/2026-09-26_cauris_uncalled/NOTE.md). Of
12 uncalled C. auris genomes, 10 are reference-consensus assemblies (CLC,
BWA-MEM) with 8.1-8.9 kb of N runs exactly over the reference idiomorph. The
locus was never absent; the assembly could not resolve it, most likely
because the isolate carries the other idiomorph and its reads cannot map
across non-homologous sequence. Reporting such a genome as plain
"not detected" reads as absence. This module finds the case and names it.

The rule, applied only to families with no reported call:

* An N block is one or more N runs merged across islands of sequence no
  longer than `GAP_MERGE_BP`, holding at least `MIN_GAP_N_BP` N bases.
  Merging is needed because the real gaps are not one run: the SCO genomes
  hold five runs split by 16-56 bp islands, the England genomes three runs
  split by 21 and 143 bp.
* The block is AT THE LOCUS when a hit of one of the family's flanking genes
  (any `flanking_*` role) at `ANCHOR_MIN_IDENTITY` or more lies inside it or
  within `ANCHOR_WINDOW_BP` of either end. Core hits do not anchor: the
  ruling anchors on flanks, and a lone short core fragment is the commonest
  background hit.

Measured on the 16 C. auris genomes of results/2026-09-26_gap_validation
(10 reference-consensus, 1 de novo break, 5 called controls). Without the
identity floor all 10 gap genomes were flagged, but 8 of them only at
117-353 bp N blocks beside flank PARALOG hits at 22.9-35.7% -- the wrong
place. The only anchors at the real MTL gap were islands of sequence inside
it, hitting PIK1 at 91.5% and PAP1 at 100%. With the floor, 2 genomes are
flagged, both at the true gap, and nothing else.

Limit, stated plainly: the C. auris MTL flanks in the roster (PAP1, OBP1,
PIK1) lie INSIDE the idiomorph, so when the whole idiomorph is N there is
usually nothing left to anchor on. A conserved gene outside the idiomorph,
curated as a flank, is what would let the other 8 be reported. In a
lineage whose true flank orthologs hit below the floor, no gap is reported:
the rule misses rather than misplaces.

It reports; it never calls. A gap report carries no idiomorph and changes no
confidence, class or not-detected reason.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from pathlib import Path

from MATPredict.detect.family_registry import FamilyKey

logger = logging.getLogger(__name__)

#: Minimum N bases in a block. 100 is the NCBI/AGP convention for a gap of
#: UNKNOWN size, so a 100-N run can hide a whole idiomorph; shorter runs are
#: normally sized gaps of a few dozen bases, too small to hold a MAT gene.
MIN_GAP_N_BP = 100

#: Longest island of real sequence that still joins two N runs into one block.
#: The measured islands inside the C. auris gaps are 16-143 bp.
GAP_MERGE_BP = 500

#: How far outside a block a flank hit may lie and still anchor it. A gene
#: that borders an idiomorph sits within a few kb of its edge.
ANCHOR_WINDOW_BP = 5_000

#: Minimum percent identity of an anchoring flank hit. Separates same-locus
#: flank hits (91.5-100% measured; the other idiomorph's flank allele hits at
#: 51-68%) from genome-wide flank paralogs (22.9-35.7% measured).
ANCHOR_MIN_IDENTITY = 50.0

_N_RUN = re.compile(r"[Nn]+")


@dataclass(frozen=True)
class AssemblyGapAtLocus:
    """One N block at a flank-anchored locus position of an uncalled family."""

    family_key: FamilyKey
    contig: str
    start: int  # 1-based, inclusive: first N of the block
    end: int  # 1-based, inclusive: last N of the block
    n_bases: int
    anchors: list[str] = field(default_factory=list)


def n_blocks(
    sequence: str, min_n: int = MIN_GAP_N_BP, merge_bp: int = GAP_MERGE_BP,
) -> list[tuple[int, int, int]]:
    """(start, end, n_bases) of every N block in `sequence`, 1-based inclusive."""
    blocks: list[list[int]] = []
    for match in _N_RUN.finditer(sequence):
        start, end = match.start() + 1, match.end()
        if blocks and start - blocks[-1][1] - 1 <= merge_bp:
            blocks[-1][1] = end
            blocks[-1][2] += end - start + 1
        else:
            blocks.append([start, end, end - start + 1])
    return [(s, e, n) for s, e, n in blocks if n >= min_n]


def _contig_sequences(genome_fasta: Path, wanted: set[str]) -> dict[str, str]:
    """The sequences of `wanted` contigs, or {} when the FASTA cannot be read."""
    from MATPredict.detect.benchmark import _open_fasta_text

    sequences: dict[str, str] = {}
    try:
        with _open_fasta_text(genome_fasta) as handle:
            name, chunks = None, []
            for line in handle:
                if line.startswith(">"):
                    if name in wanted:
                        sequences[name] = "".join(chunks)
                    name, chunks = (line[1:].split() or [""])[0], []
                elif name in wanted:
                    chunks.append(line.strip())
            if name in wanted:
                sequences[name] = "".join(chunks)
    except OSError:
        logger.warning("assembly_gap: could not read genome FASTA %s", genome_fasta)
        return {}
    return sequences


def find_gaps_at_locus(
    hits_by_family: dict[FamilyKey, list], genome_fasta: Path,
    window_bp: int = ANCHOR_WINDOW_BP, min_identity: float = ANCHOR_MIN_IDENTITY,
) -> list[AssemblyGapAtLocus]:
    """Every flank-anchored N block for the families in `hits_by_family`.

    The caller passes only families with no reported call, each with its own
    hits. The genome is read once, and only the contigs a flank hit names.
    """
    flanks = {
        key: [h for h in hits
              if h.role.startswith("flanking") and h.identity >= min_identity]
        for key, hits in hits_by_family.items()
    }
    wanted = {h.contig for hits in flanks.values() for h in hits}
    if not wanted:
        return []
    sequences = _contig_sequences(genome_fasta, wanted)
    blocks = {contig: n_blocks(seq) for contig, seq in sequences.items()}
    gaps: list[AssemblyGapAtLocus] = []
    for key, hits in flanks.items():
        for contig in sorted({h.contig for h in hits}):
            for start, end, n in blocks.get(contig, []):
                anchors = sorted({
                    h.gene_name for h in hits
                    if h.contig == contig
                    and h.end >= start - window_bp and h.start <= end + window_bp
                })
                if anchors:
                    gaps.append(AssemblyGapAtLocus(key, contig, start, end, n, anchors))
    return gaps
