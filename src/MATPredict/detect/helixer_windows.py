"""helixer-lite on locus windows only.

The goal is better gene models for the loci the search already found, not an annotation of the genome.
Each detected locus is widened by a flank on both sides, overlapping windows on a contig are merged,
helixer-lite runs on the merged windows as one FASTA, and its window-relative coordinates are mapped back
to the genome (see tracks.from_helixer_gff3). The predictions are a track; they do not change calls.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from MATPredict.detect.tracks import Track, from_helixer_gff3

DEFAULT_FLANK = 5000


def merged_windows(spans, flank, contig_lengths):
    """spans: iterable of (contig, start, end), 1-based inclusive. Returns sorted, merged
    (contig, start, end) windows, each widened by `flank` and clipped to the contig."""
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for contig, start, end in spans:
        length = contig_lengths.get(contig)
        lo = max(1, start - flank)
        hi = end + flank if not length else min(length, end + flank)
        by_contig.setdefault(contig, []).append((lo, hi))
    out = []
    for contig, ivs in sorted(by_contig.items()):
        ivs.sort()
        cur_lo, cur_hi = ivs[0]
        for lo, hi in ivs[1:]:
            if lo <= cur_hi + 1:
                cur_hi = max(cur_hi, hi)
            else:
                out.append((contig, cur_lo, cur_hi))
                cur_lo, cur_hi = lo, hi
        out.append((contig, cur_lo, cur_hi))
    return out


def write_window_fasta(windows, genome_sequences: dict[str, str], path: Path) -> None:
    """One FASTA record per window, named contig:start-end so the GFF can be mapped back."""
    with open(path, "w") as fh:
        for contig, lo, hi in windows:
            seq = genome_sequences[contig][lo - 1:hi]
            fh.write(f">{contig}:{lo}-{hi}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def run_helixer_windows(spans, genome_sequences: dict[str, str], contig_lengths: dict[str, int],
                        *, flank: int = DEFAULT_FLANK, lineage: str = "fungi", cpus: int = 4,
                        binary: str | list[str] = "helixerlite") -> Track:
    """Run helixer-lite on the merged windows around `spans` and return its gene models as a track.
    `binary` may be a command list (e.g. a pixi or venv wrapper). Raises CalledProcessError on failure."""
    windows = merged_windows(spans, flank, contig_lengths)
    if not windows:
        return Track(source="helixer")
    cmd = [binary] if isinstance(binary, str) else list(binary)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        fasta, gff = tmp / "windows.fa", tmp / "helixer.gff3"
        write_window_fasta(windows, genome_sequences, fasta)
        subprocess.run(cmd + ["-f", str(fasta), "-o", str(gff), "-l", lineage, "-c", str(cpus)],
                       check=True, capture_output=True, text=True)
        return from_helixer_gff3(gff, genome_sequences, window_offsets=True)
