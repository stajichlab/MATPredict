"""Small locus-map helpers (matplotlib only): gene arrows on tracks, BLAST-block ribbons between tracks."""
from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Okabe-Ito colour-blind-safe palette
C = {"MAT1-1": "#0072B2", "MAT1-2": "#D55E00", "both": "#CC79A7", "none": "#999999", "flank": "#BBBBBB", "ink": "#222222"}
GENE_COLOR = {"COX13": C["flank"], "APN2": C["flank"], "SLA2": C["flank"], "MAT1-1-1": C["MAT1-1"], "MAT1-1-2": C["MAT1-1"],
              "MAT1-1-3": C["MAT1-1"], "MAT1-2-1": C["MAT1-2"], "MAT1-2-4": "#E9A97B", "remnant": "#8FB8D8"}


@dataclass
class Track:
    label: str
    contig: str
    x0: int                      # window start (bp, contig coordinates)
    x1: int                      # window end
    genes: list = field(default_factory=list)   # (start, end, strand, name)
    flip: bool = False
    shift: float = 0.0           # kb, set by align_tracks
    note: str = ""
    absent: bool = False         # draw a dashed placeholder (idiomorph not found in this data type)

    def x(self, pos: float) -> float:
        base = (self.x1 - pos) if self.flip else (pos - self.x0)
        return base / 1000.0 + self.shift


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"):
            name = line[1:].split()[0]; out[name] = []
        elif name:
            out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def blast_blocks(query_seq: dict, subject_seq: dict, min_len=150, min_pid=90.0):
    """BLASTN of query contigs against subject contigs. Returns [(q, qs, qe, s, ss, se, pid)]; subject coords ordered ss<=se with 'strand'."""
    with tempfile.TemporaryDirectory() as d:
        q, s = Path(d) / "q.fa", Path(d) / "s.fa"
        q.write_text("".join(f">{k}\n{v}\n" for k, v in query_seq.items()))
        s.write_text("".join(f">{k}\n{v}\n" for k, v in subject_seq.items()))
        subprocess.run(["makeblastdb", "-in", str(s), "-dbtype", "nucl", "-out", str(Path(d) / "db")], check=True, capture_output=True)
        out = subprocess.run(["blastn", "-query", str(q), "-db", str(Path(d) / "db"), "-evalue", "1e-20", "-perc_identity", str(min_pid),
                              "-outfmt", "6 qseqid qstart qend sseqid sstart send pident length", "-max_target_seqs", "50"],
                             check=True, capture_output=True, text=True).stdout
    blocks = []
    for line in out.splitlines():
        a, qs, qe, b, ss, se, pid, ln = line.split("\t")
        if int(ln) >= min_len:
            blocks.append((a, int(qs), int(qe), b, int(ss), int(se), float(pid)))
    return blocks


def align_tracks(tracks, links):
    """links[i] are blocks between tracks[i] and tracks[i+1] as (qs, qe, ss, se). Shift track i+1 so linked ends line up."""
    for i, blocks in enumerate(links):
        if not blocks:
            continue
        a, b = tracks[i], tracks[i + 1]
        num = sum(((a.x((qs + qe) / 2) - (b.x((ss + se) / 2) - b.shift)) * (abs(qe - qs) + 1)) for qs, qe, ss, se in blocks)
        den = sum(abs(qe - qs) + 1 for qs, qe, ss, se in blocks)
        b.shift = num / den


def draw(ax, tracks, links, row_gap=1.0, gene_h=0.22, show_scale=True, label_genes=True, link_alpha=0.35):
    ys = {i: -i * row_gap for i in range(len(tracks))}
    xmin, xmax = 1e9, -1e9
    for i, t in enumerate(tracks):
        y = ys[i]
        xa, xb = sorted((t.x(t.x0), t.x(t.x1)))
        ax.plot([xa, xb], [y, y], color=C["ink"] if not t.absent else "#AAAAAA", lw=1.0, ls="-" if not t.absent else (0, (3, 3)), zorder=1)
        xmin, xmax = min(xmin, xa), max(xmax, xb)
        for (s, e, strand, name) in t.genes:
            if e < t.x0 or s > t.x1:
                continue
            s, e = max(s, t.x0), min(e, t.x1)
            fwd = (strand == "+") != t.flip
            xs, xe = sorted((t.x(s), t.x(e)))
            head = min(0.35, (xe - xs) * 0.4)
            col = GENE_COLOR.get(name, "#DDDDDD")
            if fwd:
                pts = [(xs, y - gene_h / 2), (xe - head, y - gene_h / 2), (xe, y), (xe - head, y + gene_h / 2), (xs, y + gene_h / 2)]
            else:
                pts = [(xe, y - gene_h / 2), (xs + head, y - gene_h / 2), (xs, y), (xs + head, y + gene_h / 2), (xe, y + gene_h / 2)]
            ax.add_patch(Polygon(pts, closed=True, facecolor=col, edgecolor=C["ink"], lw=0.6, zorder=3,
                                 hatch="////" if name == "remnant" else None))
            if label_genes and (xe - xs) > 0.45:
                ax.text((xs + xe) / 2, y + gene_h / 2 + 0.05, name, ha="center", va="bottom", fontsize=6.5, color=C["ink"])
        ax.text(xmin - 0.4, y, t.label + (f"\n{t.note}" if t.note else ""), ha="right", va="center", fontsize=7.5)
    for i, blocks in enumerate(links):
        a, b = tracks[i], tracks[i + 1]
        for (qs, qe, ss, se, pid) in blocks:
            ya, yb = ys[i] - gene_h / 2 - 0.02, ys[i + 1] + gene_h / 2 + 0.02
            pts = [(a.x(qs), ya), (a.x(qe), ya), (b.x(se), yb), (b.x(ss), yb)]
            ax.add_patch(Polygon(pts, closed=True, facecolor="#7FB069" if pid >= 98 else "#E8C547", alpha=link_alpha,
                                 edgecolor="none", zorder=0))
    if show_scale:
        y0 = min(ys.values()) - 0.55
        ax.plot([xmin, xmin + 2], [y0, y0], color=C["ink"], lw=1.5)
        ax.text(xmin + 1, y0 - 0.08, "2 kb", ha="center", va="top", fontsize=7)
    ax.set_xlim(xmin - 6.0, xmax + 0.5)
    ax.set_ylim(min(ys.values()) - 0.9, 0.7)
    ax.axis("off")
