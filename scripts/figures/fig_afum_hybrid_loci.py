"""Hybrid and weak-signal A. fumigatus strains: MAT1-1 and MAT1-2 loci as reference, reads-assembled contig, whole-genome assembly.
usage: fig_afum_hybrid_loci.py MAIN_RESULTS_DIR WORKTREE_ROOT DB_EUROTIALES OUT.png   (run with the pixi python)"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from locusplot import Track, blast_blocks, read_fasta, align_tracks, draw, C
from protein_genes import genes_on

M, W, DB, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
R = M / "2026-10-07_afum_reads"
GEN = Path("/bigdata/stajichlab/shared/projects/Afumigatus_pangenome/scaffolded/genomes")
pan = W / "panels/aspergillus_fumigatus"
(refn1, ref1), = read_fasta(pan / "MAT1-1.fasta").items()
(refn2, ref2), = read_fasta(pan / "MAT1-2.fasta").items()
spec = read_fasta(pan / "specific_regions.fasta")
specs = {"MAT1-1": [v for k, v in spec.items() if k.startswith("MAT1-1")][0], "MAT1-2": [v for k, v in spec.items() if k.startswith("MAT1-2")][0]}
refs = {"MAT1-1": ("A1163 MAT1-1", ref1), "MAT1-2": ("Af293 MAT1-2", ref2)}
truth = {r["strain"]: r for r in csv.DictReader(open(R / "blast/assembly_truth.tsv"), delimiter="\t")}
STRAINS = [("IFM_59359", "putative hybrid (Lofgren et al.)"), ("IFM_61407", "putative hybrid (Lofgren et al.)"),
           ("DMC_AF100-1_3", "putative hybrid; no assembly here"), ("AF100-12_2", "assembly MAT1-2 only; reads show MAT1-1 (marginal)"),
           ("F18149-Manchester", "assembly both; reads mostly MAT1-1")]


def best_contig(query_seq, subject, min_len=150, min_pid=90.0):
    bl = blast_blocks({"q": query_seq}, subject, min_len=min_len, min_pid=min_pid)
    tot = {}
    for b in bl:
        tot[b[3]] = tot.get(b[3], 0) + abs(b[5] - b[4]) + 1
    if not tot:
        return None, []
    c = max(tot, key=tot.get)
    return c, [b for b in bl if b[3] == c]


def column(strain, idi, genome_fa):
    label, ref = refs[idi]
    tracks, seqs = [], []
    rt = Track(label, "ref", 1, len(ref), genes_on({"ref": ref}, DB)["ref"], note="curated record")
    tracks.append(rt); seqs.append({"ref": ref})
    # reads-assembled contig (round 4)
    rounds = sorted(int(p.name[5]) for p in (R / "recruit" / strain).glob("round?_scaffolds.fasta"))
    nround = rounds[-1] if rounds else 0
    f = R / "recruit" / strain / f"round{nround}_scaffolds.fasta"
    rd = read_fasta(f) if rounds else {}
    c, _ = best_contig(specs[idi], rd, min_len=100, min_pid=90) if rd else (None, [])
    if c:
        t = Track(f"{strain}\nreads, round {nround}", c, 1, len(rd[c]), genes_on({c: rd[c]}, DB).get(c, []), note=f"{len(rd[c]) / 1000:.1f} kb")
    else:
        t = Track(f"{strain}\nreads, round {nround}", "none", 1, len(ref), [], absent=True, note="no contig")
    tracks.append(t); seqs.append({c: rd[c]} if c else {})
    # whole-genome assembly
    st = truth.get(strain)
    gfa = GEN / f"{strain}.sorted.fasta"
    if st and st[f"{idi}_state"] != "absent" and gfa.exists():
        g = read_fasta(gfa)
        c2, bl = best_contig(specs[idi], g, min_len=100, min_pid=90)
        if c2:
            lo, hi = min(min(b[4], b[5]) for b in bl) - 4500, max(max(b[4], b[5]) for b in bl) + 4500
            lo, hi = max(1, lo), min(len(g[c2]), hi)
            win = g[c2][lo - 1:hi]
            t = Track(f"{strain}\nassembly", c2, 1, len(win), genes_on({c2: win}, DB).get(c2, []), note=f"{c2}:{lo}-{hi}")
            tracks.append(t); seqs.append({c2: win})
        else:
            tracks.append(Track(f"{strain}\nassembly", "none", 1, len(ref), [], absent=True, note="no region")); seqs.append({})
    else:
        why = "no assembly" if not (st and gfa.exists()) else "region absent (BLAST)"
        tracks.append(Track(f"{strain}\nassembly", "none", 1, len(ref), [], absent=True, note=why)); seqs.append({})
    links = []
    for i in range(len(tracks) - 1):
        if tracks[i].absent or tracks[i + 1].absent or not seqs[i] or not seqs[i + 1]:
            links.append([]); continue
        bl = blast_blocks(seqs[i + 1], seqs[i], min_len=150, min_pid=92)
        links.append([(b[4], b[5], b[1], b[2], b[6]) for b in bl])
    for i, ls in enumerate(links):
        if ls:
            tracks[i + 1].flip = (sum(1 for l in ls if (l[0] > l[1]) != (l[2] > l[3])) > len(ls) / 2) != tracks[i].flip
    align_tracks(tracks, [[(l[0], l[1], l[2], l[3]) for l in ls] for ls in links])
    return tracks, links


summary = []
fig, axes = plt.subplots(len(STRAINS), 2, figsize=(15, 3.1 * len(STRAINS)))
for ri, (strain, note) in enumerate(STRAINS):
    for ci, idi in enumerate(("MAT1-1", "MAT1-2")):
        ax = axes[ri][ci]
        tracks, links = column(strain, idi, None)
        for t in tracks[1:]:
            kind = "reads" if "reads" in t.label else "assembly"
            summary.append((strain, idi, kind, t.label.split("\n")[1] if "\n" in t.label else "", "absent" if t.absent else t.contig,
                            "" if t.absent else t.x1 - t.x0 + 1, ",".join(sorted(g[3] for g in t.genes)), t.note))
        draw(ax, tracks, links, row_gap=1.0, show_scale=(ci == 1), link_alpha=0.4)
        ax.set_title(f"{strain}: {idi}" + (f"   ({note})" if ci == 0 else ""), loc="left", fontsize=8.5)
fig.tight_layout(); fig.savefig(OUT, dpi=130, bbox_inches="tight"); print("wrote", OUT)
with open(Path(OUT).with_suffix(".tsv"), "w") as o:
    o.write("strain\tidiomorph\tsource\tlabel\tcontig\twindow_bp\tgenes_found\tnote\n")
    for r in summary:
        o.write("\t".join(str(x) for x in r) + "\n")
