"""Fola loci: isolate 50a (genome assembly versus recruit-and-extend contigs), VSP-0947 two idiomorph contigs, extension by round.
usage: fig_fola_loci.py MAIN_RESULTS_DIR WORKTREE_ROOT OUT.png   (run with the pixi python)"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent))
from locusplot import Track, blast_blocks, read_fasta, align_tracks, draw, C
from gff_genes import gff_genes, locus_contig

M, W, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
ext = M / "2026-10-06_fola_reads_blastx/recruit_ext/50a"
r1 = M / "2026-10-06_fola_reads_blastx/recruit/50a"
genome = read_fasta("/bigdata/stajichlab/jstajich/projects/fola_reads/50a_spades/scaffolds.fasta")

def contig_track(label, fasta, report, gff, note=""):
    seqs = read_fasta(fasta)
    full, s, e, idio = locus_contig(report, "MAT1-2")
    t = Track(label, full, 1, len(seqs[full]), gff_genes(gff, full.split("_length")[0]), note=note or f"{len(seqs[full]) / 1000:.1f} kb contig")
    return t, {full: seqs[full]}

rounds = [("round 1 (blastx recruits)", r1 / "scaffolds.fasta", r1 / "detect/detection_report.yaml", r1 / "detect/detected_loci.gff3"),
          ("round 3", ext / "round3_scaffolds.fasta", ext / "detect_r3/detection_report.yaml", ext / "detect_r3/detected_loci.gff3"),
          ("round 5", ext / "round5_scaffolds.fasta", ext / "detect_r5/detection_report.yaml", ext / "detect_r5/detected_loci.gff3")]
ctr = [contig_track(*r) for r in rounds]
# genome track: window around the round-5 contig
gname = [k for k in genome if k.startswith("NODE_86_")][0]
gtrack_full = Track("whole-genome assembly\n(SPAdes, 50a)", gname, 114000, 133000, [], note="NODE_86")
gfull = locus_contig("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_50a/detect/detection_report.yaml", "MAT1-2")
gtrack_full.genes = [g for g in gff_genes("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_50a/detect/detected_loci.gff3", "NODE_86")]
tracks = [gtrack_full] + [t for t, _ in ctr]
seqs = [{gname: genome[gname]}] + [s for _, s in ctr]
links = []
for i in range(len(tracks) - 1):
    b = blast_blocks(seqs[i + 1], seqs[i], min_len=300, min_pid=95)       # lower track as query, upper as subject
    links.append([(x[4], x[5], x[1], x[2], x[6]) for x in b])
# links are (upper-track coords, lower-track coords); orient so tracks[i] is the upper
aligned = [[(l[0], l[1], l[2], l[3]) for l in ls] for ls in links]
# flip tracks whose blocks run in reverse relative to the previous track
for i, ls in enumerate(links):
    rev = sum(1 for l in ls if (l[0] > l[1]) != (l[2] > l[3])) > len(ls) / 2
    tracks[i + 1].flip = rev != tracks[i].flip if ls else False
align_tracks(tracks, aligned)

fig = plt.figure(figsize=(13.5, 10))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.3, wspace=0.34)
ax = fig.add_subplot(gs[0, :])
draw(ax, tracks, links, row_gap=1.35, link_alpha=0.4)
ax.set_title("a. Isolate 50a MAT1-2 locus: whole-genome assembly (top) and read-recruit-and-extend contigs; green ribbons = BLASTN >= 95% identity", loc="left", fontsize=9.5)

# b. VSP-0947: two idiomorph contigs against the Fola panel
ax = fig.add_subplot(gs[1, 0])
pan = W / "panels/fusarium_oxysporum_fola"
g947 = read_fasta("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/genomes/VSP-0947_AVITI_scaffolds.fasta")
rows = []
for idi, ref, tag in (("MAT1-1", "MAT1-1.fasta", "NODE_931_"), ("MAT1-2", "MAT1-2.fasta", "NODE_1613_")):
    (rn, rs), = read_fasta(pan / ref).items()
    cn = [k for k in g947 if k.startswith(tag)][0]
    rows.append((idi, rn, rs, cn, g947[cn]))
tr, lk = [], []
gff947 = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-06_fola_detect_all/runs/VSP-0947_AVITI_scaffolds/detected_loci.gff3"
for idi, rn, rs, cn, cs in rows:
    a = Track(f"{idi} reference\n({'VSP-0980' if idi == 'MAT1-1' else 'AT141'} panel)", rn, 1, len(rs), [], note=f"{len(rs) / 1000:.1f} kb")
    b = Track(f"VSP-0947 {cn.split('_length')[0]}\n(contig {len(cs) / 1000:.1f} kb, cov {float(cn.split('_cov_')[1]):.1f}x)", cn, 1, len(cs),
              gff_genes(gff947, cn.split("_length")[0]) if idi == "MAT1-1" else [])
    bl = blast_blocks({cn: cs}, {rn: rs}, min_len=100, min_pid=90)
    li = [(x[4], x[5], x[1], x[2], x[6]) for x in bl]
    b.flip = sum(1 for x in li if (x[0] > x[1]) != (x[2] > x[3])) > len(li) / 2
    align_tracks([a, b], [[(x[0], x[1], x[2], x[3]) for x in li]])
    tr += [a, b]; lk += [li, []]
lk = lk[:-1]
draw(ax, tr, lk, row_gap=1.3, show_scale=True)
ax.set_title("b. VSP-0947: each idiomorph on its own short contig, no flanks", loc="left", fontsize=9.5)
fig2_ax = fig.add_subplot(gs[1, 1])
# c. longest contig per round for the four pilot strains
for strain, colr in (("50a", C["MAT1-2"]), ("VSP-0980", C["MAT1-1"]), ("VSP-0931", "#009E73"), ("VSP-0947", C["both"])):
    f = M / f"2026-10-06_fola_reads_blastx/recruit_ext/{strain}/rounds.tsv"
    rr = [r for r in csv.DictReader(open(f), delimiter="\t") if r["round"].isdigit() and r["longest"]]
    fig2_ax.plot([int(r["round"]) for r in rr], [int(r["longest"]) / 1000 for r in rr], marker="o", color=colr, label=strain, lw=1.6)
fig2_ax.set_xlabel("recruit-and-extend round"); fig2_ax.set_ylabel("longest assembled contig (kb)"); fig2_ax.set_xticks([1, 2, 3, 4, 5])
fig2_ax.set_title("c. Contig length by round (pilot strains)", loc="left", fontsize=9.5); fig2_ax.legend(frameon=False, fontsize=8)
fig.savefig(OUT, dpi=150, bbox_inches="tight"); print("wrote", OUT)
