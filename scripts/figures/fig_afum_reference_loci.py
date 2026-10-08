"""Figure: the two curated A. fumigatus locus segments aligned (BLASTN blocks); genes from the database GFF3.
usage: fig_afum_reference_loci.py WORKTREE_ROOT DB_EUROTIALES_DIR OUT.png  (run with the pixi python so blastn is on PATH)"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import matplotlib.pyplot as plt
from locusplot import Track, blast_blocks, read_fasta, align_tracks, draw, C

ROOT, DB, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
pan = ROOT / "panels/aspergillus_fumigatus"


def segment_genes(rec, rename=None):
    off = None
    for l in open(DB / rec / "locus.gff3"):
        if l.startswith("##sequence-region"):
            off = int(l.split()[2])
    out = []
    for l in open(DB / rec / "locus.gff3"):
        f = l.rstrip("\n").split("\t")
        if len(f) > 8 and f[2] == "gene":
            name = dict(kv.split("=") for kv in f[8].split(";"))["Name"]
            if rename and name in rename:
                name = rename[name]
            out.append((int(f[3]) - off + 1, int(f[4]) - off + 1, f[6], name))
    return out


(na, sa), = read_fasta(pan / "MAT1-1.fasta").items()
(nf, sf), = read_fasta(pan / "MAT1-2.fasta").items()
ga = segment_genes("746128_a1163_MAT_MAT1-1", {"MAT1-2-1": "remnant"})
gf = segment_genes("746128_af293_MAT_MAT1-2")
blocks = blast_blocks({"A1163": sa}, {"Af293": sf}, min_len=100, min_pid=85)
rev = sum(1 for b in blocks if b[4] > b[5]) > len(blocks) / 2
ta = Track("A1163 MAT1-1", "A1163", 1, len(sa), ga, note="database record")
tf = Track("Af293 MAT1-2", "Af293", 1, len(sf), gf, flip=rev, note="database record")
links = [[(b[1], b[2], b[4], b[5], b[6]) for b in blocks]]   # keep strand: ss > se means reverse
align_tracks([ta, tf], [[(l[0], l[1], l[2], l[3]) for l in links[0]]])
fig, ax = plt.subplots(figsize=(9, 2.6))
draw(ax, [ta, tf], links, row_gap=1.25)
ax.set_title("A. fumigatus MAT1-1 (A1163) and MAT1-2 (Af293) loci: shared flanks, idiomorph-specific core and the MAT1-2-1 remnant", fontsize=8.5, loc="left")
fig.tight_layout(); fig.savefig(OUT, dpi=160); print("wrote", OUT, "blocks", len(blocks), "flipped", rev)
