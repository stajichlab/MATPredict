"""Check and write the P. fabianii MAT-locus models: CDS from genome coordinates, translation, start/stop/splice-site checks.
Models are AUGUSTUS predictions (saccharomyces_cerevisiae_S288C for MATa1 and MATalpha1, pichia_stipitis for MATa2, any model for DIC1/SLA2)
selected by agreement across species models and homology; no RNA-seq exists for the species. Usage: build_models.py GENOME_DIR"""
import subprocess, sys, re
D = sys.argv[1]
G = {"CP048715.1": "GCA_022641835.1_ASM2264183v1", "BCGI01000001.1": "GCA_001599195.1_JCM_3601_assembly_v001",
     "CM115473.1": "GCA_050575385.1_ASM5057538v1", "MPUK01000005.1": "GCA_001983305.1_ASM198330v1"}
# (record, contig, gene, strand, exons in genomic order, note)
M = [
 ("MATa:JOY008", "CP048715.1", "dic1", "-", [(628452, 629375)]),
 ("MATa:JOY008", "CP048715.1", "MATa2", "-", [(629853, 630419), (630492, 630863)]),
 ("MATa:JOY008", "CP048715.1", "MATa1", "+", [(631437, 631737), (631835, 632133)]),
 ("MATa:JOY008", "CP048715.1", "sla2", "+", [(641746, 644868)]),
 ("MATa:JCM3601", "BCGI01000001.1", "dic1", "-", [(641150, 642073)]),
 ("MATa:JCM3601", "BCGI01000001.1", "MATa2", "-", [(642551, 643117), (643190, 643561)]),
 ("MATa:JCM3601", "BCGI01000001.1", "MATa1", "+", [(644135, 644435), (644533, 644831)]),
 ("MATa:JCM3601", "BCGI01000001.1", "sla2", "+", [(654441, 657563)]),
 ("MATalpha:TBRC4498", "CM115473.1", "dic1", "-", [(631649, 632572)]),
 ("MATalpha:TBRC4498", "CM115473.1", "MATalpha1", "+", [(637540, 638454)]),
 ("MATalpha:TBRC4498", "CM115473.1", "sla2", "+", [(648035, 651157)]),
 ("MATalpha:strain65", "MPUK01000005.1", "dic1", "-", [(611410, 612333)]),
 ("MATalpha:strain65", "MPUK01000005.1", "MATalpha1", "+", [(617300, 618214)]),
 ("MATalpha:strain65", "MPUK01000005.1", "sla2", "+", [(627795, 630917)]),
]
b = "TCAG"; aa = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
C = {x + y + z: aa[16 * i + 4 * j + k] for i, x in enumerate(b) for j, y in enumerate(b) for k, z in enumerate(b)}
def rc(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))
def tr(s): return "".join(C.get(s[i:i + 3], "X") for i in range(0, len(s) - 2, 3))
def get(contig, a, b_):
    return "".join(subprocess.check_output(["samtools", "faidx", f"{D}/{G[contig]}.fa", f"{contig}:{a}-{b_}"], text=True).split("\n")[1:]).upper()
out = open("data/models.tsv", "w"); fa = open("data/models_proteins.faa", "w")
out.write("record\tcontig\tgene\tstrand\texons\tcds_nt\taa\tstart\tstop_in_frame_end\tintron_signals\n")
for rec, c, g, st, ex in M:
    seqs = [get(c, a, b_) for a, b_ in ex]
    if st == "-": seqs = [rc(s) for s in seqs[::-1]]; ex_t = ex[::-1]
    else: ex_t = ex
    cds = "".join(seqs)
    # next three bases after the last exon end (stop codon if not included)
    last = ex[0][0] if st == "-" else ex[-1][1]
    nxt = rc(get(c, last - 3, last - 1)) if st == "-" else get(c, last + 1, last + 3)
    p = tr(cds)
    has_stop = p.endswith("*")
    stop_ext = (not has_stop) and tr(nxt) == "*"
    sig = []
    for i in range(len(ex_t) - 1):  # intron between consecutive exons in transcript order
        if st == "+": a_, b_ = ex_t[i][1] + 1, ex_t[i + 1][0] - 1; intr = get(c, a_, b_)
        else: a_, b_ = ex_t[i + 1][1] + 1, ex_t[i][0] - 1; intr = rc(get(c, a_, b_))
        sig.append(f"{a_}-{b_}({len(intr)}bp):{intr[:6]}..{intr[-6:]}")
    internal = p.rstrip("*")[:-0 or None].count("*")
    out.write("\t".join([rec, c, g, st, ";".join(f"{x}-{y}" for x, y in ex_t if True), str(len(cds)), str(len(p.rstrip('*'))), p[:1] + ("(ATG)" if cds[:3] == "ATG" else "(not ATG)"),
              "stop_in_cds" if has_stop else ("stop_after" if stop_ext else f"NO_STOP(next={tr(nxt)})") + (f";internal_stops={internal}" if internal else ""), "|".join(sig)]) + "\n")
    fa.write(f">{rec.replace(':', '_')}|{g}|{c}:{ex[0][0]}-{ex[-1][1]}({st})\n{p.rstrip('*')}\n")
