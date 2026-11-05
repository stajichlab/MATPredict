"""Cut a window (clipped to the contig) around each locus of the two candidate genomes. Run on the HPCC: python3 windows.py OUTDIR"""
import gzip, sys, os
L = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
LOCI = [
 ("GCA_003707675.3_ASM370767v3", "philentoma", [("PPKZ03000004.1", 215815, 223380), ("PPKZ03000012.1", 331460, 336694)]),
 ("GCA_054949465.1_ASM5494946v1", "edaphicus", [("JBSXOT010000020.1", 28026, 44107), ("JBSXOT010000274.1", 588, 2372)]),
]
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
for gen, slug, loci in LOCI:
    seqs, name = {}, None
    for line in gzip.open(L + gen + ".fa.gz", "rt"):
        if line[0] == ">": name = line[1:].split()[0]; seqs[name] = []
        else: seqs[name].append(line.strip())
    for c, a, b in loci:
        s = "".join(seqs[c]); lo = max(1, a - 15000); hi = min(len(s), b + 15000)
        open(f"{out}/{slug}__{c}_{lo}-{hi}.fa", "w").write(f">{c}:{lo}-{hi}\n" + "\n".join(s[lo-1:hi][i:i+80] for i in range(0, hi-lo+1, 80)) + "\n")
        print(slug, c, "len", len(s), "window", lo, hi)
