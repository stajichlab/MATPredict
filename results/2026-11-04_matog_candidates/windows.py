"""Cut a +/-20 kb window (clipped to the contig) around each called MATog locus. Run on the HPCC: python3 windows.py OUTDIR"""
import gzip, sys, os
L = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/"
# genome, slug, [(contig, start, end)]
LOCI = [
 ("GCA_030570775.1_ASM3057077v1", "pignaliae", [("JAJLWO010000194.1", 3, 7582), ("JAJLWO010000212.1", 2579, 4131)]),
 ("GCA_030559035.1_ASM3055903v1", "wickerhamii", [("JAKTQX010000069.1", 1, 2378), ("JAKTQX010000058.1", 47, 3114)]),
 ("GCA_030568015.1_ASM3056801v1", "dorogensis", [("JAJLTE010000037.1", 1, 2040), ("JAJLTE010000018.1", 14828, 17761)]),
 ("GCA_030574315.1_ASM3057431v1", "siamensis", [("JAJMET010000031.1", 3, 1978), ("JAJMET010000024.1", 20, 2966)]),
 ("GCA_002008015.1_ASM200801v1", "boidinii_alpha", [("MSSE01000022.1", 111232, 119719)]),
 ("GCA_002008025.1_ASM200802v1", "boidinii_a", [("MSRY01000009.1", 302361, 310790)]),
]
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
for gen, slug, loci in LOCI:
    seqs, name = {}, None
    for line in gzip.open(L + gen + ".fa.gz", "rt"):
        if line[0] == ">": name = line[1:].split()[0]; seqs[name] = []
        else: seqs[name].append(line.strip())
    for c, a, b in loci:
        s = "".join(seqs[c]); lo = max(1, a - 20000); hi = min(len(s), b + 20000)
        fn = f"{out}/{slug}__{c}_{lo}-{hi}.fa"
        open(fn, "w").write(f">{c}:{lo}-{hi}\n" + "\n".join(s[lo-1:hi][i:i+80] for i in range(0, hi-lo+1, 80)) + "\n")
        print(slug, gen, c, "len", len(s), "window", lo, hi)
