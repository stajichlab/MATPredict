#!/usr/bin/env python3
"""Cut the HD1-like/HD2 locus region of each AMF genome allele into one reference FASTA.

For each genome, every HD1+HD2 pair on one contig (stage-1 control classes HD1, HD2) gives one reference
named <label>_a<k>: the span of both genes plus --pad bp each side (flank genes are not included).
Usage: amf_allele_refs.py BASE OUT.fa GENOME:LABEL [GENOME:LABEL ...]   (needs samtools on PATH)
"""
import os, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from homeodomain_tree import table
ANN = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation"
base, out = sys.argv[1], sys.argv[2]
PAD = 500
with open(out, "w") as fo:
    for spec in sys.argv[3:]:
        g, _, label = spec.partition(":")
        D = f"{base}/stage1_all/per_genome"
        pos = {r[0]: (r[1], int(r[2]), int(r[3]), int(r[5])) for r in table(f"{D}/{g}.genes_all.tsv")}
        cand = table(f"{D}/{g}.candidates.tsv")
        h1 = [c[0] for c in cand if c[5] == "HD1"]; h2 = [c[0] for c in cand if c[5] == "HD2"]
        k = 0
        for b in h2:
            for a in h1:
                if pos[a][0] == pos[b][0] and abs(pos[a][3] - pos[b][3]) <= 1:
                    k += 1
                    s = max(1, min(pos[a][1], pos[b][1]) - PAD); e = max(pos[a][2], pos[b][2]) + PAD
                    fa = f"{ANN}/{g}/predict_results/{g}.scaffolds.fa.gz"
                    r = subprocess.run(["samtools", "faidx", fa, f"{pos[a][0]}:{s}-{e}"], capture_output=True, text=True, check=True)
                    seq = "".join(r.stdout.splitlines()[1:])
                    fo.write(f">{label}_a{k}|{pos[a][0]}:{s}-{e}\n{seq}\n")
        print(label, g, "alleles", k)
