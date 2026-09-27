"""Step 2: cluster sexM/sexP tblastn HSPs per genome and tag which fall in a reported locus.

Source HSPs: the flank-synteny run's genome-wide tblastn of the curated Mucoromycota
MAT proteins (hits/<genome>.tbn.zst; -evalue 10, -seg no, max_target_seqs 50).
HMG cluster = sexM/sexP HSPs on one contig within 3 kb of each other.
A cluster counts as an HMG copy if its best HSP has E <= 1e-3.
Writes clusters.tsv (all copies) with the locus_id(s) each overlaps.
"""
import collections, csv, glob, io, os, subprocess
S = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_flank_synteny_ED/hits"
GAP, EMAX = 3000, 1e-3
loci = collections.defaultdict(list)
for r in csv.DictReader(open("loci.tsv"), delimiter="\t"):
    loci[r["genome"]].append(r)
genomes = sorted({os.path.basename(f)[:-8] for f in glob.glob(f"{S}/*.tbn.zst")})
out = open("clusters.tsv", "w")
out.write("cluster_id\tgenome\tcontig\tstart\tend\tbest_query\tbest_gene\tbest_evalue\tbest_bits\tn_hsp\tlocus_ids\n")
n_in = n_out = 0
for g in genomes:
    txt = subprocess.run(["zstdcat", f"{S}/{g}.tbn.zst"], capture_output=True, text=True).stdout
    hs = collections.defaultdict(list)
    for line in txt.splitlines():
        f = line.split("\t")
        gene = f[0].split("|")[2]
        if gene not in ("sexM", "sexP"):
            continue
        s, e = sorted((int(f[6]), int(f[7])))
        hs[f[1]].append((s, e, f[0], gene, float(f[8]), float(f[9])))
    k = 0
    for contig, h in hs.items():
        h.sort()
        groups, cur, end = [], [], -1
        for x in h:
            if cur and x[0] > end + GAP:
                groups.append(cur); cur = []
            cur.append(x); end = max(end, x[1])
        groups.append(cur)
        for grp in groups:
            best = min(grp, key=lambda x: (x[4], -x[5]))
            if best[4] > EMAX:
                continue
            s, e = min(x[0] for x in grp), max(x[1] for x in grp)
            ids = [L["locus_id"] for L in loci.get(g, [])
                   if L["contig"] == contig and int(L["start"]) <= e and int(L["end"]) >= s]
            k += 1
            out.write(f"{g}|h{k}\t{g}\t{contig}\t{s}\t{e}\t{best[2]}\t{best[3]}\t{best[4]:.3g}\t{best[5]}\t{len(grp)}\t{','.join(ids)}\n")
            n_in += bool(ids); n_out += not ids
print("genomes", len(genomes), "HMG copies in loci", n_in, "outside loci", n_out)
