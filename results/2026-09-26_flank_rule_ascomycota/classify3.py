"""With detect's tblastn settings (-seg no, evalue 10): for each call x own core gene,
the best in-locus hit's e-value and its rank among genome-wide loci of that gene.
usage: classify3.py HITDIR OUT.tsv"""
import csv, collections, subprocess, sys
PAD = 20000
hd, outp = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open("changed.tsv"), delimiter="\t"))
w = csv.writer(open(outp, "w"), delimiter="\t")
w.writerow(["panel", "genome", "order", "species", "mat_family", "outcome", "why", "max_core_dist_bp", "call", "gene",
            "rank_here", "n_loci_e1e-3", "here_evalue", "here_bits", "here_pid", "best_evalue", "best_where"])
cache = {}
for r in rows:
    g = r["genome"]
    if g not in cache:
        cache[g] = [l.split("\t") for l in subprocess.run(["zstdcat", f"{hd}/{g}.tsv.zst"], capture_output=True, text=True).stdout.splitlines()]
    for gene in r["core_genes"].split(","):
        hits = sorted(((float(h[9]), float(h[8]), h[1], *sorted((int(h[6]), int(h[7]))), float(h[2])) for h in cache[g]
                       if h[0].split("|")[-1] == gene), reverse=True)
        loci = []
        for bs, ev, c, lo, hi, pid in hits:
            if any(c == l[2] and abs(lo - l[3]) < 5000 for l in loci): continue
            loci.append((bs, ev, c, lo, hi, pid))
        c0, s0, e0 = r["contig"], int(r["start"]) - PAD, int(r["end"]) + PAD
        rank = next((i + 1 for i, l in enumerate(loci) if l[2] == c0 and s0 <= l[3] and l[4] <= e0), None)
        h = loci[rank - 1] if rank else None
        w.writerow([r["panel"], g, r["order"], r["species"], r["mat_family"], r["outcome"], r["why"], r["max_core_dist_bp"],
                    f'{c0}:{r["start"]}-{r["end"]}', gene, rank or "absent", sum(1 for l in loci if l[1] <= 1e-3),
                    h[1] if h else "", h[0] if h else "", h[5] if h else "", loci[0][1] if loci else "",
                    f"{loci[0][2]}:{loci[0][3]}" if loci else ""])
