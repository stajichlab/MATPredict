"""Best-hit localisation restricted to each call's OWN core gene names.
For each gene in the call: rank of the best in-locus hit among all genome-wide hits
(distinct loci, merged within 5 kb). rank 1 = locus holds the genome's best hit.
Writes classified2.tsv (one row per call x gene)."""
import csv, collections, subprocess
PAD = 20000
rows = list(csv.DictReader(open("changed.tsv"), delimiter="\t"))
w = csv.writer(open("classified2.tsv", "w"), delimiter="\t")
w.writerow(["panel", "genome", "order", "species", "mat_family", "outcome", "call", "gene", "domain",
            "rank_here", "n_loci_genomewide", "best_here_bits", "best_here_pid", "best_bits", "best_where"])
HMG = ("MAT1-2", "sexM", "sexP", "MTLa2", "MTLA2", "mata2")
cache = {}
for r in rows:
    g = r["genome"]
    if g not in cache:
        cache[g] = [l.split("\t") for l in subprocess.run(["zstdcat", f"hits/{g}.tsv.zst"], capture_output=True, text=True).stdout.splitlines()]
    for gene in r["core_genes"].split(","):
        hits = [(float(h[9]), h[1], *sorted((int(h[6]), int(h[7]))), float(h[2])) for h in cache[g]
                if h[0].split("|")[-1] == gene and float(h[8]) <= 1e-3]
        # merge into loci: best bitscore per (contig, 5kb bin cluster)
        loci = []
        for bs, c, lo, hi, pid in sorted(hits, reverse=True):
            if any(c == lc and abs(lo - llo) < 5000 for _, lc, llo, _, _ in loci): continue
            loci.append((bs, c, lo, hi, pid))
        c0, s0, e0 = r["contig"], int(r["start"]) - PAD, int(r["end"]) + PAD
        rank = next((i + 1 for i, (bs, c, lo, hi, pid) in enumerate(loci) if c == c0 and s0 <= lo and hi <= e0), None)
        here = loci[rank - 1] if rank else None
        dom = "HMG" if gene.startswith(HMG) else "other"
        w.writerow([r["panel"], g, r["order"], r["species"], r["mat_family"], r["outcome"],
                    f'{c0}:{r["start"]}-{r["end"]}', gene, dom, rank or "absent", len(loci),
                    here[0] if here else "", here[4] if here else "", loci[0][0] if loci else "",
                    f"{loci[0][1]}:{loci[0][2]}" if loci else ""])
