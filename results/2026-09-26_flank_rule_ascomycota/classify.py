"""For each changed call: is the genome-wide best tblastn hit of each core gene at this locus?
A locus 'window' = call interval padded by 20 kb. Writes classified.tsv."""
import csv, collections, io, subprocess
PAD = 20000
rows = list(csv.DictReader(open("changed.tsv"), delimiter="\t"))
out = csv.writer(open("classified.tsv", "w"), delimiter="\t")
hdr = list(rows[0].keys()) + ["genes_best_here", "genes_best_elsewhere", "best_here_detail", "verdict"]
out.writerow(hdr)
cache = {}
for r in rows:
    g = r["genome"]
    if g not in cache:
        txt = subprocess.run(["zstdcat", f"hits/{g}.tsv.zst"], capture_output=True, text=True).stdout
        best = {}
        for l in txt.splitlines():
            q, s, pid, ln, qs, qe, ss, se, ev, bs, ql = l.split("\t")
            gene = q.split("|")[-1]; bs = float(bs)
            lo, hi = sorted((int(ss), int(se)))
            if gene not in best or bs > best[gene][0]:
                best[gene] = (bs, s, lo, hi, float(pid), float(ev))
        cache[g] = best
    best = cache[g]
    c, s, e = r["contig"], int(r["start"]) - PAD, int(r["end"]) + PAD
    here, away, det = [], [], []
    for gene, (bs, sc, lo, hi, pid, ev) in sorted(best.items()):
        if ev > 1e-5: continue
        if sc == c and lo >= s and hi <= e:
            here.append(gene); det.append(f"{gene}:{bs:.0f}:{pid:.0f}%:{ev:.0e}")
        else:
            away.append(gene)
    verdict = "best_hit_here" if here else ("best_hit_elsewhere" if away else "no_core_hit_1e-5")
    out.writerow(list(r.values()) + [",".join(here), ",".join(away), ";".join(det), verdict])
