"""Collapse tblastn HSPs of flank genes into regions and write one protein per region.

Input: tblastn table (qseqid sseqid pident length qstart qend sstart send evalue bitscore sseq).
Region = HSPs of one gene on one contig within GAP bp of each other.
Region protein = HSPs of the best-scoring query in that region, ordered by qstart,
gaps and stops removed. This is used only for the reverse (back) search.
"""
import collections, sys
FLANK = {"tptA", "rnhA", "glrA", "algA", "btbA"}
GAP = 3000
tab, genome, out = sys.argv[1], sys.argv[2], sys.argv[3]
hsps = collections.defaultdict(list)
for line in open(tab):
    f = line.rstrip("\n").split("\t")
    gene = f[0].split("|")[2]
    if gene not in FLANK:
        continue
    s, e = sorted((int(f[6]), int(f[7])))
    hsps[(gene, f[1])].append((s, e, f[0], int(f[4]), float(f[9]), f[10]))
with open(out, "w") as fo:
    for (gene, contig), hs in hsps.items():
        hs.sort()
        groups, cur, end = [], [], -1
        for h in hs:
            if cur and h[0] > end + GAP:
                groups.append(cur); cur = []
            cur.append(h); end = max(end, h[1]) if cur[:-1] else h[1]
        groups.append(cur)
        for g in groups:
            byq = collections.defaultdict(list)
            for h in g:
                byq[h[2]].append(h)
            q, qh = max(byq.items(), key=lambda kv: sum(x[4] for x in kv[1]))
            bits = sum(x[4] for x in qh)
            seq = "".join(x[5] for x in sorted(qh, key=lambda x: x[3])).replace("-", "").replace("*", "")
            if len(seq) < 30:
                continue
            s, e = min(x[0] for x in g), max(x[1] for x in g)
            fo.write(f">{genome}|{gene}|{contig}|{s}|{e}|{bits:.1f}\n{seq}\n")
