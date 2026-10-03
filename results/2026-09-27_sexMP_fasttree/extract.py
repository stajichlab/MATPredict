"""Step 3: re-extract candidate proteins for one genome.

Windows: every sexM/sexP locus span +-5 kb (candidates), and each chosen non-locus
HMG copy +-5 kb (outgroups). miniprot (--trans) aligns the 15 curated sexM/sexP
proteins to each window; overlapping alignments are merged and the best-scoring
one (AS) kept per genomic interval. A window with no miniprot model falls back to
tblastn: HSPs clustered by genomic position (3 kb gap); per cluster the best-scoring
query's non-overlapping HSP chain, gaps and stops removed (method=tblastn; >= 30 aa).
Genetic code 1 for every genome here (all reports say genetic_code 1).
Usage: extract.py GENOME  -> extract2/GENOME.faa, extract2/GENOME.tsv
"""
import csv, os, subprocess, sys, gzip, collections
g = sys.argv[1]
D = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-27_sexMP_fasttree"
P = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_sexMP_phylogeny"
LIB = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes"
BIN = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
W = os.path.join(os.environ["SCRATCH"], "sexmp_ft", g); os.makedirs(W, exist_ok=True)
PAD = 5000
wins = []
for r in csv.DictReader(open(f"{D}/loci.tsv"), delimiter="\t"):
    if r["genome"] == g:
        wins.append(("cand", r["locus_id"], r["contig"], int(r["start"]), int(r["end"])))
# outgroup copies from the earlier run; drop any that a (new) locus span now covers
spans = [(w[2], w[3], w[4]) for w in wins]
for r in csv.DictReader(open(f"{P}/outgroup_clusters.tsv"), delimiter="\t"):
    if r["genome"] == g:
        c, s0, e0 = r["contig"], int(r["start"]), int(r["end"])
        if any(c == c2 and s0 <= e2 and e0 >= s2 for c2, s2, e2 in spans):
            continue
        wins.append(("nonlocus", r["cluster_id"], c, s0, e0))
need = {w[2] for w in wins}
seq, name = {}, None
with gzip.open(f"{LIB}/{g}.fa.gz", "rt") as fh:
    for l in fh:
        if l[0] == ">":
            name = l[1:].split()[0]; seq[name] = [] if name in need else None
        elif seq.get(name) is not None:
            seq[name].append(l.strip())
seq = {k: "".join(v) for k, v in seq.items() if v is not None}
wf = f"{W}/win.fna"
with open(wf, "w") as fo:
    for i, (kind, wid, c, s, e) in enumerate(wins):
        a = max(0, s - 1 - PAD); b = min(len(seq[c]), e + PAD)
        fo.write(f">w{i}\n{seq[c][a:b]}\n")
mp = subprocess.run([f"{BIN}/miniprot", "--trans", "-t1", wf, f"{D}/refs_sexMP.faa"],
                    capture_output=True, text=True).stdout.splitlines()
alns = collections.defaultdict(list)
for i, l in enumerate(mp):
    if l.startswith("##STA"):
        p = mp[i - 1].split("\t")
        AS = int([x for x in p if x.startswith("AS:i:")][0][5:])
        qlen, qs, qe = int(p[1]), int(p[2]), int(p[3])
        alns[p[5]].append(dict(ts=int(p[7]), te=int(p[8]), AS=AS, q=p[0], cov=(qe - qs) / qlen,
                               prot=l.split("\t")[1].rstrip("*").replace("*", "X")))
rows, faa = [], []
nofb = []
for i, (kind, wid, c, s, e) in enumerate(wins):
    a = alns.get(f"w{i}", [])
    a.sort(key=lambda x: -x["AS"])
    kept = []
    for x in a:
        if all(x["te"] <= y["ts"] or x["ts"] >= y["te"] for y in kept):
            kept.append(x)
    if not kept:
        nofb.append(i); continue
    for j, x in enumerate(kept):
        cid = f"{wid}|m{j}"
        rows.append([cid, g, kind, wid, "miniprot", x["q"], x["AS"], round(x["cov"], 3), len(x["prot"])])
        faa.append((cid, x["prot"]))
if nofb:
    sub = f"{W}/fb.fna"
    lines = open(wf).read().split("\n")
    keep = {f">w{i}" for i in nofb}
    with open(sub, "w") as fo:
        for k in range(0, len(lines) - 1, 2):
            if lines[k] in keep: fo.write(lines[k] + "\n" + lines[k + 1] + "\n")
    subprocess.run([f"{BIN}/makeblastdb", "-in", sub, "-dbtype", "nucl", "-out", f"{W}/fb"], capture_output=True)
    tb = subprocess.run([f"{BIN}/tblastn", "-query", f"{D}/refs_sexMP.faa", "-db", f"{W}/fb", "-evalue", "10",
                         "-seg", "no", "-outfmt", "6 qseqid sseqid qstart qend sstart send evalue bitscore sseq"],
                        capture_output=True, text=True).stdout.splitlines()
    by = collections.defaultdict(lambda: collections.defaultdict(list))
    for l in tb:
        f = l.split("\t"); by[f[1]][f[0]].append(f)
    for i in nofb:
        kind, wid = wins[i][0], wins[i][1]
        q = by.get(f"w{i}")
        if not q:
            rows.append([f"{wid}|none", g, kind, wid, "none", "", 0, 0, 0]); continue
        # cluster all HSPs of the window by genomic position (3 kb gap); one protein per cluster
        allh = sorted((min(int(x[4]), int(x[5])), max(int(x[4]), int(x[5])), x) for hs in q.values() for x in hs)
        clusters, cur, end = [], [], -1
        for h in allh:
            if cur and h[0] > end + 3000:
                clusters.append(cur); cur = []
            cur.append(h); end = max(end, h[1])
        clusters.append(cur)
        for j, cl in enumerate(clusters):
            byq = collections.defaultdict(list)
            for _, _, x in cl: byq[x[0]].append(x)
            best, hs = max(byq.items(), key=lambda kv: sum(float(x[7]) for x in kv[1]))
            hs.sort(key=lambda x: -float(x[7])); use = []
            for h in hs:
                if all(int(h[3]) < int(u[2]) or int(h[2]) > int(u[3]) for u in use): use.append(h)
            use.sort(key=lambda x: int(x[2]))
            prot = "".join(h[8] for h in use).replace("-", "").replace("*", "X")
            emin = min(float(h[6]) for h in use)
            if len(prot) < 30: continue
            cid = f"{wid}|t{j}"
            rows.append([cid, g, kind, wid, "tblastn", best, round(sum(float(x[7]) for x in use), 1), f"E={emin:.2g}", len(prot)])
            faa.append((cid, prot))
with open(f"{D}/extract2/{g}.faa", "w") as fo:
    for cid, p in faa: fo.write(f">{cid}\n{p}\n")
with open(f"{D}/extract2/{g}.tsv", "w") as fo:
    for r in rows: fo.write("\t".join(map(str, r)) + "\n")
subprocess.run(["rm", "-rf", W])
