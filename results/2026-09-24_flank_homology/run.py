"""SLA2/APN2 homology case for Pezizales, on annotated proteomes.

For each annotated Pezizales proteome: blastp the Morchella importuna (Pezizales,
KY782629/KY782630) and A. fumigatus (curated record) SLA2, APN2, MAT1-1-1 and
MAT1-2-1 proteins; keep the best hit per gene, its E-value, and the E-value of
the next-best DIFFERENT protein (a paralog gap). The best hit is checked
reciprocally against the A. fumigatus Af293 proteome. Positions come from the
assembly's GFF, so SLA2/APN2 can be placed relative to the MAT gene.
"""
import collections, gzip, os, re, subprocess, sys, glob
D = os.path.dirname(os.path.abspath(__file__))
BLAST = "/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin"
REC = "/bigdata/stajichlab/jstajich/projects/MATPredict/.claude/worktrees/polish-scope-cuts/db/Ascomycota"

def rd(path):
    s, k = {}, None
    op = gzip.open if path.endswith(".gz") else open
    for l in op(path, "rt"):
        if l.startswith(">"):
            k = l[1:].split()[0]; s[k] = []
        else:
            s[k].append(l.strip())
    return {k: "".join(v) for k, v in s.items()}

def sh(cmd):
    return subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True).stdout

# --- queries
q = {}
m29, m30 = rd(f"{D}/cds_KY782629.1.faa"), rd(f"{D}/cds_KY782630.1.faa")
pick = lambda s, pid: next(v for k, v in s.items() if pid in k)
q["Mimp|SLA2"] = pick(m29, "AVI60811")
q["Mimp|APN2"] = pick(m29, "AVI60802") + pick(m29, "AVI60803")   # APN2 split into two CDS
q["Mimp|MAT1-2-1"] = pick(m29, "AVI60809")
q["Mimp|MAT1-1-1"] = pick(m30, "AVI60816")
afu = {}
for rec in ["Eurotiales/746128_a1163_MAT_MAT1-1", "Eurotiales/746128_af293_MAT_MAT1-2"]:
    for k, v in rd(f"{REC}/{rec}/proteins.faa").items():
        name = re.search(r"name=([^|]+)", k).group(1)
        if name in ("SLA2", "APN2", "MAT1-1-1", "MAT1-2-1"):
            afu.setdefault(name, v)
for n, v in afu.items():
    q[f"Afum|{n}"] = v
with open(f"{D}/queries.faa", "w") as fh:
    for k, v in q.items():
        fh.write(f">{k}\n{v}\n")

# --- Af293 proteome for the reciprocal check; which Af293 protein IS SLA2/APN2
if not os.path.exists(f"{D}/afu.pdb"):
    sh(f"zcat {D}/afu.faa.gz > {D}/afu.faa && {BLAST}/makeblastdb -in {D}/afu.faa -dbtype prot -out {D}/afu >/dev/null")
afu_id = {}
for line in sh(f"{BLAST}/blastp -query {D}/queries.faa -db {D}/afu -evalue 1e-20 -outfmt '6 qseqid sseqid bitscore' -max_target_seqs 1").splitlines():
    qq, s, b = line.split("\t")
    if qq.startswith("Afum|"):
        afu_id.setdefault(qq.split("|")[1], s)
print("Af293 orthologs:", afu_id, file=sys.stderr)

def gff_pos(path):
    pos = {}
    for l in gzip.open(path, "rt"):
        if l.startswith("#"): continue
        f = l.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        m = re.search(r"protein_id=([^;]+)", f[8])
        if not m: continue
        p = m.group(1); a, b = int(f[3]), int(f[4])
        c, lo, hi = pos.get(p, (f[0], a, b))
        pos[p] = (c, min(lo, a), max(hi, b))
    return pos

rows = []
for faa in sorted(glob.glob(f"{D}/proteomes/*.faa.gz")):
    g = os.path.basename(faa)[:-7]
    db = f"{D}/work/{g}"
    os.makedirs(f"{D}/work", exist_ok=True)
    if not os.path.exists(db + ".pdb"):
        sh(f"zcat {faa} > {db}.faa && {BLAST}/makeblastdb -in {db}.faa -dbtype prot -out {db} >/dev/null")
    hits = collections.defaultdict(list)
    for line in sh(f"{BLAST}/blastp -query {D}/queries.faa -db {db} -evalue 1e-3 -outfmt '6 qseqid sseqid pident length qlen evalue bitscore' -max_target_seqs 10").splitlines():
        qq, s, pid, ln, ql, ev, bits = line.split("\t")
        hits[qq].append((s, float(pid), int(ln), int(ql), float(ev), float(bits)))
    pos = gff_pos(f"{D}/proteomes/{g}.gff.gz")
    seqs = rd(f"{db}.faa")
    for qq, hs in hits.items():
        # best per subject protein
        best = {}
        for h in hs:
            if h[0] not in best or h[5] > best[h[0]][5]:
                best[h[0]] = h
        ranked = sorted(best.values(), key=lambda h: -h[5])
        top = ranked[0]; nxt = ranked[1] if len(ranked) > 1 else None
        # reciprocal: top hit vs Af293
        with open(f"{D}/work/rbh.faa", "w") as fh:
            fh.write(f">{top[0]}\n{seqs[top[0]]}\n")
        r = sh(f"{BLAST}/blastp -query {D}/work/rbh.faa -db {D}/afu -evalue 1e-5 -outfmt '6 sseqid evalue' -max_target_seqs 1").split("\n")[0].split("\t")
        gene = qq.split("|")[1]
        rbh = (r[0] == afu_id.get(gene)) if r and r[0] else False
        rows.append(dict(genome=g, query=qq, gene=gene, hit=top[0], pident=top[1], cov=round(100*top[2]/top[3]),
                         evalue=top[4], next_evalue=(nxt[4] if nxt else None), rbh_afu=rbh, pos=pos.get(top[0])))

with open(f"{D}/besthits.tsv", "w") as fh:
    fh.write("genome\tquery\thit\tpident\tqcov\tevalue\tnext_evalue\trbh_afu\tcontig\tstart\tend\n")
    for r in rows:
        c, a, b = r["pos"] or ("", "", "")
        fh.write(f"{r['genome']}\t{r['query']}\t{r['hit']}\t{r['pident']}\t{r['cov']}\t{r['evalue']:.1e}\t"
                 f"{'' if r['next_evalue'] is None else f'{r['next_evalue']:.1e}'}\t{r['rbh_afu']}\t{c}\t{a}\t{b}\n")
print(f"wrote besthits.tsv ({len(rows)} rows)", file=sys.stderr)
