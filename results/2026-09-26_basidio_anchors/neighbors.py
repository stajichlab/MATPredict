"""Candidate HD- and PR-locus anchors outside Agaricomycotina.
For each annotated genome: blastp the curated HD-class (and receptor-class)
proteins against the proteome; the best HD hit marks the HD locus. Collect
every protein whose gene lies within +-WIN bp of it. Then blastp neighbour
sets across genomes and count, per neighbour, in how many OTHER genomes a
homolog is also an HD-locus neighbour."""
import os, sys, subprocess, collections, re, glob
D = os.path.dirname(os.path.abspath(__file__)); A = f'{D}/annot'
WIN = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
CLASS = sys.argv[2] if len(sys.argv) > 2 else 'HD'
BIN = '/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin'
def run(c): subprocess.run(c, shell=True, check=True, env={**os.environ, 'PATH': BIN + ':' + os.environ['PATH']})
q = f'{D}/q_{CLASS}.faa'
with open(q, 'w') as o:
    keep = False
    for l in open(f'{D}/queries.faa'):
        if l.startswith('>'): keep = l[1:].startswith(CLASS + '|')
        if keep: o.write(l)
genomes = sorted(os.path.basename(p) for p in glob.glob(f'{A}/GC*') if os.path.isdir(p))
nb = {}; names = {}
for g in genomes:
    dd = glob.glob(f'{A}/{g}/ncbi_dataset/data/*/')[0]
    gff, faa = dd + 'genomic.gff', dd + 'protein.faa'
    loc = {}
    for l in open(gff):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS': continue
        m = re.search(r'protein_id=([^;]+)', f[8])
        if not m: continue
        pid = m.group(1); s, e = int(f[3]), int(f[4])
        prod = re.search(r'product=([^;]+)', f[8])
        if pid in loc: loc[pid] = (f[0], min(loc[pid][1], s), max(loc[pid][2], e), loc[pid][3])
        else: loc[pid] = (f[0], s, e, prod.group(1) if prod else '')
    db = f'{A}/{g}/prot'
    if not os.path.exists(db + '.pdb'): run(f'makeblastdb -in {faa} -dbtype prot -out {db} >/dev/null')
    out = subprocess.run(f'{BIN}/blastp -query {q} -db {db} -evalue 1e-5 -outfmt "6 qseqid sseqid pident bitscore" -max_target_seqs 5',
                         shell=True, capture_output=True, text=True).stdout.split('\n')
    hits = [l.split('\t') for l in out if l]
    if not hits: print(g, 'no', CLASS, 'hit'); continue
    best = max(hits, key=lambda h: float(h[3]))
    hdset = {h[1] for h in hits}
    c, s, e, p = loc[best[1]]
    near = [pid for pid, (cc, ss, ee, pp) in loc.items() if cc == c and ee >= s - WIN and ss <= e + WIN and pid not in hdset]
    print(f'{g}\tbest {CLASS} {best[1]} ({best[0].split("|")[1]}, {float(best[2]):.0f}%) at {c}:{s}-{e}; {len(near)} neighbours')
    nb[g] = near
    for pid in near: names[(g, pid)] = loc[pid][3]
    from Bio import SeqIO
    seqs = {r.id: r for r in SeqIO.parse(faa, 'fasta')}
    SeqIO.write([seqs[p] for p in near if p in seqs], f'{A}/{g}/nb_{CLASS}.faa', 'fasta')
# pool and all-vs-all
pool = f'{A}/pool_{CLASS}.faa'
with open(pool, 'w') as o:
    for g in nb:
        for l in open(f'{A}/{g}/nb_{CLASS}.faa'):
            o.write(l.replace('>', f'>{g}__', 1) if l.startswith('>') else l)
run(f'makeblastdb -in {pool} -dbtype prot -out {A}/pool_{CLASS} >/dev/null')
out = subprocess.run(f'{BIN}/blastp -query {pool} -db {A}/pool_{CLASS} -evalue 1e-10 -outfmt "6 qseqid sseqid pident bitscore"',
                     shell=True, capture_output=True, text=True).stdout.split('\n')
share = collections.defaultdict(set)
for l in out:
    if not l: continue
    a, b, pid, bits = l.split('\t')
    ga, gb = a.split('__')[0], b.split('__')[0]
    if ga != gb: share[a].add(gb)
rows = sorted(((len(v), k) for k, v in share.items()), reverse=True)
print(f'\nneighbours (+-{WIN} bp of best {CLASS} hit) shared with other genomes\' {CLASS} neighbourhoods:')
seen = set()
for n, k in rows:
    g, pid = k.split('__', 1)
    print(f'{n}\t{g}\t{pid}\t{names.get((g, pid), "")[:70]}\t{sorted(share[k])}')
