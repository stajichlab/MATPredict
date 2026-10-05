"""Order and orientation of APC5, SLA2, COX13, APN2 (and any HMG gene) around SLA2/APN2.

Queries are the NC1011 proteins (KAI0195599.1 APC5-like, KAI0195601.1 SLA2, KAI0195602.1 COX13,
KAI0195603.1 APN2). Each genome's proteome is searched with diamond; the best hit per query gives
coordinates from the NCBI GFF. Output: one row per genome with the gene order along the contig.
"""
import subprocess, sys, os, io, zipfile, urllib.request, collections, csv
D = sys.argv[1]
Q = {'APC5': 'KAI0195599.1', 'SLA2': 'KAI0195601.1', 'COX13': 'KAI0195602.1', 'APN2': 'KAI0195603.1'}
W = f'{D}/synteny'; os.makedirs(W, exist_ok=True)
rows = [l.rstrip('\n').split('\t') for l in open(f'{D}/data/tree_genomes.tsv')][1:]

def fetch(acc):
    z = f'{W}/{acc}.zip'
    if not os.path.exists(z):
        url = f'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download?include_annotation_type=GENOME_GFF&include_annotation_type=PROT_FASTA&filename=x.zip'
        urllib.request.urlretrieve(url, z)
    zf = zipfile.ZipFile(z)
    names = zf.namelist()
    g = [n for n in names if n.endswith('genomic.gff')]; p = [n for n in names if n.endswith('protein.faa')]
    if not g or not p: return None
    open(f'{W}/{acc}.gff', 'wb').write(zf.read(g[0])); open(f'{W}/{acc}.faa', 'wb').write(zf.read(p[0]))
    return True

# query proteins from NC1011
nc = 'GCA_022453505.1'
fetch(nc)
seqs = {}; k = None
for l in open(f'{W}/{nc}.faa'):
    if l[0] == '>': k = l[1:].split()[0]; seqs[k] = []
    else: seqs[k].append(l.strip())
with open(f'{W}/queries.faa', 'w') as o:
    for n, a in Q.items(): o.write(f'>{n}\n{"".join(seqs[a])}\n')

def cds_coords(gff):
    c = {}
    for l in open(gff):
        if l[0] == '#': continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'CDS': continue
        a = dict(x.split('=', 1) for x in f[8].split(';') if '=' in x)
        n = a.get('Name') or a.get('protein_id')
        if not n: continue
        s, e = int(f[3]), int(f[4])
        if n in c: c[n][1] = min(c[n][1], s); c[n][2] = max(c[n][2], e)
        else: c[n] = [f[0], s, e, f[6]]
    return c

out = open(f'{D}/synteny_order.tsv', 'w'); out.write('accession\tlabel\tgroup\torder_along_contig\tcontig\tspan_kb\tn_genes_found\n')
detail = open(f'{D}/synteny_detail.tsv', 'w'); detail.write('accession\tlabel\tgene\tprotein\tcontig\tstart\tend\tstrand\tpident\tevalue\n')
allrows = [(nc, 'Xylaria_flabelliformis_NC1011', 'XYL_ours')] + [(a, l, g) for a, l, g in rows if a != nc]
for acc, label, grp in allrows:
    try:
        if not fetch(acc): out.write(f'{acc}\t{label}\t{grp}\tno_annotation\t\t\t0\n'); continue
    except Exception as e:
        out.write(f'{acc}\t{label}\t{grp}\tdownload_failed\t\t\t0\n'); continue
    subprocess.run(['diamond', 'makedb', '--in', f'{W}/{acc}.faa', '-d', f'{W}/{acc}', '--quiet'], check=True)
    r = subprocess.run(['diamond', 'blastp', '-q', f'{W}/queries.faa', '-d', f'{W}/{acc}', '--evalue', '1e-5', '--max-target-seqs', '5',
                        '--outfmt', '6', 'qseqid', 'sseqid', 'pident', 'evalue', 'bitscore', '--quiet', '--sensitive'], capture_output=True, text=True)
    best = {}
    for l in r.stdout.splitlines():
        q, s, pid, ev, bs = l.split('\t'); bs = float(bs)
        if q not in best or bs > best[q][3]: best[q] = (s, float(pid), float(ev), bs)
    c = cds_coords(f'{W}/{acc}.gff'); found = []
    for q, (s, pid, ev, bs) in best.items():
        if s in c:
            ct, st, en, sd = c[s]; found.append((q, s, ct, st, en, sd, pid, ev))
            detail.write(f'{acc}\t{label}\t{q}\t{s}\t{ct}\t{st}\t{en}\t{sd}\t{pid}\t{ev}\n')
    # contig carrying most genes; order along it
    cnt = collections.Counter(f[2] for f in found)
    if not cnt: out.write(f'{acc}\t{label}\t{grp}\tnone_found\t\t\t0\n'); continue
    ct = cnt.most_common(1)[0][0]; on = sorted([f for f in found if f[2] == ct], key=lambda f: f[3])
    order = ' '.join(f'{f[0]}({f[5]})' for f in on)
    span = (max(f[4] for f in on) - min(f[3] for f in on)) / 1000
    out.write(f'{acc}\t{label}\t{grp}\t{order}\t{ct}\t{span:.1f}\t{len(found)}\n')
    out.flush()
print('done')
