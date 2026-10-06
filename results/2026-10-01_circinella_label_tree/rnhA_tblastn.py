"""tblastn of the five Circinella-group queries against rnhA +-20 kb in every group genome."""
import csv, subprocess, sys
S = sys.argv[1]
sys.argv = [sys.argv[0]]
exec(open('extract.py').read().split('rows = list(')[0])
G = [r for r in csv.DictReader(open('genomes.tsv'), delimiter='\t') if r['group'] == 'circinella_group']
lt = {r['gid']: r for r in csv.DictReader(open('loci.tsv'), delimiter='\t')}
out = open('rnhA_tblastn.tsv', 'w')
out.write('gid\tlabel\tquery\tpident\tqstart\tqend\tqlen\tsstart\tsend\tevalue\tbits\tdist_rnhA_kb\n')
for g in G:
    rn = [m for m in mp_rows(g['gid']) if m['gene'] == 'rnhA']
    if not rn:
        out.write(f"{g['gid']}\t{g['label']}\tNO_RNHA\n"); continue
    m = max(rn, key=lambda m: float(m['qcov']) * float(m['identity']))
    c, s, e = m['contig'], int(m['start']), int(m['end'])
    cs = read_contig(g['genome'], c)
    lo, hi = max(0, s - 20000), min(len(cs), e + 20000)
    reg = f"{S}/{g['gid']}.fa"
    open(reg, 'w').write(f">{c}\n{cs[lo:hi]}\n")
    p = subprocess.run([f'{E}/tblastn', '-query', 'circ_queries2.faa', '-subject', reg, '-evalue', '1e-3', '-seg', 'no',
                        '-outfmt', '6 qseqid pident qstart qend qlen sstart send evalue bitscore'], capture_output=True, text=True)
    best = {}
    for line in p.stdout.splitlines():
        f = line.split('\t')
        if f[0] not in best or float(f[8]) > float(best[f[0]][8]):
            best[f[0]] = f
    if not best:
        out.write(f"{g['gid']}\t{g['label']}\tNONE\n")
    for q, f in best.items():
        a, b = lo + min(int(f[5]), int(f[6])), lo + max(int(f[5]), int(f[6]))
        d = max(0, max(a, s) - min(b, e)) / 1000
        out.write('\t'.join([g['gid'], g['label'], q.split('|')[0]] + f[1:5] + [str(a), str(b), f[7], f[8], f'{d:.1f}']) + '\n')
out.close()
