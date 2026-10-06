"""Model the rnhA-adjacent sex gene in every Circinella-group genome.

miniprot (--trans) of the five Circinella-group queries (circ_queries2.faa)
on rnhA +-20 kb; keep the best model by AS; score with the classifier.
Also report the annotated mRNA overlapping the model. Writes rnhA_models.tsv/.faa.
"""
import csv, subprocess, sys
S = sys.argv[1]
sys.argv = [sys.argv[0]]
exec(open('extract.py').read().split('rows = list(')[0])
G = [r for r in csv.DictReader(open('genomes.tsv'), delimiter='\t') if r['group'] == 'circinella_group']
out = open('rnhA_models.tsv', 'w'); fa = open('rnhA_models.faa', 'w')
out.write('gid\tlabel\tquery\tqcov\tAS\tcontig\tstart\tend\tdist_rnhA_kb\tlen\tsexM\tsexP\tP1\tannot_id\tannot_len\tannot_sexM\tannot_sexP\tannot_P1\n')
for g in G:
    rn = [m for m in mp_rows(g['gid']) if m['gene'] == 'rnhA']
    if not rn:
        continue
    m = max(rn, key=lambda m: float(m['qcov']) * float(m['identity']))
    c, s, e = m['contig'], int(m['start']), int(m['end'])
    cs = read_contig(g['genome'], c)
    lo, hi = max(0, s - 20000), min(len(cs), e + 20000)
    reg = f"{S}/{g['gid']}_m.fa"
    open(reg, 'w').write(f">{c}\n{cs[lo:hi]}\n")
    p = subprocess.run([f'{E}/miniprot', '--trans', '-p', '0.2', '-N', '10', reg, 'circ_queries2.faa'], capture_output=True, text=True)
    best, cur = None, None
    for line in p.stdout.splitlines():
        if not line.startswith('#') and '\t' in line:
            f = line.split('\t')
            ms = [x for x in f if x.startswith('AS:i:')]
            cur = {'query': f[0].split('|')[0], 'qlen': int(f[1]), 'qs': int(f[2]), 'qe': int(f[3]),
                   'ts': lo + int(f[7]), 'te': lo + int(f[8]), 'AS': int(ms[0][5:]) if ms else 0}
        elif line.startswith('##STA') and cur:
            cur['protein'] = line.split('\t')[1].strip().replace('*', '')
            if best is None or cur['AS'] > best['AS']:
                best = cur
            cur = None
    if not best:
        out.write(f"{g['gid']}\t{g['label']}\tNONE\n"); continue
    sc = score(best['protein'])
    d = max(0, max(best['ts'], s) - min(best['te'], e)) / 1000
    aid, aseq = annotated(g.get('gff'), g.get('proteins'), c, best['ts'], best['te'])
    asc = score(aseq) if aseq else {}
    out.write('\t'.join(map(str, [g['gid'], g['label'], best['query'], round((best['qe'] - best['qs']) / best['qlen'], 2), best['AS'],
              c, best['ts'], best['te'], f'{d:.1f}', len(best['protein']), sc['sexM'], sc['sexP'], sc['P1'],
              aid or '', len(aseq) if aseq else '', asc.get('sexM', ''), asc.get('sexP', ''), asc.get('P1', '')])) + '\n')
    fa.write(f">{g['gid']}\n{best['protein']}\n")
    if aseq:
        fa.write(f">{g['gid']}|annot|{aid}\n{aseq.replace('*','')}\n")
out.close(); fa.close()
