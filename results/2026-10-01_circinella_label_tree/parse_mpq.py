"""Parse genome-wide miniprot hits of the three Circinella-derived queries.

For each genome and query type: the best hit by AS, with its translated
protein, location, distance to the best rnhA hit, and classifier scores.
Hits of the sexM-type and sexMstrong queries at one locus are merged (best AS).
Writes mpq_hits.tsv and mpq_best.faa.
"""
import csv, sys
sys.argv = [sys.argv[0]]
exec(open('extract.py').read().split('rows = list(')[0])
lt = {r['gid']: r for r in csv.DictReader(open('loci.tsv'), delimiter='\t')}
G = {r['gid']: r for r in csv.DictReader(open('genomes.tsv'), delimiter='\t')}
out = open('mpq_hits.tsv', 'w')
fa = open('mpq_best.faa', 'w')
out.write('gid\tlabel\tqtype\tquery\tcontig\tstart\tend\tAS\tqcov\tident\tlen\tdist_rnhA_kb\tsexM\tsexP\tP1\n')
for gid, g in G.items():
    p = Path(f'mpq/{gid}.paf')
    if not p.exists():
        continue
    hits, cur = [], None
    for line in open(p):
        if not line.startswith('#') and '\t' in line:
            f = line.rstrip('\n').split('\t')
            tags = {x[:2]: x[5:] for x in f[12:] if x[2:5] in (':i:', ':f:')}
            cur = {'query': f[0], 'qlen': int(f[1]), 'qs': int(f[2]), 'qe': int(f[3]), 'contig': f[5],
                   'start': int(f[7]), 'end': int(f[8]), 'AS': int(tags.get('AS', 0)),
                   'ident': round(int(f[9]) / max(1, int(f[10])), 3)}
        elif line.startswith('##STA') and cur:
            cur['protein'] = line.split('\t')[1].strip().replace('*', '')
            hits.append(cur); cur = None
    l = lt.get(gid, {})
    for qt in ('sexPtype', 'sexMtype'):
        hs = [h for h in hits if (qt == 'sexPtype') == h['query'].startswith('CIRC_sexPtype')]
        if not hs:
            continue
        h = max(hs, key=lambda x: x['AS'])
        d = ''
        if l.get('rnhA_contig') == h['contig'] and l.get('rnhA_start'):
            d = round(max(0, max(h['start'], int(l['rnhA_start'])) - min(h['end'], int(l['rnhA_end']))) / 1000, 1)
        sc = score(h['protein'])
        qcov = round((h['qe'] - h['qs']) / h['qlen'], 2)
        out.write('\t'.join(map(str, [gid, g['label'], qt, h['query'].split('|')[0], h['contig'], h['start'], h['end'], h['AS'],
                                      qcov, h['ident'], len(h['protein']), d, sc['sexM'], sc['sexP'], sc['P1']])) + '\n')
        fa.write(f'>{gid}|{qt}\n{h["protein"]}\n')
out.close(); fa.close()
