"""Per clade genome: distance from the MAT ORF (six-frame scan) to the best SLA2 / APN2 / COX13 hits (and to GAP1 / BGL1 / ORF1 when flank2/ exists).
usage: flank_join.py FLANK_DIR FLANK2_DIR -> clade_flank.tsv"""
import csv, os, re, sys
from collections import Counter
fd, fd2 = sys.argv[1:3]
gates = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_gates.tsv'), delimiter='\t'))
def loc(w):
    m = re.match(r'(.+):(\d+)-(\d+)([+-])', w)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None
def read(d, g):
    p = f'{d}/{g}.tsv'
    if not os.path.exists(p): return {}
    return {r['gene']: r for r in csv.DictReader(open(p), delimiter='\t')}
rows = []
for r in gates:
    if r['clade'] != 'clade': continue
    l = loc(r['scan_where'])
    fl = read(fd, r['genome']); fl2 = read(fd2, r['genome'])
    out = dict(genome=r['genome'], genus=r['genus'], species=r['species'], strain=r['strain'], gene=r['gene'], scan_score=r['scan_score'], stage=r['stage'],
               mat_contig=l[0], mat_start=l[1], mat_end=l[2], n50=r['n50'])
    clen = None
    for tag, d in (('', fl), ('2', fl2)):
        for gn in ('SLA2', 'APN2', 'COX13', 'GAP1', 'BGL1', 'ORF1'):
            h = d.get(gn)
            if not h: continue
            if clen is None and h['contig'] == l[0]: clen = int(h['contig_len'])
            same = h['contig'] == l[0]
            dist = min(abs(int(h['start']) - l[2]), abs(int(h['end']) - l[1])) if same else ''
            out[gn + '_id'] = h['pident']; out[gn + '_qcov'] = h['qcov']; out[gn + '_same'] = same; out[gn + '_dist'] = dist
    rows.append(out)
cols = ['genome', 'genus', 'species', 'strain', 'gene', 'scan_score', 'stage', 'mat_contig', 'mat_start', 'mat_end', 'n50'] + [f'{g}_{k}' for g in ('SLA2', 'APN2', 'COX13', 'GAP1', 'BGL1', 'ORF1') for k in ('id', 'qcov', 'same', 'dist')]
with open('results/2026-10-18_cochliobolus_recall/clade_flank.tsv', 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(rows)
def num(x):
    try: return float(x)
    except: return None
import statistics as st
for g in ('SLA2', 'APN2', 'COX13', 'GAP1', 'BGL1'):
    d = [num(r.get(g + '_dist')) for r in rows if r.get(g + '_dist') not in (None, '')]
    ids = [num(r.get(g + '_id')) for r in rows if r.get(g + '_id')]
    if ids:
        print(g, 'n with hit', len(ids), 'median id', st.median(ids), '| same contig', len(d), 'median dist', st.median(d) if d else None,
              'within25kb', sum(1 for x in d if x <= 25000), 'within100kb', sum(1 for x in d if x <= 100000), 'within300kb', sum(1 for x in d if x <= 300000))
print(len(rows), 'gene rows')
