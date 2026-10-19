import csv
from collections import defaultdict
rows = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_flank.tsv'), delimiter='\t'))
def num(x):
    try: return float(x)
    except: return None
t = defaultdict(lambda: defaultdict(int))
for r in rows:
    g = r['genus']; t[g]['n'] += 1
    gap = r['GAP1_same'] == 'True' and (num(r['GAP1_id']) or 0) >= 50 and (num(r['GAP1_dist']) or 1e9) <= 25000
    bgl = r['BGL1_same'] == 'True' and (num(r['BGL1_id']) or 0) >= 50 and (num(r['BGL1_dist']) or 1e9) <= 25000
    t[g]['gap'] += gap; t[g]['bgl'] += bgl; t[g]['either'] += (gap or bgl); t[g]['both'] += (gap and bgl)
    r['_gap'] = gap; r['_bgl'] = bgl
for g, d in t.items(): print(g, dict(d))
# genes whose neither flank is on the MAT contig within 25 kb
for r in rows:
    if not (r['_gap'] or r['_bgl']):
        print('NO ADJ', r['genus'], r['species'][-16:], r['strain'][:10], r['gene'], r['stage'], 'n50', r['n50'], 'GAP', r['GAP1_id'], r['GAP1_same'], r['GAP1_dist'], 'BGL', r['BGL1_id'], r['BGL1_same'], r['BGL1_dist'])
