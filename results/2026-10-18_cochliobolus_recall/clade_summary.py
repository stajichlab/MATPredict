import csv, sys
from collections import Counter, defaultdict
rows = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_classify.tsv'), delimiter='\t'))
def f(x):
    try: return float(x)
    except: return None
def cls(r):
    i = f(r['best_core_id'])
    if i is None: return 'no_core_hit'
    if i >= 85: return 'core>=85'
    if i >= 70: return 'core70-85'
    if i >= 55: return 'core55-70'
    return 'core<55'
tab = defaultdict(Counter)
for r in rows:
    tab[(r['clade'], r['genus'])][(r['status'], cls(r))] += 1
for k in sorted(tab):
    print(k, sum(tab[k].values()), dict(tab[k]))
print()
print('clade no-calls with core>=85 and what stopped them')
for r in rows:
    if r['clade'] == 'clade' and r['status'] != 'called':
        i = f(r['best_core_id'])
        if i is not None and i >= 70:
            print(r['genus'], r['species'][:28], r['strain'][:12], r['genome'][:15], 'id', r['best_core_id'], 'genes', r['best_core_genes'], 'adm', r['best_core_admitted'],
                  'len', r['best_core_len'], 'n50', r['n50'], 'flank', r['flank_id'], r['flank_dist'])
print()
for r in rows:
    if r['clade'] == 'clade' and r['status'] == 'called':
        print('CALLED', r['genome'], r['species'], r['best_core_id'], r['best_core_genes'], r['best_core_admitted'])
