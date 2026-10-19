import csv
from collections import Counter
rows = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_joined.tsv'), delimiter='\t'))
vals = Counter()
for r in rows:
    if r['clade'] != 'clade': continue
    print(r['genus'][:5], r['species'][-14:], r['strain'][:9], r['status'][:3], 'a', r['alpha_prof_best'], 'apf', r['alpha_pfam_best'], 'h', r['hmg_prof_best'], 'cl', r['best_core_id'], r['best_core_genes'], r['cluster_is_scan_orf'], r['n_flank70_within300kb'])
