"""Control sample for the roster-prototype: called genomes outside the clade (seed 7): 30 Pleosporales (non-clade) and 8 from each of the other 8 largest orders."""
import csv, random
from collections import defaultdict
random.seed(7)
rows = list(csv.DictReader(open('results/2026-10-18_dothideomycetes_audit/genome_audit_final.tsv'), delimiter='\t'))
clade = {'Bipolaris', 'Curvularia', 'Exserohilum', 'Ophiosphaerella'}
by = defaultdict(list)
for r in rows:
    if r['status'] == 'called' and r['genus'] not in clade:
        by[r['order']].append(r)
orders = sorted(by, key=lambda o: -len(by[o]))
pick = []
for o in orders[:9]:
    n = 30 if o == 'Pleosporales' else 8
    pick += random.sample(by[o], min(n, len(by[o])))
with open('results/2026-10-18_cochliobolus_recall/control_names_taxid.tsv', 'w') as f:
    for r in pick:
        f.write(f"{r['genome']}\t{r['taxid']}\n")
print(len(pick), [(o, len(by[o])) for o in orders[:9]])
