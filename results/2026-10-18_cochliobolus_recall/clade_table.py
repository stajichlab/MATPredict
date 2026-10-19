import csv
from collections import Counter
r=list(csv.DictReader(open('results/2026-10-18_dothideomycetes_audit/genome_audit_final.tsv'),delimiter='\t'))
GEN=('Bipolaris','Curvularia','Exserohilum','Ophiosphaerella','Cochliobolus','Setosphaeria')
fam=sorted(set(x['family'] for x in r if x['genus'] in GEN))
print(fam)
c=Counter((x['genus'],x['status'],x['idiomorphs_called']) for x in r if x['genus'] in GEN)
for k,v in sorted(c.items()):print(k,v)
print(list(r[0].keys()))
for x in r:
  if x['genus'] in GEN and (x['n_detected']!='0' or 'C4' in x['genome'] or 'C5' in x['genome']):
    print({k:x[k] for k in ('genome','status','n_detected','idiomorphs','suppressed_unpolished','verdict','gw_class','primary_win_gene')})
