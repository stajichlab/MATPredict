import csv, collections, statistics as st
R=list(csv.DictReader(open('gene_audit.tsv'),delimiter='\t'))
for r in R:
    for k in ('identity','bitscore','alpha','hmg','aa'): r[k]=float(r[k]) if r[k] not in ('','None') else 0.0
    r['modelled']=r['modelled']=='True'
C=collections.Counter((r['run'],r['gene']) for r in R); print(sorted(C.items()))
print('lost_to', collections.Counter((r['run'],r['gene']) for r in R if r['lost_to'] not in ('','None')))
print('refs'); 
for k,v in sorted(collections.Counter((r['run'],r['gene'],r['ref']) for r in R).items()): print(' ',k,v)
def q(v): v=sorted(v); return [round(v[int(p*(len(v)-1))],1) for p in (0,.1,.25,.5,.75,.9,1)] if v else []
print()
for run in ('v060','records'):
  for g in ('MAT1-1-1','MAT1-2-1','MAT1-1-3','MAT1-1-4'):
    for mod in (True,False):
      s=[r for r in R if r['run']==run and r['gene']==g and r['modelled']==mod]
      if not s: continue
      print(run,g,'modelled' if mod else 'unmodelled',len(s),'aa',q([r['aa'] for r in s]),'alpha',q([r['alpha'] for r in s]),'hmg',q([r['hmg'] for r in s]),'bits',q([r['bitscore'] for r in s]),'id',q([r['identity'] for r in s]))
