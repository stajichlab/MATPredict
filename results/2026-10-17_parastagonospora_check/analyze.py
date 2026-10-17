import csv, collections, gzip, os
T=lambda f: list(csv.DictReader(open(f),delimiter='\t'))
loci=T('locus_summary.tsv'); aud=T('gene_audit.tsv')
asm={r['acc']:r for r in T('para_asm.tsv')}; bs={r['biosample']:r for r in T('para_bs.tsv')}
sp={r['genome']:r['species'] for r in T('../2026-10-03_ascomycota_v060/genomes.tsv') if r['species'].startswith('Parastagonospora')}
def f(x):
    try: return float(x)
    except: return 0.0
# per-run locus counts
for run in ('v060','records'):
    L=[r for r in loci if r['run']==run]
    print(run,'genomes',len({r['genome'] for r in L}),'with locus',len({r['genome'] for r in L if r['nloci']!='0'}),'multi-locus genomes',len({r['genome'] for r in L if int(r['nloci'])>1}))
    print(' labels',collections.Counter(r['label'] for r in L if r['nloci']!='0'),'class',collections.Counter(r['lclass'] for r in L if r['nloci']!='0'),'conf',collections.Counter(r['conf'] for r in L if r['nloci']!='0'))
# resolutions
for run in ('v060','records'):
    print(run,'resolutions',collections.Counter(r['res'].split('(')[0]+'('+r['res'].split('(')[1].split(',')[0] if r['res'] else '' for r in loci if r['run']==run and r['nloci']!='0'))
    print(run,'candidates', collections.Counter(r['cands'].count(':') for r in loci if r['run']==run and r['nloci']!='0'))
# strong gene classification in records run
G=collections.defaultdict(list)
for r in aud:
    for k in ('identity','bitscore','alpha','hmg','aa'): r[k]=f(r[k])
    G[(r['run'],r['genome'],r['locus'])].append(r)
def own(r): return r['alpha'] if r['gene']=='MAT1-1-1' else r['hmg'] if r['gene']=='MAT1-2-1' else None
def strong(r,gn): 
    if r['gene']!=gn: return False
    return r['modelled'] in ('True',True) and r['aa']>=300 and own(r)>=150
res=collections.Counter(); both=[]
rows=[]
lab={}
for (run,g,li),es in G.items():
    lbl=es[0]['label']
    a=[e for e in es if e['gene']=='MAT1-1-1']; h=[e for e in es if e['gene']=='MAT1-2-1']
    sa=any(strong(e,'MAT1-1-1') for e in a); sh=any(strong(e,'MAT1-2-1') for e in h)
    wa=[e for e in a if not strong(e,'MAT1-1-1')]; wh=[e for e in h if not strong(e,'MAT1-2-1')]
    cls=('both_strong' if sa and sh else 'alpha_strong' if sa else 'hmg_strong' if sh else 'neither_strong')
    res[(run,lbl,cls,'weak_alpha' if wa else '-', 'weak_hmg' if wh else '-')]+=1
    rows.append((run,g,li,lbl,cls))
    if run=='records': lab[g]=(lbl,cls)
for k,v in sorted(res.items()): print(k,v)
# weak cross-hit stats in records
for gn,ownf in (('MAT1-1-1','alpha'),('MAT1-2-1','hmg')):
    w=[e for e in aud if e['run']=='records' and e['gene']==gn and not strong(e,gn)]
    print('records weak',gn,len(w),'modelled',sum(e['modelled'] in ('True',True) for e in w),'aa max',max(e['aa'] for e in w),'bits max',max(e['bitscore'] for e in w),'own max',max(own(e) for e in w), 'other-profile max', max((e['hmg'] if gn=='MAT1-1-1' else e['alpha']) for e in w), 'ids',min(e['identity'] for e in w),max(e['identity'] for e in w))
    big=[e for e in w if own(e)>=30 or e['bitscore']>=36]
    for e in big: print('  notable',e['genome'],e['gene'],e['status'],e['aa'],e['bitscore'],e['alpha'],e['hmg'],e['identity'])
# genome-level summary
cnt=collections.Counter(v for v in lab.values()); print('records label x class',cnt)
# species x label
c=collections.Counter((sp[g],lab.get(g,('none',''))[0]) for g in sp); print(c)
# export genome table
meta={}
for g in sp:
    a=asm.get('_'.join(g.split('_')[:2]),{}); b=bs.get(a.get('biosample'),{})
    meta[g]=dict(strain=b.get('strain') or b.get('isolate') or a.get('strain') or a.get('isolate'),geo=b.get('geo_loc_name'),host=b.get('host'),date=b.get('collection_date'),biosample=a.get('biosample'),bioproject=a.get('bioproject'),n50=a.get('contigs'))
# mash clusters
d={}
for l in gzip.open('para_mash.tsv.gz','rt'):
    x,y,z,*_=l.split('\t'); d[(os.path.basename(x)[:-6],os.path.basename(y)[:-6])]=float(z)
gs=sorted({a for a,b in d}); par={x:x for x in gs}
def fnd(x):
    while par[x]!=x: par[x]=par[par[x]]; x=par[x]
    return x
T0=0.001
for (a,b),z in d.items():
    if a<b and z<T0: par[fnd(a)]=fnd(b)
cl=collections.defaultdict(list)
for g in gs: cl[fnd(g)].append(g)
print('clusters at mash<%s: %d from %d; sizes'%(T0,len(cl),len(gs)),sorted((len(v) for v in cl.values()),reverse=True)[:12])
mixed=0; clone_rows=[]
for k,v in cl.items():
    ls={lab.get(g,('none',''))[0] for g in v}
    if len(v)>1 and len(ls-{'none'})>1: mixed+=1
    for g in v: clone_rows.append((g,k,len(v)))
print('multi-member clusters with both labels:',mixed)
# representative per cluster: first by assembly N50 not needed; take sorted first with a call
rep={}
for k,v in cl.items():
    called=[g for g in sorted(v) if lab.get(g,('none',''))[0] in ('MAT1-1','MAT1-2')]
    rep[k]=called[0] if called else None
rc=collections.Counter(lab[g][0] for g in rep.values() if g); print('clone-corrected (one per cluster):',rc, 'n',sum(rc.values()))
# per species and region
def region(g):
    s=meta[g]['geo'] or ''
    return 'WA' if 'Western Australia' in s else (s.split(':')[0] or 'unknown')
tab=collections.Counter((sp[g],region(g),lab.get(g,('none',''))[0]) for g in sp)
for k,v in sorted(tab.items()): print(k,v)
# clone-corrected per region/species
tabc=collections.Counter((sp[g],region(g),lab[g][0]) for g in rep.values() if g)
print('clone-corrected by region'); 
for k,v in sorted(tabc.items()): print(k,v)
# per-genome table
with open('genomes_table.tsv','w') as fh:
    w=csv.writer(fh,delimiter='\t'); w.writerow(['genome','species','strain','region','geo','host','date','biosample','bioproject','records_label','locus_class','v060_label','mash_cluster','cluster_size'])
    v0={r['genome']:r['label'] for r in loci if r['run']=='v060' and r['nloci']!='0'}
    cm={g:(k,n) for g,k,n in clone_rows}
    for g in sorted(sp):
        m=meta[g]; w.writerow([g,sp[g],m['strain'],region(g),m['geo'],m['host'],m['date'],m['biosample'],m['bioproject'],lab.get(g,('none',''))[0],lab.get(g,('',''))[1],v0.get(g,'none'),cm[g][0][:15],cm[g][1]])
