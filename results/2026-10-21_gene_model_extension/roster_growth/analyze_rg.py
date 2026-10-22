import pandas as pd, numpy as np, glob
S='/private/tmp/claude-501/-Users-jstajich-projects-MATPredict/519c0d9b-c67c-487f-bb14-2f7c2115d61e/scratchpad'
def load(path_list,radius=None):
    d=pd.concat([pd.read_csv(p,sep='\t') for p in path_list],ignore_index=True); d=d[d.model_found=='yes'].copy()
    if radius: d=d[d.radius==radius]
    for c in ['gap_before_aa','delta_aa','ext_aa','model_aa','q_end','truth_introns_after_edge']: d[c]=pd.to_numeric(d[c],errors='coerce')
    return d
base=pd.concat([load([f'{S}/s0b/s0b_rescue2.tsv'],'genus'),load([f'{S}/heldout/heldout_rescue.tsv'])],ignore_index=True)
aug=load(sorted(glob.glob('rg_rescue_*.tsv')))
tr=pd.concat([pd.read_csv(f'{S}/s0b/truth.tsv',sep='\t'),pd.read_csv(f'{S}/heldout/truth_heldout.tsv',sep='\t')],ignore_index=True); tr['key']=tr.record+'|'+tr.gene
nex=tr.set_index('key').nexons
ri_base=pd.read_csv(f'{S}/s0b/ref_introns_all.tsv',sep='\t'); ri_aug=pd.read_csv('ref_introns_aug.tsv',sep='\t')
refs={};k=None
for l in open(f'{S}/s0b/refs.faa'):
    if l[0]=='>': k=l[1:].strip(); refs[k]=''
    else: refs[k]+=l.strip()
lens={}
for h,s in refs.items(): lens.setdefault(h.split('|')[1],[]).append(len(s))
rng={g:(np.percentile(v,10),np.percentile(v,90)) for g,v in lens.items()}
def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; dd=1+z*z/n; c=(p+z*z/(2*n))/dd; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/dd; return max(0,c-h),min(1,c+h)
def prep(d,ri):
    d=d.copy(); d['nexons']=d.key.map(nex); d=d[(d.gap_before_aa>=30)&(d.ext_state=='stop')].copy()
    idx={(r.record,r.gene):r for r in ri.itertuples()}
    def f(r):
        x=idx.get((r.best_ref,r.gene))
        if x is None: return pd.Series(dict(ref_known=False,ref_after=np.nan,ref_nexons=np.nan))
        ints=[float(v) for v in x.intron_aa.split(';')] if isinstance(x.intron_aa,str) else []
        return pd.Series(dict(ref_known=True,ref_after=sum(v>r.q_end+2 for v in ints),ref_nexons=int(x.nexons)))
    d=pd.concat([d,d.apply(f,axis=1)],axis=1)
    d['total']=d.model_aa+d.ext_aa; d['inrange']=[(a>=rng[g][0] and a<=rng[g][1]) for a,g in zip(d.total,d.gene)]
    d['len_rule']=(d.ext_aa>=200)|((d.ext_aa>=30)&(d.ext_aa<200)&d.inrange); d['guard_clear']=d.ref_known&(d.ref_after==0)
    d['good']=d.delta_aa.abs()<=5; d['early']=d.delta_aa<-20; d['over5']=d.delta_aa>5; d['bad']=d.early|d.over5
    d['y']=(d.truth_introns_after_edge>=1); return d
b=prep(base,ri_base); a=prep(aug,ri_aug)
common=set(zip(b.key,b.frac))&set(zip(a.key,a.frac)); b=b[[ (k,f) in common for k,f in zip(b.key,b.frac)]]; a=a[[ (k,f) in common for k,f in zip(a.key,a.frac)]]
print('rows compared (same genes and truncations):',len(b),' genes',b.key.nunique())
def guard(d,name):
    k=d[d.ref_known]; fl=k.ref_after>=1; y=k.y
    tp=int((fl&y).sum()); fp=int((fl&~y).sum()); fn=int((~fl&y).sum()); tn=int((~fl&~y).sum())
    print(f'{name:34s} reference known {len(k)}/{len(d)}; predicting "intron after edge": flagged {int(fl.sum())}, true positives {tp}, false positives {fp}, missed {fn}  precision {tp/max(tp+fp,1):.2f} recall {tp/max(tp+fn,1):.2f}')
    m=(k.ref_nexons==k.nexons); print(f'{"":34s} best-hit reference exon count equals the deposited gene exon count in {int(m.sum())} of {len(k)} rows ({100*m.mean():.0f}%)')
    from sklearn.metrics import roc_auc_score
    try: print(f'{"":34s} AUC of the flag for "intron after edge": {roc_auc_score(y,fl):.3f}')
    except Exception as e: pass
guard(b,'baseline (55 references)'); guard(a,'augmented (161 references)')
def stat(name,d,mask):
    g=d[mask]; n=len(g); bad=int(g.bad.sum()); e=int(g.early.sum()); o=int(g.over5.sum()); lo,hi=wilson(bad,n); ng=g.key.nunique(); pg=g.groupby('key').bad.any()
    print(f'  {name:34s} accepted {n:3d}  early {e:2d} overshoot {o:2d}  wrong {100*bad/max(n,1):4.1f}% (CI {100*lo:.0f}-{100*hi:.0f}%)  genes {ng:3d} genes-wrong {int(pg.sum()):2d} ({100*pg.mean() if len(pg) else 0:4.1f}%)  correct {int(g.good.sum())}')
for nm,d in (('BASELINE (55 refs, leave-genus-out)',b),('AUGMENTED (161 refs, leave-genus-out)',a)):
    print('\n'+nm)
    stat('A length rule',d,d.len_rule); stat('B length rule + guard clear',d,d.len_rule&d.guard_clear); stat("D' B + ref <= 2 exons",d,d.len_rule&d.guard_clear&(d.ref_nexons<=2)); stat('F  B + ref >= 3 exons',d,d.len_rule&d.guard_clear&(d.ref_nexons>=3))
print('\nbest-hit reference now comes from the added deposits in',int(a.best_ref.astype(str).str.startswith('HO_').sum()),'of',len(a),'rows')
a.to_csv('augmented_rows.tsv',sep='\t',index=False); b.to_csv('baseline_rows.tsv',sep='\t',index=False)
