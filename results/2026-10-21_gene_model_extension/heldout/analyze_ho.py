import pandas as pd, numpy as np
pd.set_option('display.width',200); pd.set_option('display.max_rows',100)
d=pd.read_csv('heldout_rescue.tsv',sep='\t'); print('rows',len(d)); print(d.model_found.value_counts().to_dict())
m=d[d.model_found=='yes'].copy()
for c in ['gap_before_aa','delta_aa','ext_aa','model_aa','q_end','truth_introns_after_edge','truth_aa']: m[c]=pd.to_numeric(m[c],errors='coerce')
tr=pd.read_csv('truth_heldout.tsv',sep='\t'); tr['key']=tr.record+'|'+tr.gene
m=m.merge(tr[['key','nexons','species']].rename(columns={'species':'sp'}),on='key',how='left')
print('strand mismatches (no end comparison):',int(m.gap_before_aa.isna().sum())); m=m[m.gap_before_aa.notna()]
t=m[(m.gap_before_aa>=30)&(m.ext_state=='stop')].copy()
print('genes with a model:',m.key.nunique(),' truncated rows reaching a stop:',len(t),' genes:',t.key.nunique())
# rules (frozen)
refs={};k=None
for l in open('refs.faa'):
    if l[0]=='>': k=l[1:].strip(); refs[k]=''
    else: refs[k]+=l.strip()
lens={}
for h,s in refs.items(): lens.setdefault(h.split('|')[1],[]).append(len(s))
rng={g:(np.percentile(v,10),np.percentile(v,90)) for g,v in lens.items()}
print('length ranges (10th-90th pct):',{g:tuple(round(x) for x in r) for g,r in rng.items()})
ri=pd.read_csv('ref_introns_all.tsv',sep='\t'); ri['best_ref']=ri.record
def refinfo(r):
    x=ri[(ri.record==r.best_ref)&(ri.gene==r.gene)]
    if len(x)==0: return pd.Series(dict(ref_known=False,ref_after=np.nan,ref_nexons=np.nan))
    s=x.iloc[0].intron_aa; ints=[float(v) for v in s.split(';')] if isinstance(s,str) else []
    return pd.Series(dict(ref_known=True,ref_after=sum(v>r.q_end+2 for v in ints),ref_nexons=int(x.iloc[0].nexons)))
t=pd.concat([t,t.apply(refinfo,axis=1)],axis=1)
t['total']=t.model_aa+t.ext_aa; t['inrange']=[(a>=rng[g][0] and a<=rng[g][1]) for a,g in zip(t.total,t.gene)]
t['len_rule']=(t.ext_aa>=200)|((t.ext_aa>=30)&(t.ext_aa<200)&t.inrange)
t['guard_clear']=t.ref_known&(t.ref_after==0)
t['good']=t.delta_aa.abs()<=5; t['early']=t.delta_aa<-20; t['over5']=t.delta_aa>5
def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; dd=1+z*z/n; c=(p+z*z/(2*n))/dd; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/dd; return max(0,c-h),min(1,c+h)
def stat(name,mask):
    g=t[mask]; n=len(g); k=int(g.early.sum()); lo,hi=wilson(k,n); ng=g.key.nunique(); nk=g[g.early].key.nunique()
    # per-gene rate (a gene counts as failing if any of its accepted rows is an early stop)
    pg=g.groupby('key').early.any(); lo2,hi2=wilson(int(pg.sum()),len(pg))
    print(f'{name:62s} rows {n:3d} early {k:2d} ({100*k/max(n,1):4.1f}%, CI {100*lo:.0f}-{100*hi:.0f}%)  genes {ng:3d} genes-with-early {nk:2d} per-gene {100*pg.mean() if len(pg) else 0:4.1f}% (CI {100*lo2:.0f}-{100*hi2:.0f}%)  orders {g.order.nunique():2d}  correct(<=5aa) {int(g.good.sum()):3d}  overshoot>5 {int(g.over5.sum())}')
print()
stat('E every extension >= 30 aa',t.ext_aa>=30)
stat('A length rule only',t.len_rule)
stat('B length rule + guard clear',t.len_rule&t.guard_clear)
stat("D' B + best-hit reference has <=2 exons",t.len_rule&t.guard_clear&(t.ref_nexons<=2))
stat('F  B + best-hit reference has >=3 exons',t.len_rule&t.guard_clear&(t.ref_nexons>=3))
print('\noverall: exact within 5 aa',int(t.good.sum()),'of',len(t),' early',int(t.early.sum()),' overshoot>5',int(t.over5.sum()),' max overshoot',t.delta_aa.max())
print('by deposited exon count (all truncated rows):'); print(t.groupby('nexons').apply(lambda g: pd.Series(dict(genes=g.key.nunique(),rows=len(g),good=int(g.good.sum()),early=int(g.early.sum())))))
print('\nD-prime rows by order:'); dp=t[t.len_rule&t.guard_clear&(t.ref_nexons<=2)]; print(dp.groupby('order').apply(lambda g: pd.Series(dict(rows=len(g),genes=g.key.nunique(),early=int(g.early.sum())))).to_string())
print('\nguard on all rows with known ref:'); k=t[t.ref_known]; fl=k.ref_after>=1
print(' flagged',int(fl.sum()),'early among flagged',int((fl&k.early).sum()),' not flagged early',int((~fl&k.early).sum()),'not flagged good',int((~fl&k.good).sum()),'of',int((~fl).sum()))
t.to_csv('heldout_truncated_rows.tsv',sep='\t',index=False)
# near-complete alignments
n=m[(m.frac==1.0)&(m.gap_before_aa<30)]; print('\nnear-complete alignments (100% refs, gap<30 aa): rows',len(n),' extended (ext>0):',int((n.ext_aa>0).sum()),' exact/within5:',int((n.delta_aa.abs()<=5).sum()))
