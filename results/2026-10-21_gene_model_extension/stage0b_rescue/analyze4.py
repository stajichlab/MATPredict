import pandas as pd, numpy as np
pd.set_option('display.width',220); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',100)
d=pd.read_csv('s0b_rescue2.tsv',sep='\t'); m=d[d.model_found=='yes'].copy()
for c in ['gap_before_aa','delta_aa','ext_aa','model_aa','q_end','truth_introns_after_edge']: m[c]=pd.to_numeric(m[c],errors='coerce')
ri=pd.read_csv('ref_introns_all.tsv',sep='\t').set_index('record')
def ref_after(r):
    if r.best_ref not in ri.index: return np.nan
    s=ri.loc[r.best_ref,'intron_aa']
    if not isinstance(s,str): return 0
    return sum(float(x)>r.q_end+2 for x in s.split(';'))
m['ref_known']=m.best_ref.isin(ri.index)
m['ref_introns_after']=m.apply(ref_after,axis=1)
t=m[(m.gap_before_aa>=30)&(m.ext_state=='stop')].copy()
t['good']=t.delta_aa.abs()<=5
t['under']=t.delta_aa<-20
print('truncated rows with a stop:',len(t),' best reference has exon info:',int(t.ref_known.sum()),f'({100*t.ref_known.mean():.0f}%)')
print('\n--- ORACLE: true introns after the edge in the deposited gene vs outcome ---')
print(pd.crosstab(t.truth_introns_after_edge,[t.good]).rename(columns={True:'exact/within5',False:'not'}))
print(pd.crosstab(t.truth_introns_after_edge,t.under).rename(columns={True:'undershoot',False:'not undershoot'}))
print('\n--- GUARD: introns the best-hit REFERENCE has after the edge (rows whose reference has exon info) ---')
k=t[t.ref_known]
print(pd.crosstab(k.ref_introns_after,k.good).rename(columns={True:'exact/within5',False:'not'}))
print(pd.crosstab(k.ref_introns_after,k.under).rename(columns={True:'undershoot',False:'not undershoot'}))
k=k.assign(flag=k.ref_introns_after>=1)
tp=int((k.flag&k.under).sum()); fp=int((k.flag&~k.under).sum()); fn=int((~k.flag&k.under).sum()); tn=int((~k.flag&~k.under).sum())
print(f'\nguard flag = reference has >=1 intron after the edge, target = extension undershoots by >20 aa  (n={len(k)})')
print(f' flagged & undershoot={tp}  flagged & fine={fp}  not flagged & undershoot={fn}  not flagged & fine={tn}')
print(f' precision={tp/(tp+fp):.2f}  recall={tp/(tp+fn):.2f}  | not-flagged rows exact/within5: {int((~k.flag&k.good).sum())} of {int((~k.flag).sum())}')
print(f' flagged rows exact/within5 (would be held back though correct): {int((k.flag&k.good).sum())} of {int(k.flag.sum())}')
# agreement between reference-based and truth-based counts
print('\nreference intron count after edge vs truth intron count after edge:'); print(pd.crosstab(k.ref_introns_after,k.truth_introns_after_edge))
# combine with the length screen from before
refs={};kk=None
for l in open('refs.faa'):
    if l[0]=='>': kk=l[1:].strip(); refs[kk]=''
    else: refs[kk]+=l.strip()
lens={}
for h,s in refs.items(): lens.setdefault(h.split('|')[1],[]).append(len(s))
rng={g:(np.percentile(v,10),np.percentile(v,90)) for g,v in lens.items()}
t['total']=t.model_aa+t.ext_aa
t['inrange']=[(a>=rng[g][0] and a<=rng[g][1]) for a,g in zip(t.total,t.gene)]
t['len_rule']=((t.ext_aa>=200)|((t.ext_aa>=30)&(t.ext_aa<200)&t.inrange))
t['guard_ok']=np.where(t.ref_known,t.ref_introns_after==0,np.nan)
print('\n--- COMBINED: length rule AND (reference known with no intron after edge) ---')
for name,mask in [('length rule only',t.len_rule),('length rule & reference known & no intron after edge',t.len_rule&(t.guard_ok==1)),
                  ('length rule & NOT guard-clear (ref known, intron after, or unknown)',t.len_rule&~(t.guard_ok==1))]:
    x=t[mask]; print(f'{name:70s} n={len(x):3d} good={int(x.good.sum()):3d} undershoot={int(x.under.sum()):3d} overshoot>5={int((x.delta_aa>5).sum())}')
print('\nreference-unknown rows:',int((~t.ref_known).sum()),'  of those good',int(t[~t.ref_known].good.sum()),'undershoot',int(t[~t.ref_known].under.sum()))
t.to_csv('truncated_rows4.tsv',sep='\t',index=False)
