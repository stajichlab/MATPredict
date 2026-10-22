import pandas as pd, numpy as np
pd.set_option('display.width',220); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',100)
t=pd.read_csv('truncated_rows.tsv',sep='\t')
refs={};k=None
for l in open('refs.faa'):
    if l[0]=='>': k=l[1:].strip(); refs[k]=''
    else: refs[k]+=l.strip()
lens={}
for h,s in refs.items(): lens.setdefault(h.split('|')[1],[]).append(len(s))
rng={g:(np.percentile(v,10),np.percentile(v,90),np.median(v),max(v)) for g,v in lens.items()}
print({g:tuple(round(x) for x in r) for g,r in rng.items()})
t=t[t.ext_state=='stop'].copy()
t['total_aa']=t.model_aa+t.ext_aa
t['lo']=t.gene.map(lambda g:rng[g][0]); t['hi']=t.gene.map(lambda g:rng[g][1]); t['maxlen']=t.gene.map(lambda g:rng[g][3])
t['in_range']=(t.total_aa>=t.lo)&(t.total_aa<=t.hi)
t['below']=t.total_aa<t.lo; t['above_guard']=t.total_aa>1.25*t.maxlen
t['outcome']=np.where(t.exact,'exact',np.where(t.delta_aa.abs()<=5,'within5',np.where(t.delta_aa<-20,'undershoot','other')))
print('\noutcome by in-range / below the 10th percentile of curated lengths (all stop-found truncated rows):')
print(pd.crosstab(t.outcome,[t.below]))
print('\nof rows total length inside [10th,90th] percentile:',int(t.in_range.sum()),'; outcome:'); print(t[t.in_range].outcome.value_counts().to_dict())
print('rows below range:',int(t.below.sum()),'; outcome:',t[t.below].outcome.value_counts().to_dict())
print('rows above 1.25x longest:',int(t.above_guard.sum()))
t['tier']=pd.cut(t.ext_aa,[-1,29,199,10000],labels=['<30','30-199','>=200'])
print('\nTIER x in_range -> outcome counts'); print(pd.crosstab([t.tier,t.in_range],t.outcome))
# rule: accept as complete_orf_supported = ext>=200 or (30<=ext<200 and total in range)
t['rule_accept']=((t.ext_aa>=200)&(t.total_aa<=t.hi*1.25))|((t.ext_aa>=30)&(t.ext_aa<200)&t.in_range)
print('\nrule (>=200 aa ORF-only, or 30-199 aa AND total length in curated 10-90th pct):')
print(pd.crosstab(t.rule_accept,t.outcome)); 
a=t[t.rule_accept]; print('accepted:',len(a),'exact or within5:',int(a.outcome.isin(['exact','within5']).sum()),'undershoot:',int((a.outcome=='undershoot').sum()),'overshoot>20:',int((a.delta_aa>20).sum()))
r=t[~t.rule_accept]; print('provisional/rejected:',len(r),'of which exact/within5:',int(r.outcome.isin(['exact','within5']).sum()),'undershoot:',int((r.outcome=='undershoot').sum()))
# by frac
print('\naccepted-rule performance by frac:'); print(t.groupby('frac').apply(lambda g: pd.Series(dict(n=len(g),accepted=int(g.rule_accept.sum()),acc_good=int((g.rule_accept&g.outcome.isin(['exact','within5'])).sum()),acc_bad=int((g.rule_accept&(g.outcome=='undershoot')).sum())))))
# by order
print('\nexact rate by order (all truncated stop-found):'); print(t.groupby('order').apply(lambda g: pd.Series(dict(n=len(g),exact=int(g.exact.sum()),undershoot=int((g.outcome=='undershoot').sum())))))
# does truth have an intron after the model end? (nexons vs undershoot)
t.to_csv('truncated_rows2.tsv',sep='\t',index=False)
