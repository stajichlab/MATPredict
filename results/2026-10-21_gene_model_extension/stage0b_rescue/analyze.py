import pandas as pd, numpy as np
pd.set_option('display.width',220); pd.set_option('display.max_columns',30); pd.set_option('display.max_rows',200)
d=pd.read_csv('s0b_rescue.tsv',sep='\t')
print('rows',len(d)); print(d.model_found.value_counts().to_dict())
m=d[d.model_found=='yes'].copy()
for c in ['gap_before_aa','delta_aa','ext_aa','model_aa','truth_aa']: m[c]=pd.to_numeric(m[c],errors='coerce')
m['exact']=m.exact.astype(str)=='True'
m['strand_ok']=m.gap_before_aa.notna()
print('\nmodels found per frac/radius'); print(m.groupby(['radius','frac']).size().unstack())
print('strand mismatch (excluded from end comparison):',int((~m.strand_ok).sum()))
m=m[m.strand_ok]
# alignment-truncated: model ends >=30 aa before the true stop
m['trunc']=m.gap_before_aa>=30
print('\nalignment ends >=30 aa before the true end, by frac/radius:'); print(m.groupby(['radius','frac']).trunc.agg(['sum','count']).unstack(0))
t=m[m.trunc].copy()
t['ext_found']=t.ext_state=='stop'
t['abs_delta']=t.delta_aa.abs()
def summ(g):
    return pd.Series(dict(n=len(g),stop_found=int(g.ext_found.sum()),exact=int(g.exact.sum()),within5=int((g.abs_delta<=5).sum()),within20=int((g.abs_delta<=20).sum()),
        overshoot_gt20=int((g.delta_aa>20).sum()),undershoot_lt_m20=int((g.delta_aa<-20).sum()),median_delta=g.delta_aa.median()))
print('\n=== truncated alignments: does extension to the next stop recover the true end? ==='); print(t.groupby(['radius','frac']).apply(summ).round(1))
print('\nall truncated pooled:'); print(summ(t).round(1).to_dict())
# by extension length tier
t['tier']=pd.cut(t.ext_aa,[-1,29,199,10000],labels=['<30','30-199','>=200'])
print('\nby extension length tier (truncated rows with a stop):'); print(t[t.ext_found].groupby('tier').apply(summ).round(1))
# donor
print('\nby donor GT at model end:'); print(t[t.ext_found].groupby('donor_GT').apply(summ).round(1))
# per truth intron-after-edge: truth gene nexons
tr=pd.read_csv('truth.tsv',sep='\t')[['record','gene','nexons']]; t['record']=t.key.str.split('|').str[0]
t=t.merge(tr,on=['record','gene'],how='left')
# false extension on complete models (frac 1.0, alignment reaches the end)
c=m[(m.frac==1.0)&(m.gap_before_aa<30)]
print('\n=== complete / near-complete alignments (frac 1.0, gap<30 aa) ===')
print(c.groupby('radius').apply(lambda g: pd.Series(dict(n=len(g),closed_by_alignment=int((g.closed_by_alignment==True).sum()),ext_stop_found=int((g.ext_state=='stop').sum()),extended_gt0=int((g.ext_aa>0).sum()),exact=int(g.exact.sum()),
     delta_gt5=int((g.delta_aa.abs()>5).sum())))))
t.to_csv('truncated_rows.tsv',sep='\t',index=False)
