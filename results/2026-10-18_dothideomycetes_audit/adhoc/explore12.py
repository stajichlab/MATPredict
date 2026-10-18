import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100)
D = sys.argv[1]
L = pd.read_csv(D + '/loci_audit_flagged.tsv', sep='\t', low_memory=False)
c = L[L.called]
x = c[c.prof_order_call.notna()]
print('order-out profile scored loci', len(x), x.order.value_counts().to_dict())
print(pd.crosstab(x.order, x.prof_order_call.where(x.prof_order_call == x.idiomorph, 'other:' + x.prof_order_call.astype(str))))
print('genus-out used:', c.prof_set.astype(str).str.startswith('genus:').sum(), c[c.prof_set.astype(str).str.startswith('genus:')].prof_set.value_counts().to_dict())
g = c[c.prof_set.astype(str).str.startswith('genus:')]
print(pd.crosstab(g.prof_set, g.vote_vs_profile))
print('full-profile loci', (c.prof_set == 'full').sum())
