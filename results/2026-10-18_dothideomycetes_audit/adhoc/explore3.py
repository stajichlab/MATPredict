import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 200)
D = sys.argv[1]
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
G = G[(G.gw_scanned.astype(str) == 'True') & G.primary_idiomorph.notna()].copy()
G['hmg_at'] = pd.to_numeric(G.gw_hmg_best_at_locus, errors='coerce')
G['alpha_at'] = pd.to_numeric(G.gw_alpha_best_at_locus, errors='coerce')
m2 = G[G.primary_idiomorph == 'MAT1-2']
t = m2.groupby('genus').agg(n=('genome', 'size'), hmg_orf_at_locus=('hmg_at', lambda s: s.notna().sum()), q10=('hmg_at', lambda s: s.quantile(.1)),
                            med=('hmg_at', 'median'), q90=('hmg_at', lambda s: s.quantile(.9)))
print(t.sort_values('n', ascending=False).head(25))
m1 = G[G.primary_idiomorph == 'MAT1-1']
t = m1.groupby('genus').agg(n=('genome', 'size'), alpha_orf_at_locus=('alpha_at', lambda s: s.notna().sum()), q10=('alpha_at', lambda s: s.quantile(.1)),
                            med=('alpha_at', 'median'), q90=('alpha_at', lambda s: s.quantile(.9)))
print(t.sort_values('n', ascending=False).head(25))
