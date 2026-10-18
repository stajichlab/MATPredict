import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100); pd.set_option('display.max_colwidth', 80)
D = sys.argv[1]
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
c = G[G.primary_idiomorph.notna()]
print(pd.crosstab(c.primary_idiomorph, c.domain_final))
print(c[c.domain_final == 'no'][['genome', 'species', 'primary_idiomorph', 'primary_confidence', 'primary_class', 'primary_own_id', 'primary_own_cov', 'verdict']].to_string())
print(pd.crosstab(c.primary_vs_profile, c.primary_idiomorph))
print(c[c.primary_vs_profile == 'profile_undetermined'].groupby('genus').size().sort_values(ascending=False).head(8))
print('scan second_gene by verdict'); print(pd.crosstab(c.second_gene_scan, c.verdict_class))
print('flip: v0.6.0 vs now among genomes in both');
x = G[G.in_v060 == True]
print(pd.crosstab(x.v060_idiomorphs.fillna('none'), x.idiomorphs_called.fillna('none')))
print('gained/lost vs records run')
print(pd.crosstab(G.records_run_idiomorphs.fillna('none'), G.idiomorphs_called.fillna('none')))
