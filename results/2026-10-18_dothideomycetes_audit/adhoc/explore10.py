import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100)
D = sys.argv[1]
L = pd.read_csv(D + '/loci_audit_flagged.tsv', sep='\t', low_memory=False)
c = L[L.called]
for k in ('F_optional_wins', 'F_own_core_missing', 'F_own_partial_edge', 'F_own_partial_interior', 'F_low_identity', 'F_profile_disagree'):
    print(k, int(c[k].astype(str).eq('True').sum()))
print(c[c.F_own_partial_edge.astype(str) == 'True'].groupby('genus').size().sort_values(ascending=False).head(8))
print('own_cov quantiles by order')
print(c.groupby('order').own_cov.describe()[['count', '25%', '50%', '75%']].round(0).sort_values('count', ascending=False).head(8))
print('own_full by order'); print(c.groupby('order').own_full.apply(lambda s: s.astype(str).eq('True').sum()).sort_values(ascending=False).head(6))
print(c.genes_found.str.contains('SLA2').sum(), 'with SLA2;', c.genes_found.str.contains('APN2').sum(), 'APN2;', c.genes_found.str.contains('COX13').sum(), 'COX13')
print('loci with any flank:', (c.flank_found.notna()).sum())
print(c.n_loci_genome.value_counts().head())
