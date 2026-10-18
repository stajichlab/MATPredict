import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 200); pd.set_option('display.max_colwidth', 60)
D = sys.argv[1]
L = pd.read_csv(D + '/loci_audit_flagged.tsv', sep='\t', low_memory=False)
c = L[L.called]
x = c[c.core_lost_to_optional.astype(str) == 'True']
print('core lost to optional: loci', len(x), 'genomes', x.genome.nunique())
print(x.groupby(['genus', 'species', 'idiomorph']).size())
print(x.res_pairs.value_counts().head(8))
y = c[c.win_optional.astype(str) == 'True']
print('winning gene optional:'); print(y[['genome', 'species', 'idiomorph', 'confidence', 'locus_class', 'win_gene', 'win_bits', 'genes_found', 'own_present']].to_string())
z = c[c.locus_class == 'idiomorph_gene_only']
print('gene_only by genus'); print(z.groupby(['genus', 'idiomorph']).size().sort_values(ascending=False).head(12))
print(z[['own_id', 'own_cov']].describe().T)
print('partial_locus domain/profile'); print(pd.crosstab(c.locus_class, [c.vote_vs_profile]))
print(pd.crosstab(c.locus_class, [c.domain_support]))
print('own_id<40:'); w = c[c.own_id < 40]; print(w.groupby(['genus', 'idiomorph']).size().sort_values(ascending=False).head(10))
print('confidence x idiomorph'); print(pd.crosstab(c.idiomorph, c.confidence))
