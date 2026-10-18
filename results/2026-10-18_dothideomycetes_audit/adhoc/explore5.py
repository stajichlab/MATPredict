import sys
import pandas as pd
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 60); pd.set_option('display.max_rows', 200); pd.set_option('display.max_colwidth', 40)
D = sys.argv[1]
L = pd.read_csv(D + '/loci_audit_flagged.tsv', sep='\t', low_memory=False)
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
two = G[G.idiomorphs_called.astype(str).str.contains(r'\+')]
print(two.groupby(['order', 'genus']).size())
x = L[L.genome.isin(two.genome) & L.called & ~L.genus.isin(['Friedmanniomyces', 'Elasticomyces'])]
print(x[['genome', 'species', 'contig', 'start', 'idiomorph', 'confidence', 'locus_class', 'genes_found', 'own_id', 'own_cov', 'prof_call', 'domain_support', 'n_loci_genome']].sort_values(['genome', 'start']).to_string())
