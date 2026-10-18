import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100); pd.set_option('display.max_colwidth', 100)
D = sys.argv[1]
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
print(G.verdict.value_counts())
c = G[G.primary_idiomorph.notna()]
print('second gene pfam x scan'); print(pd.crosstab(c.second_gene_pfam, c.second_gene_scan))
print('domain_final'); print(c.domain_final.value_counts())
print('primary conf x verdict'); print(pd.crosstab(c.primary_confidence, c.verdict_class))
print('class x verdict'); print(pd.crosstab(c.primary_class, c.verdict_class))
print('scan baseline: genera with >=20 called genomes and zero strong/weak second gene')
t = c.groupby('genus').agg(n=('genome', 'size'), sg=('second_gene_scan', lambda s: (s.isin(['strong', 'weak'])).sum()))
print(t[(t.n >= 20)].sort_values('n', ascending=False))
print('both idiomorph loci genomes', c.idiomorphs_called.astype(str).str.contains(r'\+').sum())
z = G[G.genome.isin(['GCF_000219625.1_MYCGR_v2.0'])]
print(G[G.status == 'no_report'][['genome', 'species']])
M = pd.read_csv(D + '/genome_list.tsv', sep='\t', dtype=str)
print(M[(M.class_by_csv == 'True') & (M.taxid_under_147541 == 'False')][['ASM_ACCESSION', 'NCBI_TAXID', 'SPECIES_IN']])
print(M[M.has_fasta == 'False'][['ASM_ACCESSION', 'SPECIES_IN']])
