import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100); pd.set_option('display.max_columns', 30)
D = sys.argv[1]
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
G['order'] = G['order'].fillna('(blank)')
c = G[G.primary_idiomorph.notna()]
big = ['Pleosporales', 'Botryosphaeriales', 'Mycosphaerellales', 'Dothideales', 'Venturiales', 'Cladosporiales', 'Myriangiales', 'unclassified Dothideomycetes order', 'Trypetheliales', 'Capnodiales']
s = G[~G.order.isin(big)]
print('small orders: genomes', len(s), 'orders', s.order.nunique(), 'report', (s.status != 'no_report').sum(), 'called', s.primary_idiomorph.notna().sum(),
      'MAT1-1', (s.primary_idiomorph == 'MAT1-1').sum(), 'MAT1-2', (s.primary_idiomorph == 'MAT1-2').sum(), 'two', s.idiomorphs_called.astype(str).str.contains(r'\+').sum(),
      'notdet', (s.status == 'not_detected').sum(), 'v060', (s.v060_idiomorphs.astype(str).str.len().gt(3)).sum(),
      'conf', (s.verdict_class == 'confirmed').sum(), 'confl', (s.verdict_class == 'conflicting').sum(), 'unv', (s.verdict_class == 'unverifiable').sum())
print(pd.crosstab(c.primary_idiomorph, c.domain_final))
base = ['Alternaria', 'Parastagonospora', 'Venturia', 'Zymoseptoria', 'Ascochyta', 'Corynespora', 'Hortaea', 'Cercospora', 'Fulvia', 'Dothistroma', 'Macrophomina']
b = c[c.genus.isin(base) | (c.species == 'Pyrenophora teres')]
print('baseline calls', len(b), 'MAT1-1', (b.primary_idiomorph == 'MAT1-1').sum(), 'MAT1-2', (b.primary_idiomorph == 'MAT1-2').sum(),
      'second gene at locus (scan weak/strong):', b.second_gene_scan.isin(['strong', 'weak']).sum(), 'scanned', (b.gw_scanned.astype(str) == 'True').sum())
print(b.groupby('genus').size().to_dict())
grp = pd.read_csv(D + '/tables/group_genus.tsv', sep='\t')
print('genera with >=1 conflicting genome:', (grp.conflicting > 0).sum(), ' genomes', grp.conflicting.sum())
st = pd.read_csv(D + '/genome_status.tsv', sep='\t', low_memory=False)
print(st.two_idiomorphs.value_counts().head(3).to_dict(), st.zygosity.notna().sum())
m = c[c.genus == 'Macrophomina']
print(pd.crosstab(m.primary_idiomorph, m.gw_class))
print(m[m.primary_idiomorph == 'MAT1-2'].gw_hmg_best_at_locus.describe())
cov = pd.read_csv(D + '/biosample/biosample_dothideo_coverage.tsv', sep='\t')
print('biosample fetched', cov.fetched.sum(), 'of', len(cov))
t = G[G.genus.isin(['Friedmanniomyces', 'Elasticomyces', 'Cryomyces', 'Recurvomyces', 'Meristemomyces'])]
print('two-loci in black fungi', t.idiomorphs_called.astype(str).str.contains(r'\+').sum(), t.genus.value_counts().to_dict())
print('capitalensis profile undetermined', ((G.species == 'Phyllosticta capitalensis') & (G.primary_vs_profile == 'profile_undetermined')).sum())
