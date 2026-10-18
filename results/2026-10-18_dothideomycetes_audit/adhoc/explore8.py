import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100); pd.set_option('display.max_colwidth', 100)
D = sys.argv[1]
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
L = pd.read_csv(D + '/loci_audit_flagged.tsv', sep='\t', low_memory=False)
m = G[(G.genus == 'Macrophomina') & G.primary_idiomorph.notna()]
print(pd.crosstab(m.primary_idiomorph, m.domain_final))
mm = L[(L.genus == 'Macrophomina') & L.called]
print(mm.groupby('idiomorph')[['own_id', 'own_cov', 'win_bits', 'prof_alpha_on_alpha', 'prof_hmg_on_hmg']].median())
two = G[G.idiomorphs_called.astype(str).str.contains(r'\+')]
x = L[L.genome.isin(two.genome) & L.called]
mc = x.groupby('genome').own_cov.min()
print('two-idiomorph genomes:', len(two), ' with a locus of own-gene coverage < 35%:', int((mc < 35).sum()), ' < 25%:', int((mc < 25).sum()))
print(x[x.own_cov < 35].groupby(['genus']).size().sort_values(ascending=False).head(10))
# Aureobasidium split
a = G[(G.genus == 'Aureobasidium') & G.primary_idiomorph.notna()]
print(pd.crosstab(a.species.where(a.species.isin(['Aureobasidium pullulans', 'Aureobasidium melanogenum', 'Aureobasidium subglaciale']), 'other'), a.primary_idiomorph))
print(a.gw_class.value_counts())
# Pyrenophora tritici-repentis details
p = L[(L.species == 'Pyrenophora tritici-repentis') & L.called]
print(p.confidence.value_counts(), p.locus_class.value_counts(), p.genes_found.value_counts().head(5))
# Botryosphaeria dothidea details
b = L[(L.species == 'Botryosphaeria dothidea') & L.called]
print(b[['own_id', 'own_cov', 'h_id', 'h_cov', 'h_bits', 'prof_hmg_on_hmg', 'prof_alpha_on_alpha']].median())
bg = G[(G.species == 'Botryosphaeria dothidea') & G.primary_idiomorph.notna()]
print(pd.to_numeric(bg.gw_hmg_orf_len_at_locus, errors='coerce').describe())
print(pd.to_numeric(bg.gw_hmg_best_at_locus, errors='coerce').describe())
print(G[G.gw_class.notna()].groupby('genus').gw_hmg_orf_len_at_locus.median().loc[['Aureobasidium', 'Botryosphaeria', 'Pyrenophora', 'Didymella', 'Friedmanniomyces']])
