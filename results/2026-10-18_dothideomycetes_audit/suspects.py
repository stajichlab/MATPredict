"""Evidence table for the named suspect and reference groups.  python suspects.py DIR -> tables/suspect_groups.tsv"""
import sys
from pathlib import Path
import pandas as pd
pd.set_option('display.width', 300); pd.set_option('display.max_columns', 60); pd.set_option('display.max_colwidth', 60)
D = Path(sys.argv[1])
G = pd.read_csv(D / 'genome_audit_final.tsv', sep='\t', low_memory=False)
L = pd.read_csv(D / 'loci_audit_flagged.tsv', sep='\t', low_memory=False)
c = L[L.called]
prim = c.assign(rc=c.confidence.map({'high': 0, 'medium': 1, 'low': 2}).fillna(3),
                rl=c.locus_class.map({'mat_locus': 0, 'partial_locus': 1, 'idiomorph_gene_only': 2}).fillna(3)) \
    .sort_values(['rc', 'rl', 'win_bits'], ascending=[True, True, False]).drop_duplicates('genome').set_index('genome')
for k in ('own_id', 'own_cov', 'win_gene', 'locus_class', 'confidence'):
    G['p_' + k] = G.genome.map(prim[k])
G['hmg_at'] = pd.to_numeric(G.get('gw_hmg_best_at_locus'), errors='coerce')
G['alpha_at'] = pd.to_numeric(G.get('gw_alpha_best_at_locus'), errors='coerce')

groups = [
    ('Aureobasidium (all)', lambda d: d.genus == 'Aureobasidium', 'homothallic (Gostincar 2014, PMID 24984952)'),
    ('Pyrenophora tritici-repentis', lambda d: d.species == 'Pyrenophora tritici-repentis', 'homothallic (PMIDs 20373969, 29050833)'),
    ('Pyrenophora teres', lambda d: d.species == 'Pyrenophora teres', 'heterothallic (Rau 2005, not re-opened)'),
    ('Botryosphaeria dothidea', lambda d: d.species == 'Botryosphaeria dothidea', 'not established here (Nagel 2018, PMID 29530630, abstract only)'),
    ('Macrophomina phaseolina', lambda d: d.species == 'Macrophomina phaseolina', 'not established here (Nagel 2018, abstract only)'),
    ('Macrophomina other', lambda d: (d.genus == 'Macrophomina') & (d.species != 'Macrophomina phaseolina'), ''),
    ('Friedmanniomyces', lambda d: d.genus == 'Friedmanniomyces', 'none found'),
    ('Elasticomyces', lambda d: d.genus == 'Elasticomyces', 'none found'),
    ('Cryomyces/Recurvomyces/Meristemomyces', lambda d: d.genus.isin(['Cryomyces', 'Recurvomyces', 'Meristemomyces']), 'none found'),
    ('Didymella pinodes', lambda d: d.species == 'Didymella pinodes', 'homothallic (Woudenberg 2012, PMID 22014305; as Peyronellaea pinodes)'),
    ('Didymella pinodella', lambda d: d.species == 'Didymella pinodella', 'heterothallic (PMID 22014305; as Pey. pinodella)'),
    ('Didymella other', lambda d: (d.genus == 'Didymella') & ~d.species.isin(['Didymella pinodes', 'Didymella pinodella']), 'mixed (PMID 22014305)'),
    ('Stemphylium', lambda d: d.genus == 'Stemphylium', 'homothallic species known (from memory / repo note, Inderbitzin 2005; not re-opened)'),
    ('Phyllosticta citricarpa', lambda d: d.species == 'Phyllosticta citricarpa', 'heterothallic (repo note: KT708823/KX280782)'),
    ('Phyllosticta capitalensis', lambda d: d.species == 'Phyllosticta capitalensis', 'homothallic (repo note KT708824; unpublished deposit)'),
    ('Cercospora', lambda d: d.genus == 'Cercospora', 'heterothallic + off-locus fragments (Bolton 2014, PMID 24216224)'),
    ('Cladosporium', lambda d: d.genus == 'Cladosporium', 'none found'),
    ('Alternaria alternata', lambda d: d.species == 'Alternaria alternata', 'asexual; one idiomorph per isolate (PMIDs 26666175, 23441968)'),
    ('Alternaria other', lambda d: (d.genus == 'Alternaria') & (d.species != 'Alternaria alternata'), ''),
    ('Parastagonospora', lambda d: d.genus == 'Parastagonospora', 'heterothallic (Bennett 2003, PMID 12948511)'),
    ('Zymoseptoria', lambda d: d.genus == 'Zymoseptoria', 'heterothallic (Waalwijk 2002, PMID 11929216)'),
    ('Venturia', lambda d: d.genus == 'Venturia', 'heterothallic (PMID 19689797; V. effusa PMID 29381450)'),
    ('Corynespora', lambda d: d.genus == 'Corynespora', 'none found'),
    ('Ascochyta', lambda d: d.genus == 'Ascochyta', 'heterothallic (rabiei BioSample; PMID 22014305 for Didymellaceae)'),
    ('Plenodomus', lambda d: d.genus == 'Plenodomus', 'heterothallic (Cozijnsen 2003, PMID 12679880)'),
    ('Diplodia', lambda d: d.genus == 'Diplodia', 'D. sapinea heterothallic (Bihon 2014, PMID 24220137)'),
    ('Fulvia', lambda d: d.genus == 'Fulvia', 'MAT1-1 or MAT1-2 per strain (Stergiopoulos 2007, PMID 17178244)'),
    ('Hortaea', lambda d: d.genus == 'Hortaea', 'none found'),
    ('Bipolaris', lambda d: d.genus == 'Bipolaris', 'heterothallic (Turgeon 1993, PMID 8479433)'),
    ('Curvularia/Exserohilum/Ophiosphaerella', lambda d: d.genus.isin(['Curvularia', 'Exserohilum', 'Ophiosphaerella']), ''),
    ('Neofusicoccum', lambda d: d.genus == 'Neofusicoccum', ''),
    ('Sphaerulina', lambda d: d.genus == 'Sphaerulina', ''),
]
rows = []
for name, f, lit in groups:
    d = G[f(G)]
    called = d[d.primary_idiomorph.notna()]
    sc = called[called.gw_scanned.astype(str) == 'True'] if 'gw_scanned' in called else called.iloc[0:0]
    rows.append(dict(group=name, genomes=len(d), called=len(called), MAT11=int((called.primary_idiomorph == 'MAT1-1').sum()),
                     MAT12=int((called.primary_idiomorph == 'MAT1-2').sum()),
                     vote_profile_agree=int((called.primary_vs_profile == 'agree').sum()),
                     vote_profile_disagree=int((called.primary_vs_profile == 'disagree').sum()),
                     profile_undetermined=int((called.primary_vs_profile == 'profile_undetermined').sum()),
                     pfam_both_genes_strong=int((called.second_gene_pfam == 'strong').sum()),
                     scanned=len(sc),
                     scan_second_gene_strong=int((sc.second_gene_scan == 'strong').sum()),
                     scan_second_gene_weak=int((sc.second_gene_scan == 'weak').sum()),
                     scan_none=int((sc.second_gene_scan == 'none').sum()),
                     core_set_aside_for_optional=int((called.primary_core_lost_to_optional.astype(str) == 'True').sum()),
                     med_own_id=round(called.p_own_id.median(), 1) if len(called) else None,
                     med_own_cov=round(called.p_own_cov.median(), 1) if len(called) else None,
                     confirmed=int((d.verdict_class == 'confirmed').sum()), conflicting=int((d.verdict_class == 'conflicting').sum()),
                     unverifiable=int((d.verdict_class == 'unverifiable').sum()), no_call=int((d.verdict_class == 'no_call').sum()),
                     literature=lit))
T = pd.DataFrame(rows)
T.to_csv(D / 'tables' / 'suspect_groups.tsv', sep='\t', index=False)
print(T.drop(columns=['literature']).to_string())
