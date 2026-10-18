"""Species with a skewed idiomorph ratio: is the skew carried by one BioProject (clonal sampling) or spread across several?
python explore6.py DIR -> tables/species_skew_by_bioproject.tsv"""
import sys
from pathlib import Path
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 200)
D = Path(sys.argv[1])
G = pd.read_csv(D / 'genome_audit_final.tsv', sep='\t', low_memory=False)
B = pd.read_csv(D / 'bioproject.tsv', sep='\t')
G = G.merge(B[['genome', 'bioproject']], on='genome', how='left')
sp = pd.read_csv(D / 'tables' / 'species_ratio.tsv', sep='\t')
skew = sp[(sp.n >= 8) & (sp.p_binom_1to1 < 0.05)]
rows = []
for s in skew.species:
    x = G[(G.species == s) & G.primary_idiomorph.notna()]
    t = x.groupby('bioproject').primary_idiomorph.agg(n='size', m1=lambda v: (v == 'MAT1-1').sum(), m2=lambda v: (v == 'MAT1-2').sum()).sort_values('n', ascending=False)
    for bp, r in t.head(6).iterrows():
        rows.append(dict(species=s, bioproject=bp, n=int(r.n), MAT11=int(r.m1), MAT12=int(r.m2), species_n=len(x), n_bioprojects=len(t)))
out = pd.DataFrame(rows)
out.to_csv(D / 'tables' / 'species_skew_by_bioproject.tsv', sep='\t', index=False)
print(skew.to_string())
print(out.to_string())
