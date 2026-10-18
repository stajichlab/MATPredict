"""Paralog HMG-box check for MAT1-2 calls: is the HMG ORF at the called locus the best-scoring HMG ORF in the genome (MAT1-2-1 profile)?
python paralog_check.py DIR -> tables/hmg_paralog_check.tsv"""
import sys
from pathlib import Path
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100)
D = Path(sys.argv[1])
G = pd.read_csv(D / 'genome_audit_final.tsv', sep='\t', low_memory=False)
x = G[(G.primary_idiomorph == 'MAT1-2') & (G.gw_scanned.astype(str) == 'True')].copy()
x['hmg_locus'] = pd.to_numeric(x.gw_hmg_best_at_locus, errors='coerce')
x['hmg_best'] = pd.to_numeric(x.gw_hmg_best, errors='coerce')
x['gap'] = x['hmg_best'] - x['hmg_locus']
x['no_hmg_orf_at_locus'] = x['hmg_locus'].isna()
x['elsewhere_higher_by_10'] = x['gap'] >= 10
x['elsewhere_higher_by_30'] = x['gap'] >= 30
t = x.groupby(['order', 'genus']).agg(n=('genome', 'size'), no_hmg_orf_at_locus=('no_hmg_orf_at_locus', 'sum'),
                                      elsewhere_higher_by_10=('elsewhere_higher_by_10', 'sum'),
                                      elsewhere_higher_by_30=('elsewhere_higher_by_30', 'sum'),
                                      med_locus=('hmg_locus', 'median'), med_best=('hmg_best', 'median')).reset_index()
t.sort_values('n', ascending=False).to_csv(D / 'tables' / 'hmg_paralog_check.tsv', sep='\t', index=False)
print('MAT1-2 primary calls scanned:', len(x), ' no HMG ORF at locus:', int(x.no_hmg_orf_at_locus.sum()),
      ' off-locus HMG ORF higher by >=10 bits:', int(x.elsewhere_higher_by_10.sum()), ' by >=30:', int(x.elsewhere_higher_by_30.sum()))
print(t.sort_values('elsewhere_higher_by_10', ascending=False).head(15).to_string())
