import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_columns', 40); pd.set_option('display.max_rows', 200)
D = sys.argv[1]
G = pd.read_csv(D + '/genome_audit_verdict.tsv', sep='\t', low_memory=False)
G = G[G.gw_scanned.astype(str) == 'True']
G = G[G.primary_idiomorph.notna()]
G['at_hmg'] = pd.to_numeric(G.gw_hmg_best_at_locus, errors='coerce')
G['at_alpha'] = pd.to_numeric(G.gw_alpha_best_at_locus, errors='coerce')
for gen in sys.argv[2:]:
    x = G[G.genus == gen]
    for idio, col in (('MAT1-1', 'at_hmg'), ('MAT1-2', 'at_alpha')):
        y = x[x.primary_idiomorph == idio]
        if len(y):
            print(gen, idio, 'n', len(y), col, 'present', y[col].notna().sum(), 'quantiles',
                  y[col].quantile([.1, .5, .9]).round(0).tolist() if y[col].notna().any() else None)
