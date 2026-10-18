import sys
import pandas as pd
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 100)
D = sys.argv[1]
g = pd.read_csv(D + '/tables/group_genus.tsv', sep='\t')
print(g.groupby('label').agg(genera=('genus', 'size'), genomes=('genomes', 'sum'), called=('called', 'sum'), confirmed=('confirmed', 'sum'), conflicting=('conflicting', 'sum'), unverifiable=('unverifiable', 'sum')))
print(len(g), 'genera')
f = pd.read_csv(D + '/tables/group_family.tsv', sep='\t')
print(f.groupby('label').agg(families=('family', 'size'), genomes=('genomes', 'sum')))
print(g[g.label.str.startswith('mixed')][['order', 'genus', 'genomes', 'called', 'confirmed', 'conflicting', 'unverifiable']].to_string())
print(g[g.label.str.startswith('unverifiable')][['order', 'genus', 'genomes', 'called', 'confirmed', 'conflicting', 'unverifiable']].to_string())
G = pd.read_csv(D + '/genome_audit_final.tsv', sep='\t', low_memory=False)
print(G.genus.nunique(), 'genera;', G.order.fillna('(blank)').nunique(), 'orders;', G.species.nunique(), 'species')
