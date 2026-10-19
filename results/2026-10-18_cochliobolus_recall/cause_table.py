"""Genome-level cause table for the 101 uncalled Cochliobolus-clade genomes. Writes cause_by_genome.tsv and prints the genus x cause counts."""
import csv
from collections import Counter, defaultdict
gates = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_gates.tsv'), delimiter='\t'))
flank = {(r['genome'], r['gene']): r for r in csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_flank.tsv'), delimiter='\t')}
cls = {r['genome']: r for r in csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_genomes.tsv'), delimiter='\t')}
def num(x):
    try: return float(x)
    except: return None
by = defaultdict(list)
for g in gates:
    if g['clade'] == 'clade': by[g['genome']].append(g)
rows = []
for genome, c in cls.items():
    if c['clade'] != 'clade': continue
    gs = by.get(genome, [])
    n50 = None
    st = c['status']
    if not gs:
        cause = 'no MAT-box ORF found by scan'
        adj = ''
        idn = ''
    else:
        n50 = num(gs[0]['n50'])
        stages = {g['stage'] for g in gs}
        ids = [num(g['cluster_identity']) for g in gs if num(g['cluster_identity']) is not None]
        idn = max(ids) if ids else ''
        adjs = []
        for g in gs:
            f = flank[(genome, g['gene'])]
            ok = lambda k: f.get(k + '_same') == 'True' and (num(f.get(k + '_id')) or 0) >= 50 and (num(f.get(k + '_dist')) or 1e9) <= 25000
            adjs.append(ok('GAP1') or ok('BGL1'))
        adj = 'GAP1/BGL1 beside MAT gene' if all(adjs) else 'MAT gene on a contig without GAP1/BGL1' if not any(adjs) else 'mixed'
        if st == 'called': cause = 'called (baseline)'
        elif 'withheld' in stages and 'evidence_floor' not in stages: cause = 'modelled_gene_bar (gene in an admitted cluster with another gene)'
        elif 'evidence_floor' in stages and 'withheld' not in stages: cause = 'evidence floor (one gene in the cluster)'
        else: cause = 'mixed (floor + bar)'
    rows.append(dict(genome=genome, genus=c['genus'], species=c['species'], strain=c['strain'], status=st, cause=cause, n_scan_genes=len(gs),
                     genes=','.join(g['gene'] for g in gs), best_cluster_identity=idn, adjacency=adj, n50=n50))
with open('results/2026-10-18_cochliobolus_recall/cause_by_genome.tsv', 'w') as f:
    w = csv.DictWriter(f, list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
t = Counter((r['genus'], r['cause']) for r in rows)
for k in sorted(t): print(k, t[k])
print()
t2 = Counter((r['genus'], r['adjacency']) for r in rows if r['status'] != 'called')
for k in sorted(t2): print('adjacency', k, t2[k])
n = Counter((r['genus'], 'N50<100kb' if (r['n50'] or 0) < 100000 else 'N50 100kb-1Mb' if (r['n50'] or 0) < 1e6 else 'N50>=1Mb') for r in rows if r['status'] != 'called')
for k in sorted(n): print('n50', k, n[k])
print('total uncalled', sum(1 for r in rows if r['status'] != 'called'))
