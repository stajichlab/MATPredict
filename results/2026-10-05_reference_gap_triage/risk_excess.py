"""Excess misses (not called on a good assembly, beyond a 3% background) and the families holding them. Usage: risk_excess.py <worktree> <out_prefix> (python3.12 with pyarrow)"""
import collections, csv, sys
import pyarrow.parquet as pq
W, PFX = sys.argv[1:3]
T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/tables/'
busco = {r['ASMID']: r['complete_pct'] for r in pq.read_table(T + 'busco_genome.parquet', columns=['ASMID', 'complete_pct']).to_pylist()}
asm = {r['ASMID']: r for r in pq.read_table(T + 'asm_stats.parquet', columns=['ASMID', 'contig_count', 'N50_bp']).to_pylist()}
def poor(g):
    b = busco.get(g); a = asm.get(g)
    if b is None and a is None: return True
    return (b is not None and b < 70) or (a is not None and (a['N50_bp'] < 20000 or a['contig_count'] > 5000))
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
out = open(PFX + '_summary.txt', 'w')
def P(s=''): out.write(s + '\n'); print(s)
BG = 0.03
for ph, gp in (('Ascomycota', f'{W}/results/2026-10-03_ascomycota_v060/genomes.tsv'), ('Basidiomycota', f'{W}/results/2026-10-03_basidiomycota_v060/genomes.tsv')):
    G = rd(gp); by = collections.defaultdict(list)
    for r in G: by[r['order'] or 'order not assigned'].append(r)
    exc = {}; tot = 0
    for o, v in by.items():
        good = [r for r in v if not poor(r['genome'])]; ng = sum(r['status'] != 'called' for r in good)
        e = max(0.0, ng - BG * len(good)); exc[o] = (e, ng, len(good), len(v)); tot += e
    P(f'\n## {ph}: excess misses (not called on a good assembly, beyond a {int(BG*100)}% background) = {tot:.0f} genomes ({100*tot/len(G):.1f}% of {len(G)})')
    big = sorted(exc.items(), key=lambda kv: -kv[1][0])
    cum = 0
    P('| order | excess misses | not called, good assembly | good assemblies | genomes | cumulative share of all excess |'); P('|---|---|---|---|---|---|')
    for o, (e, ng, gd, n) in big[:14]:
        cum += e; P(f'| {o} | {e:.0f} | {ng} | {gd} | {n} | {100*cum/tot:.0f}% |')
    P(f'orders with any excess: {sum(1 for e in exc.values() if e[0] >= 1)}; orders carrying 80% of the excess: ' + str(next(i + 1 for i in range(len(big)) if sum(x[1][0] for x in big[:i + 1]) >= 0.8 * tot)))
    if ph == 'Ascomycota':
        for o in ('Saccharomycetales', 'Pleosporales', 'Pichiales', 'Dipodascales', 'Chaetothyriales', 'Serinales'):
            fam = collections.Counter(); famn = collections.Counter()
            for r in by[o]:
                famn[r['family'] or 'family not assigned'] += 1
                if r['status'] != 'called' and not poor(r['genome']): fam[r['family'] or 'family not assigned'] += 1
            P(f'\n{o} not called on a good assembly, by family: ' + ', '.join(f'{f} {n} of {famn[f]}' for f, n in fam.most_common(8)))
