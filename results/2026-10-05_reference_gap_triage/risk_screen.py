"""Screen the v0.6.0 BFD campaigns for results that could be missed or wrong. Usage: risk_screen.py <worktree> <out_prefix>  (python3.12 with pyarrow)
1. Uncalled genomes split by assembly quality (BFD BUSCO and N50): poor assemblies explain a miss, good ones are reference gaps or true absence.
2. Idiomorph-ratio skew within species (Ascomycota, one idiomorph per genome): a 1:1 ratio is expected for a heterothallic species sampled as a population;
   strong skew flags possible mislabelling (the Parastagonospora case was 178:4 before records, 91:91 after). A screen, not a verdict: clonal lineages and homothallic species also skew.
3. Weak calls: low or medium confidence, partial locus or gene only, fallback-routed."""
import collections, csv, math, statistics as st, sys
import pyarrow.parquet as pq
W, PFX = sys.argv[1:3]
T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/tables/'
busco = {r['ASMID']: r['complete_pct'] for r in pq.read_table(T + 'busco_genome.parquet', columns=['ASMID', 'complete_pct']).to_pylist()}
asm = {r['ASMID']: r for r in pq.read_table(T + 'asm_stats.parquet', columns=['ASMID', 'contig_count', 'N50_bp', 'total_length_bp']).to_pylist()}
def rd(p): return list(csv.DictReader(open(p), delimiter='\t'))
out = open(PFX + '_summary.txt', 'w')
def P(s=''): out.write(s + '\n'); print(s)
def poor(g):
    b = busco.get(g); a = asm.get(g)
    if b is None and a is None: return None
    return (b is not None and b < 70) or (a is not None and (a['N50_bp'] < 20000 or a['contig_count'] > 5000))
res = {}
for ph, gp, lp in (('Ascomycota', f'{W}/results/2026-10-03_ascomycota_v060/genomes.tsv', f'{W}/results/2026-10-03_ascomycota_v060/loci.tsv'),
                   ('Basidiomycota', f'{W}/results/2026-10-03_basidiomycota_v060/genomes.tsv', f'{W}/results/2026-10-03_basidiomycota_v060/loci.tsv')):
    G = rd(gp); L = rd(lp); lbg = collections.defaultdict(list)
    for r in L: lbg[r['genome']].append(r)
    P(f'\n## {ph}: {len(G)} genomes')
    # 1. misses by assembly quality
    for status in ('called',):
        cats = collections.Counter()
        for r in G:
            q = poor(r['genome']); cats[(('poor' if q else 'good' if q is False else 'no stats'), r['status'] == 'called')] += 1
        for q in ('good', 'poor', 'no stats'):
            n = cats[(q, True)] + cats[(q, False)]
            if n: P(f'assembly {q}: {n} genomes, called {cats[(q, True)]} ({100*cats[(q, True)]/n:.1f}%), not called {cats[(q, False)]}')
    P('\nNot called, by order (orders with 20 or more genomes): genomes | not called | not called with a good assembly | % of good-assembly genomes not called | routing')
    by = collections.defaultdict(list)
    for r in G: by[r['order'] or 'order not assigned'].append(r)
    rows = []
    for o, v in by.items():
        if len(v) < 20: continue
        good = [r for r in v if poor(r['genome']) is False]; ng = [r for r in good if r['status'] != 'called']
        nc = [r for r in v if r['status'] != 'called']
        fb = 100 * sum(r['routing'] == 'phylum_fallback' for r in v) / len(v)
        rows.append((o, len(v), len(nc), len(ng), 100 * len(ng) / max(1, len(good)), fb))
    P('| order | genomes | not called | not called, good assembly | % good-assembly not called | % fallback |'); P('|---|---|---|---|---|---|')
    for r in sorted(rows, key=lambda r: -r[3])[:16]: P(f'| {r[0]} | {r[1]} | {r[2]} | {r[3]} | {r[4]:.1f}% | {r[5]:.0f}% |')
    allgood = [r for r in G if poor(r['genome']) is False]
    P(f'\nTotal not called in good assemblies: {sum(r["status"] != "called" for r in allgood)} of {len(allgood)} good-assembly genomes')
    # baseline: orders with a curated record and lineage routing and >=90% called
    res[ph] = (G, lbg)
    # 3. weak calls
    conf = collections.Counter(); lc = collections.Counter(); gweak = set(); gall = set()
    for r in L:
        conf[r['confidence']] += 1; lc[r['locus_class']] += 1
    gen_fb = {r['genome'] for r in G if r['routing'] == 'phylum_fallback' and r['status'] == 'called'}
    for g, ls in lbg.items():
        if all(x['confidence'] == 'low' or x['locus_class'] in ('idiomorph_gene_only', 'partial_locus') for x in ls): gweak.add(g)
    P(f'\nLoci: confidence {dict(conf)}; locus class {dict(lc)}')
    P(f'Genomes whose every call is low confidence, partial or gene only: {len(gweak)}; called through the phylum fallback: {len(gen_fb)}')
# 2. skew screen, Ascomycota
G, lbg = res['Ascomycota']; sp = collections.defaultdict(lambda: [0, 0, 0])
for r in G:
    ls = lbg.get(r['genome'], [])
    ids = {x['idiomorph'] for x in ls if x['family_called'] == 'Ascomycota:MAT'}
    if not ids: continue
    k = r['species']
    if ids == {'MAT1-1'}: sp[k][0] += 1
    elif ids == {'MAT1-2'}: sp[k][1] += 1
    else: sp[k][2] += 1
order_of = {r['species']: r['order'] for r in G}
def binom_p(a, n):  # two-sided tail p for min count a out of n at 0.5
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(0, a + 1)) / 2 ** n)
sk = []
for s, (a, b, c) in sp.items():
    n = a + b
    if n < 10: continue
    mn = min(a, b)
    if binom_p(mn, n) < 1e-4 and mn / n <= 0.15: sk.append((s, order_of.get(s, ''), a, b, c, binom_p(mn, n)))
P(f'\n## Idiomorph skew within species (Ascomycota:MAT, species with 10 or more single-idiomorph genomes; minority share 15% or less, two-sided binomial p < 1e-4 against 1:1)')
P(f'species screened: {sum(1 for s, v in sp.items() if v[0] + v[1] >= 10)}; flagged: {len(sk)}; genomes in flagged species: {sum(a + b for _, _, a, b, _, _ in sk)}')
P('| species | order | MAT1-1 | MAT1-2 | both |'); P('|---|---|---|---|---|')
for s in sorted(sk, key=lambda x: -(x[2] + x[3]))[:25]: P(f'| {s[0]} | {s[1]} | {s[2]} | {s[3]} | {s[4]} |')
