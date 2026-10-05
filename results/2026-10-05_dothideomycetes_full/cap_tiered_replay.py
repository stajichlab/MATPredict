"""Cost and safety of a tiered polish rule: polish every admitted cluster with best identity >= T, then top up to 6 from the rest in the current rank.
Truth = loci called in v0.6.0, the full run with records, or the cap-off re-run (same definition as cap_rank_replay.py).
Usage: cap_tiered_replay.py <worktree> <out_prefix>
"""
import collections, csv, glob, json, os, statistics as st, sys
import yaml
try: L = yaml.CSafeLoader
except AttributeError: L = yaml.SafeLoader
W, PFX = sys.argv[1:3]
O = f'{W}/results/2026-10-05_dothideomycetes_full'; T = f'{W}/results/2026-10-05_dothideo_cap0_lost'; OLD = f'{W}/results/2026-10-03_ascomycota_v060'
genomes = [l.split('\t')[0] for l in open(f'{O}/list.tsv')]
meta = {r['genome']: r for r in csv.DictReader(open(f'{OLD}/genomes.tsv'), delimiter='\t')}
truth = collections.defaultdict(list)
for r in csv.DictReader(open(f'{OLD}/loci.tsv'), delimiter='\t'): truth[r['genome']].append((r['contig'], int(r['start']), int(r['end'])))
def add(pattern):
    for p in glob.glob(pattern):
        g = os.path.basename(os.path.dirname(p)); rep = yaml.load(open(p), Loader=L)
        for x in rep.get('detected') or []: truth[g].append((x['contig'], x['start'], x['end']))
add(f'{O}/small_*/runs/*/detection_report.yaml'); add(f'{T}/out/runs/*/detection_report.yaml')
cur = lambda r: (-r['gene_count'], -r['best_identity'], -r['hit_count'], r['cluster_start'])
THR = (45, 50, 60); cost = {t: 0 for t in THR}; cost6 = 0; costall = 0; lost = {t: [] for t in THR}; nstrong = {t: [] for t in THR}; weak_true = []
for g in genomes:
    f = glob.glob(f'{O}/small_*/runs/{g}/evidence_diagnostics.jsonl')
    if not f: continue
    rows = [x for x in (json.loads(l) for l in open(f[0])) if x['kind'] == 'evidence' and x['family'] == 'Ascomycota:MAT' and x['admitted']]
    cost6 += min(6, len(rows)); costall += len(rows)
    order = sorted(rows, key=cur)
    for t in THR:
        strong = [r for r in rows if r['best_identity'] >= t]; nstrong[t].append(len(strong))
        chosen = {id(r) for r in strong} | {id(r) for r in order[:6]}
        # top-up: strong plus the best of the rest until 6 polished in total
        rest = [r for r in order if r['best_identity'] < t]; chosen = {id(r) for r in strong} | {id(r) for r in rest[:max(0, 6 - len(strong))]}
        cost[t] += len(chosen)
        for contig, a, b in truth.get(g, []):
            best = None
            for r in rows:
                if r['contig'] != contig: continue
                ov = min(b, r['cluster_end']) - max(a, r['cluster_start'])
                if ov > 0 and (best is None or ov > best[0]): best = (ov, r)
            if best and id(best[1]) not in chosen: lost[t].append((g, best[1]['best_identity']))
            if best and t == 50 and best[1]['best_identity'] < 50:
                weak_true.append((g, meta[g]['species'][:30], meta[g]['order'], best[1]['best_identity'], best[1]['gene_count'], order.index(best[1]) + 1, len(strong), len(rows)))
out = open(PFX + '_summary.txt', 'w')
def P(s=''): out.write(s + '\n'); print(s)
P('Tiered rule: polish every admitted cluster with best identity >= T, then top up to 6 polished in total from the rest (current rank)')
P(f'cap 6 today: {cost6} clusters polished; no cap: {costall}')
P('| T | true loci lost | clusters polished | vs cap 6 | genomes with more than 6 strong clusters | max strong clusters in a genome |'); P('|---|---|---|---|---|---|')
for t in THR: P(f'| {t}% | {len(lost[t])} | {cost[t]} | {cost[t]/cost6:.2f}x | {sum(n > 6 for n in nstrong[t])} | {max(nstrong[t])} |')
P(''); P('Distribution of the number of admitted clusters with best identity >= 50% per genome: ' + str(dict(sorted(collections.Counter(min(n, 10) for n in nstrong[50]).items()))) + ' (10 = 10 or more)')
P(''); P('True loci whose cluster has best identity < 50% (the divergent ones a floor would cut): genome | species | order | identity | genes | rank under current rule | strong clusters in genome | admitted')
for w in sorted(weak_true, key=lambda x: x[3]): P('  ' + ' | '.join(str(x) for x in w))
