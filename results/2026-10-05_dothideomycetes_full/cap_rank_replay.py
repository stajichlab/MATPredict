"""Replay the polish-cap ranking on the Dothideomycetes full run and test alternative culling rules.

Truth = loci called in v0.6.0, in the full run with records, or in the cap-off re-run. For each genome the admitted candidate clusters
(Ascomycota:MAT evidence rows with admitted=true) are ranked under several rules; a true locus is KEPT if its cluster is in the top `cap`.
Usage: cap_rank_replay.py <worktree> <out_prefix>
"""
import collections, csv, glob, json, os, statistics as st, sys
import yaml
try: L = yaml.CSafeLoader
except AttributeError: L = yaml.SafeLoader
W, PFX = sys.argv[1:3]
O = f'{W}/results/2026-10-05_dothideomycetes_full'; T = f'{W}/results/2026-10-05_dothideo_cap0_lost'; OLD = f'{W}/results/2026-10-03_ascomycota_v060'
genomes = [l.split('\t')[0] for l in open(f'{O}/list.tsv')]
truth = collections.defaultdict(list)
for r in csv.DictReader(open(f'{OLD}/loci.tsv'), delimiter='\t'):
    truth[r['genome']].append((r['contig'], int(r['start']), int(r['end']), 'v060'))
def add_reports(pattern, tag):
    for p in glob.glob(pattern):
        g = os.path.basename(os.path.dirname(p)); rep = yaml.load(open(p), Loader=L)
        for x in rep.get('detected') or []: truth[g].append((x['contig'], x['start'], x['end'], tag))
add_reports(f'{O}/small_*/runs/*/detection_report.yaml', 'new'); add_reports(f'{T}/out/runs/*/detection_report.yaml', 'cap0')
lost = {l.strip() for l in open(f'{T}/lost.txt')}
SCHEMES = {
    'current: genes, identity, hits': lambda r: (-r['gene_count'], -r['best_identity'], -r['hit_count']),
    'identity first': lambda r: (-r['best_identity'], -r['gene_count'], -r['hit_count']),
    'strong hit (>=60%) first, then current': lambda r: (-(r['best_identity'] >= 60), -r['gene_count'], -r['best_identity'], -r['hit_count']),
    'strong hit (>=50%) first, then current': lambda r: (-(r['best_identity'] >= 50), -r['gene_count'], -r['best_identity'], -r['hit_count']),
}
CAPS = [6, 10, 15, 20, 30, 10**6]
res = {s: {c: 0 for c in CAPS} for s in SCHEMES}; ntruth = 0; notadm = 0
nadm = []; true_ident = []; noise_ident = []; cost = {c: 0 for c in CAPS}; validate = [0, 0]; lostrank = []
for g in genomes:
    f = glob.glob(f'{O}/small_*/runs/{g}/evidence_diagnostics.jsonl')
    if not f: continue
    rows = [x for x in (json.loads(l) for l in open(f[0])) if x['kind'] == 'evidence' and x['family'] == 'Ascomycota:MAT' and x['admitted']]
    nadm.append(len(rows))
    for c in CAPS: cost[c] += min(c, len(rows))
    # validate: replay of the current rank should put exactly the polish_capped=false rows in the top 6
    srt = sorted(rows, key=lambda r: (SCHEMES['current: genes, identity, hits'](r), r['cluster_start']))
    top6 = {id(r) for r in srt[:6]}; unc = {id(r) for r in rows if not r['polish_capped']}
    validate[0] += len(top6 & unc); validate[1] += max(len(top6), len(unc))
    # map truth loci to clusters (largest overlap on the same contig)
    seen = set()
    for contig, a, b, tag in truth.get(g, []):
        best = None
        for r in rows:
            if r['contig'] != contig: continue
            ov = min(b, r['cluster_end']) - max(a, r['cluster_start'])
            if ov > 0 and (best is None or ov > best[0]): best = (ov, r)
        if best is None: notadm += 1; continue
        r = best[1]
        if id(r) in seen: continue
        seen.add(id(r)); ntruth += 1; true_ident.append(r['best_identity'])
        for s, key in SCHEMES.items():
            order = sorted(rows, key=lambda x: (key(x), x['cluster_start'])); rank = next(i for i, x in enumerate(order) if x is r) + 1
            for c in CAPS:
                if rank > c: res[s][c] += 1
            if s.startswith('current') and g in lost: lostrank.append(rank)
    truers = {id(r) for r in rows for contig, a, b, tag in truth.get(g, []) if r['contig'] == contig and min(b, r['cluster_end']) - max(a, r['cluster_start']) > 0}
    noise_ident += [r['best_identity'] for r in rows if id(r) not in truers]
out = open(PFX + '_summary.txt', 'w')
def P(s=''): out.write(s + '\n'); print(s)
P(f'genomes with diagnostics {len(nadm)}; true loci mapped to an admitted cluster {ntruth}; true loci with no admitted cluster {notadm}')
P(f'admitted clusters per genome: median {st.median(nadm)}, p90 {sorted(nadm)[int(.9*len(nadm))]}, max {max(nadm)}; genomes with more than 6: {sum(n > 6 for n in nadm)} ({100*sum(n > 6 for n in nadm)/len(nadm):.0f}%)')
P(f'replay check: top-6 of the current rank agrees with the clusters actually polished in {100*validate[0]/validate[1]:.1f}% of cluster slots')
P(f'rank of the true cluster in the 15 lost genomes under the current rank: {sorted(lostrank)}')
P('')
P('Calls (true loci) lost, by ranking rule and cap:'); P('| rule | ' + ' | '.join('cap ' + (str(c) if c < 10**6 else 'none') for c in CAPS) + ' |'); P('|---|' + '---|' * len(CAPS))
for s in SCHEMES: P(f'| {s} | ' + ' | '.join(str(res[s][c]) for c in CAPS) + ' |')
P(''); P('Polishing work (clusters polished, relative to cap 6): ' + ', '.join(f'cap {c if c < 10**6 else "none"}: {cost[c]} ({cost[c]/cost[6]:.2f}x)' for c in CAPS))
P(''); P('Identity floor test (best identity of the cluster): true clusters below t lost vs noise clusters removed')
ti = sorted(true_ident); ni = sorted(noise_ident)
P(f'true clusters {len(ti)}: min {ti[0]:.1f}, 5th pct {ti[int(.05*len(ti))]:.1f}, median {st.median(ti):.1f}; noise clusters {len(ni)}: median {st.median(ni):.1f}, 95th pct {ni[int(.95*len(ni))]:.1f}')
P('| floor t | true clusters below t | noise clusters below t |'); P('|---|---|---|')
for t in (30, 35, 40, 45, 50, 60):
    P(f'| {t}% | {sum(x < t for x in ti)} ({100*sum(x < t for x in ti)/len(ti):.1f}%) | {sum(x < t for x in ni)} ({100*sum(x < t for x in ni)/len(ni):.1f}%) |')
