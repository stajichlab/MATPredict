import csv, glob, os, sys, collections, statistics as st, yaml
try: L = yaml.CSafeLoader
except AttributeError: L = yaml.SafeLoader
W = sys.argv[1]; O = f'{W}/results/2026-10-05_dothideomycetes_full'
ch = list(csv.DictReader(open(f'{O}/compare_vs_v060_genome_changes.tsv'), delimiter='\t'))
par = [r for r in ch if r['species'].startswith('Parastagonospora') and r['change'] == 'idiomorph_changed']
gen = {r['genome']: r for r in csv.DictReader(open(f'{W}/results/2026-10-03_ascomycota_v060/genomes.tsv'), delimiter='\t')}
allpar = [g for g, r in gen.items() if r['species'].startswith('Parastagonospora')]
old = collections.defaultdict(list)
for r in csv.DictReader(open(f'{W}/results/2026-10-03_ascomycota_v060/loci.tsv'), delimiter='\t'): old[r['genome']].append(r)
print('Parastagonospora genomes in the class:', len(allpar), '| idiomorph changed:', len(par))
oc = collections.Counter(old[g][0]['idiomorph'] if old[g] else 'none' for g in allpar)
nc = collections.Counter(); margins = []; genes = collections.Counter(); ids = collections.Counter()
for g in allpar:
    f = glob.glob(f'{O}/small_*/runs/{g}/detection_report.yaml')
    if not f: nc['no report'] += 1; continue
    rep = yaml.load(open(f[0]), Loader=L); det = rep.get('detected') or []
    if not det: nc['none'] += 1; continue
    x = det[0]; nc[x['idiomorph']] += 1
    c = sorted((y['score'] for y in (x.get('idiomorph_candidates') or [])), reverse=True)
    if len(c) > 1: margins.append(c[0] - c[1])
    genes[tuple(sorted(set(x['genes_found'])))] += 1
print('v0.6.0 idiomorph across all Parastagonospora genomes:', dict(oc))
print('with records:', dict(nc))
if margins: print('new-call margin (best minus second candidate score): median %.1f  min %.1f  n=%d' % (st.median(margins), min(margins), len(margins)))
print('gene sets in new calls:', dict(genes.most_common(4)))
og = collections.Counter(tuple(sorted(set(old[g][0]['genes_found'].split('|')))) for g in allpar if old[g])
print('gene sets in v0.6.0 calls:', dict(og.most_common(4)))
