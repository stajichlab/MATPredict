"""Compare the full Dothideomycetes run on main with the Dothideomycete records (run-84361bc) with the v0.6.0 campaign (no records).

Usage: compare_v060.py <new_out_dir> <ascomycota_v060_dir> <out_prefix>
Reads <new_out_dir>/*/runs/*/detection_report.yaml; baseline from genomes.tsv and loci.tsv of the v0.6.0 campaign.
"""
import collections, csv, glob, os, statistics as st, sys
import yaml
try: Loader = yaml.CSafeLoader
except AttributeError: Loader = yaml.SafeLoader
NEW, OLD, PFX = sys.argv[1:4]
genomes = [l.split('\t')[0] for l in open(f'{NEW}/list.tsv')]
gmeta = {r['genome']: r for r in csv.DictReader(open(f'{OLD}/genomes.tsv'), delimiter='\t')}
old_loci = collections.defaultdict(list)
for r in csv.DictReader(open(f'{OLD}/loci.tsv'), delimiter='\t'): old_loci[r['genome']].append(r)
new = {}
for p in glob.glob(f'{NEW}/*/runs/*/detection_report.yaml'):
    g = os.path.basename(os.path.dirname(p))
    rep = yaml.load(open(p), Loader=Loader)
    new[g] = [x for x in (rep.get('detected') or [])]
flank = lambda genes: ('both' if 'APN2' in genes and 'SLA2' in genes else 'APN2 only' if 'APN2' in genes else 'SLA2 only' if 'SLA2' in genes else 'neither')
rows = []; rank = {'low': 0, 'medium': 1, 'high': 2}
for g in genomes:
    o = old_loci.get(g, []); n = new.get(g)
    om = gmeta.get(g, {})
    rows.append(dict(genome=g, order=om.get('order', ''), species=om.get('species', ''), old_status=om.get('status', ''),
                     new_report=n is not None, old_n=len(o), new_n=len(n or []),
                     old_idio='+'.join(sorted({x['idiomorph'] for x in o})) or 'none', new_idio='+'.join(sorted({x['idiomorph'] for x in (n or [])})) or 'none',
                     old_conf=max((rank[x['confidence']] for x in o), default=-1), new_conf=max((rank[x['confidence']] for x in (n or [])), default=-1)))
have = [r for r in rows if r['new_report'] and r['old_status'] in ('called', 'uncalled')]
out = open(PFX + '_summary.md', 'w')
def P(s=''): out.write(s + '\n'); print(s)
P(f'# Dothideomycetes full run with records vs v0.6.0 campaign')
P(f'genomes listed {len(genomes)}; reports in new run {sum(r["new_report"] for r in rows)}; compared (report in both) {len(have)}')
oc = sum(r['old_n'] > 0 for r in have); nc = sum(r['new_n'] > 0 for r in have)
P(f'genomes with a call: v0.6.0 {oc} ({100*oc/len(have):.1f}%), with records {nc} ({100*nc/len(have):.1f}%)')
P(f'loci: v0.6.0 {sum(r["old_n"] for r in have)}, with records {sum(r["new_n"] for r in have)}')
gain = [r for r in have if r['old_n'] == 0 and r['new_n'] > 0]; lost = [r for r in have if r['old_n'] > 0 and r['new_n'] == 0]
both = [r for r in have if r['old_n'] > 0 and r['new_n'] > 0]
chg = [r for r in both if r['old_idio'] != r['new_idio']]
up = sum(r['new_conf'] > r['old_conf'] for r in both); down = sum(r['new_conf'] < r['old_conf'] for r in both)
P(f'gained a call {len(gain)}; lost every call {len(lost)}; called in both {len(both)} (idiomorph set changed {len(chg)}; best confidence up {up}, down {down}; locus count changed {sum(r["old_n"] != r["new_n"] for r in both)})')
P('')
P('| order | genomes | called v0.6.0 | called with records | gained | lost | idiomorph changed |'); P('|---|---|---|---|---|---|---|')
bo = collections.defaultdict(list)
for r in have: bo[r['order'] or 'order not assigned'].append(r)
for o, v in sorted(bo.items(), key=lambda x: -len(x[1]))[:12]:
    P(f'| {o} | {len(v)} | {sum(r["old_n"]>0 for r in v)} | {sum(r["new_n"]>0 for r in v)} | {sum(r["old_n"]==0 and r["new_n"]>0 for r in v)} | {sum(r["old_n"]>0 and r["new_n"]==0 for r in v)} | {sum(r["old_n"]>0 and r["new_n"]>0 and r["old_idio"]!=r["new_idio"] for r in v)} |')
P('')
old_l = [x for g in (r['genome'] for r in have) for x in old_loci.get(g, [])]; new_l = [x for g in (r['genome'] for r in have) for x in new.get(g, [])]
cl = lambda L, k: collections.Counter(x[k] for x in L)
P(f'locus class v0.6.0: {dict(cl(old_l, "locus_class"))}'); P(f'locus class with records: {dict(cl(new_l, "locus_class"))}')
P(f'confidence v0.6.0: {dict(cl(old_l, "confidence"))}'); P(f'confidence with records: {dict(cl(new_l, "confidence"))}')
fo = collections.Counter(flank(x['genes_found']) for x in old_l); fn = collections.Counter(flank((x.get('genes_found') or [])) for x in new_l)
P(f'flank genes (APN2/SLA2) v0.6.0: {dict(fo)}'); P(f'flank genes with records: {dict(fn)}')
sz = lambda L: st.median(int(x['end']) - int(x['start']) for x in L) / 1000
P(f'median locus length kb: v0.6.0 {sz(old_l):.1f}, with records {sz(new_l):.1f}')
with open(PFX + '_genome_changes.tsv', 'w') as o:
    o.write('genome\torder\tspecies\tchange\told_idiomorph\tnew_idiomorph\told_loci\tnew_loci\n')
    for r in gain + lost + chg:
        c = 'gained' if r in gain else 'lost' if r in lost else 'idiomorph_changed'
        o.write(f'{r["genome"]}\t{r["order"]}\t{r["species"]}\t{c}\t{r["old_idio"]}\t{r["new_idio"]}\t{r["old_n"]}\t{r["new_n"]}\n')
P(''); P(f'no report in new run: {[r["genome"] for r in rows if not r["new_report"]][:10]}')
