"""Per-genome comparison of the three detect arms:
base = f81dad1 (no new records, no anchors); noanchor = 1403399 with MIP1/beta_fg
excluded from search; anchor = 1403399."""
import yaml, os, collections, statistics
D = os.path.dirname(os.path.abspath(__file__))
pilot = [l.rstrip('\n').split('\t') for l in open(f'{D}/pilot.tsv')]
sub = {'Agaricomycetes':'Agarico','Dacrymycetes':'Agarico','Tremellomycetes':'Agarico','Ustilaginomycetes':'Ustilago',
       'Malasseziomycetes':'Ustilago','Exobasidiomycetes':'Ustilago','Microbotryomycetes':'Puccinio','Pucciniomycetes':'Puccinio',
       'Mixiomycetes':'Puccinio','Wallemiomycetes':'Wallemio'}
ARMS = ['base', 'noanchor', 'anchor']
def rep(arm, asm):
    p = f'{D}/detect_{arm}/runs/{asm}/detection_report.yaml'
    if not os.path.exists(p): return None, None
    r = yaml.safe_load(open(p)); w = open(f'{D}/detect_{arm}/runs/{asm}/wall_seconds').read().strip()
    return r, int(w)
def summ(r):
    if r is None: return 'NO REPORT'
    det = r.get('detected') or []
    hi = [d for d in det if d['confidence'] == 'high']
    fams = collections.Counter(d['family'].split(':')[1] for d in det)
    anch = sum(1 for d in det if {'MIP1', 'beta_fg'} & set(d.get('genes_found', [])))
    return f"{r['routing_mode']}; {len(det)} det ({len(hi)} high); anchored {anch}; {dict(fams)}"
tot = {a: collections.Counter() for a in ARMS}; walls = {a: collections.defaultdict(list) for a in ARMS}
for asm, taxid, cls, order, fam, sp in pilot:
    s = sub.get(cls, '?')
    print(f'== {s} {order} {sp}')
    for a in ARMS:
        r, w = rep(a, asm)
        print(f'   {a:9s} {w}s  {summ(r)}')
        if r is not None:
            det = r.get('detected') or []
            tot[a][(s, 'genomes')] += 1
            tot[a][(s, 'called')] += bool(det)
            tot[a][(s, 'high')] += any(d['confidence'] == 'high' for d in det)
            tot[a][(s, 'HD_high')] += any(d['confidence'] == 'high' and d['family'].endswith(':HD') for d in det)
            tot[a][(s, 'anchored')] += any({'MIP1', 'beta_fg'} & set(d.get('genes_found', [])) for d in det)
            walls[a][s].append(w)
print('\narm\tsubphylum\tgenomes\tcalled\twith_high\tHD_high\tanchored\tmedian_wall_s\ttotal_wall_s')
for a in ARMS:
    for s in ['Agarico', 'Ustilago', 'Puccinio', 'Wallemio']:
        t = tot[a]
        print(f"{a}\t{s}\t{t[(s,'genomes')]}\t{t[(s,'called')]}\t{t[(s,'high')]}\t{t[(s,'HD_high')]}\t{t[(s,'anchored')]}\t"
              f"{statistics.median(walls[a][s]) if walls[a][s] else '-'}\t{sum(walls[a][s])}")
