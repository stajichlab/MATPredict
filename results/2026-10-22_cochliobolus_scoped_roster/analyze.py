"""Compare variants. usage: analyze.py runs.tsv scope.tsv clade.txt control94.txt clade_gates.tsv audit_final.tsv OUTDIR"""
import csv, re, sys, os
PROTO = os.environ.get('PROTO', 'proto')
from collections import Counter, defaultdict
runs_p, scope_p, clade_p, ctrl_p, gates_p, audit_p, outdir = sys.argv[1:8]
os.makedirs(outdir, exist_ok=True)
R = defaultdict(dict)
for r in csv.DictReader(open(runs_p), delimiter='\t'):
    R[r['variant']][r['genome']] = r
scope = {r['genome']: r['in_pleosporaceae'] == 'True' for r in csv.DictReader(open(scope_p), delimiter='\t')}
clade = open(clade_p).read().split(); ctrl = open(ctrl_p).read().split()
aud = {r['genome']: r for r in csv.DictReader(open(audit_p), delimiter='\t')}
def mat_calls(r):
    out = []
    for c in (r['calls'].split(';') if r and r['calls'] else []):
        f = c.split('|')
        out.append(dict(family=f[0], contig=f[1], start=int(f[2]), end=int(f[3]), idiomorph=f[4], confidence=f[5], cls=f[6], genes=f[7]))
    return out
def short(r):
    return ';'.join(f"{c['family'].split(':')[1]}:{c['idiomorph']}/{c['confidence']}/{c['cls']}/{c['genes']}" for c in mat_calls(r)) or 'none'
L = []
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); L.append(s)
for v in R: log('variant', v, 'genomes', len(R[v]), Counter(r['status'] for r in R[v].values()))
# ---- truth: scan ORFs + C4/C5 strains
exp = defaultdict(dict); meta = {}
for g in csv.DictReader(open(gates_p), delimiter='\t'):
    idm = 'MAT1-1' if g['gene'] == 'MAT1-1-1' else 'MAT1-2'
    m = re.match(r'(.+):(\d+)-(\d+)', g['scan_where'])
    exp[g['genome']][idm] = (m.group(1), int(m.group(2)), int(m.group(3)))
    meta[g['genome']] = (g['genus'], g['species'], g['strain'])
strain_truth = {'GCA_000338975.1_CocheC5_3': 'MAT1-1', 'GCA_028891585.1_CocheC5_4m': 'MAT1-1',
                'GCF_000354255.1_CocheC4_1': 'MAT1-2', 'GCF_028858645.1_CocheC4_6': 'MAT1-2'}
def evaluate(g, r):
    e = exp.get(g, {}); cs = [c for c in mat_calls(r) if c['family'] == 'Ascomycota:MAT']
    if not cs: return 'no_call'
    bad = []
    for c in cs:
        if c['idiomorph'] not in e: bad.append(c['idiomorph']); continue
        ctg, s, en = e[c['idiomorph']]
        if not (c['contig'] == ctg and c['end'] >= s - 500 and c['start'] <= en + 500): bad.append(c['idiomorph'] + '@elsewhere')
    if bad: return 'DISAGREES_WITH_SCAN:' + ','.join(bad)
    return 'agrees_with_scan' + ('' if {c['idiomorph'] for c in cs} == set(e) else '(subset)')
rows = []
for g in clade:
    b, p, s = R['base'].get(g), R[PROTO].get(g), R['ssg'].get(g)
    gen, sp, strain = meta.get(g, ('', '', ''))
    rows.append(dict(genome=g, genus=gen, species=sp, strain=strain, in_scope=scope.get(g), audit_status=aud.get(g, {}).get('status'),
                     scan=','.join(sorted(exp.get(g, {}))), strain_truth=strain_truth.get(g, ''),
                     base=short(b), proto=short(p), proto_eval=evaluate(g, p) if p else 'NO_RUN', ssg=short(s), ssg_eval=evaluate(g, s) if s else 'NO_RUN'))
with open(f'{outdir}/clade_recovery.tsv', 'w') as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
unc = [r for r in rows if r['audit_status'] != 'called']
log('clade genomes', len(rows), 'audit-uncalled', len(unc), 'base-uncalled among them', sum(1 for r in unc if r['base'] == 'none'))
for var in ('proto', 'ssg'):
    called = [r for r in unc if r[var] != 'none']
    log(f'{var}: called {len(called)} of {len(unc)};', 'by genus', Counter((r['genus'], r[var] != 'none') for r in unc))
    log(f'{var}: eval', Counter(r[f'{var}_eval'] for r in unc))
    log(f'{var}: confidence/class', Counter(tuple(sorted(c.split('/')[1] + '/' + c.split('/')[2] for c in r[var].split(';'))) for r in called))
    log(f'{var}: strain truth', [(r['genome'], r['strain_truth'], r[var]) for r in rows if r['strain_truth']])
# ---- base vs proto
def compare(va, vb, gs, name):
    res = Counter(); ch = []
    for g in gs:
        a, b = R[va].get(g), R[vb].get(g)
        if not a or not b or a['status'] == 'no_report' or b['status'] == 'no_report':
            res['missing_run'] += 1; continue
        if a['report_sha'] == b['report_sha']: res['report_identical'] += 1
        elif a['calls_sha'] == b['calls_sha']: res['calls_identical_report_differs'] += 1
        else:
            res['CALLS_CHANGED'] += 1; ch.append((g, short(a), short(b)))
    log(name, len(gs), dict(res))
    return ch
allg = sorted(set(R['base']) | set(R[PROTO]))
out_scope = [g for g in allg if not scope.get(g)]
in_scope = [g for g in allg if scope.get(g)]
PREV = os.environ.get('PREV', 'proto')
ch_v1v2 = compare(PREV, PROTO, allg, PREV + ' vs ' + PROTO) if PROTO != PREV else []
ch_out = compare('base', PROTO, out_scope, 'OUTSIDE scope base vs proto')
ch_in = compare('base', PROTO, [g for g in in_scope if g not in clade], 'INSIDE scope (non-clade) base vs proto')
ch_inc = compare('base', PROTO, [g for g in in_scope if g in clade], 'INSIDE scope clade base vs proto')
ch_ctrl = compare('base', PROTO, ctrl, 'control94 base vs proto')
sub = sorted(set(R['off']))
compare('base', 'off', sub, 'policy-off branch vs main')
compare('base', 'baserep', sorted(set(R['baserep'])), 'main replicate vs main')
ch_ssg = compare('base', 'ssg', sorted(set(R['ssg'])), 'ssg vs base (its set)')
ch_ssg_ctrl = compare('base', 'ssg', ctrl, 'control94 base vs ssg')
def classify(a, b):
    def parse(s):
        return [] if s == 'none' else [x.split('/') for x in s.split(';')]
    pa, pb = parse(a), parse(b)
    if not pa and pb: return 'NEW_CALL'
    if pa and not pb: return 'LOST_CALL'
    if sorted(x[0] for x in pa) != sorted(x[0] for x in pb): return 'IDIOMORPH_SET_CHANGED'
    if sorted(x[1] for x in pa) != sorted(x[1] for x in pb): return 'CONFIDENCE_CHANGED'
    if sorted(x[2] for x in pa) != sorted(x[2] for x in pb): return 'CLASS_CHANGED'
    return 'GENES_OR_COORDS_ONLY'
for name, ch in (('v1_vs_v2', ch_v1v2), ('outside', ch_out), ('inside_nonclade', ch_in), ('inside_clade', ch_inc), ('ssg', ch_ssg)):
    with open(f'{outdir}/changes_{name}.tsv', 'w') as f:
        f.write('genome\tgenus\tkind\tbefore\tafter\n')
        for g, a, b in ch:
            f.write(f"{g}\t{aud.get(g, {}).get('genus', '')}\t{classify(a, b)}\t{a}\t{b}\n")
    log(f'changes {name}', Counter(classify(a, b) for _, a, b in ch))
open(f'{outdir}/summary.txt', 'w').write('\n'.join(L) + '\n')
