"""Before (full Basidiomycota run, run-ad1f865, old DB) vs after (run-293640d,
new rustHD/redPR records) on the pilot genomes. Writes compare.tsv and prints
a summary. Usage: python3 compare.py"""
import csv, glob, os
import yaml

D = '/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_puccinio_curation'
B = '/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-26_basidiomycota_full'
# genomes whose own assembly supplied a curated record (self-hit, record radius)
SELF = {'GCF_000149925.1_ASM14992v1': 'Pgt CRL 75-36-700-3', 'GCF_000204055.1_v1.0': 'Mlp 98AG31',
        'GCA_000151525.2_P_triticina_1_1_V2': 'Pt 1-1 BBBD Race 1'}


def report(path):
    if not os.path.exists(path):
        return None
    r = yaml.safe_load(open(path))
    calls = []
    for x in r.get('detected') or []:
        genes = ';'.join(f"{g['gene']}:{round(g['identity'] or 0, 1)}:{g['status']}" for g in x.get('gene_evidence', []))
        calls.append(f"{x['family'].split(':')[1]}/{x['idiomorph']}/{x['confidence']}/{x['locus_class']}[{genes}]")
    return {'routing': r.get('routing_mode'), 'calls': calls,
            'withheld': len(r.get('suppressed_loci') or [])}


def wall(d):
    p = os.path.join(d, 'wall_seconds')
    return open(p).read().strip() if os.path.exists(p) else ''


rows = []
for line in open(f'{D}/pilot.tsv'):
    g, taxid, order, species, size, *_ = line.rstrip('\n').split('\t')
    bdir = (glob.glob(f'{B}/*/runs/{g}/') or [''])[0]
    adir = (glob.glob(f'{D}/after_*/runs/{g}/') or [''])[0]
    b = report(os.path.join(bdir, 'detection_report.yaml')) if bdir else None
    a = report(os.path.join(adir, 'detection_report.yaml')) if adir else None
    rows.append({'genome': g, 'order': order, 'species': species, 'size_mb': round(int(size) / 1e6),
                 'self_record': SELF.get(g, ''),
                 'before_routing': b and b['routing'], 'before_called': bool(b and b['calls']),
                 'before_calls': ' | '.join(b['calls']) if b else 'NO REPORT',
                 'before_wall_s': wall(bdir) if bdir else '',
                 'after_routing': a and a['routing'], 'after_called': bool(a and a['calls']),
                 'after_calls': ' | '.join(a['calls']) if a else 'NO REPORT',
                 'after_withheld': a and a['withheld'], 'after_wall_s': wall(adir) if adir else ''})
with open(f'{D}/compare.tsv', 'w', newline='') as fo:
    w = csv.DictWriter(fo, fieldnames=list(rows[0]), delimiter='\t')
    w.writeheader()
    w.writerows(rows)
for order in ('Sporidiobolales', 'Pucciniales'):
    rs = [r for r in rows if r['order'] == order]
    nself = [r for r in rs if not r['self_record']]
    print(f"{order}: n={len(rs)} before called {sum(r['before_called'] for r in rs)}, "
          f"after called {sum(r['after_called'] for r in rs)}; excluding self-record genomes "
          f"(n={len(nself)}): before {sum(r['before_called'] for r in nself)}, after {sum(r['after_called'] for r in nself)}")
for r in rows:
    print(f"{r['genome'][:34]:35} {r['species'][:26]:27} {r['size_mb']:>5} Mb  "
          f"B[{r['before_routing']}] {r['before_calls'][:70]}\n{'':70}A[{r['after_routing']}] {r['after_calls'][:150]} "
          f"withheld={r['after_withheld']} wall {r['before_wall_s']}->{r['after_wall_s']}{'  SELF' if r['self_record'] else ''}")
