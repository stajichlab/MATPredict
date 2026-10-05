import glob, os, sys, csv, yaml
try: L = yaml.CSafeLoader
except AttributeError: L = yaml.SafeLoader
W = sys.argv[1]; T = f'{W}/results/2026-10-05_dothideo_cap0_lost'
old = {}
for r in csv.DictReader(open(f'{W}/results/2026-10-03_ascomycota_v060/loci.tsv'), delimiter='\t'): old.setdefault(r['genome'], []).append(r)
sp = {r['genome']: r['species'] for r in csv.DictReader(open(f'{W}/results/2026-10-03_ascomycota_v060/genomes.tsv'), delimiter='\t')}
rec = back = 0
for g in [l.strip() for l in open(f'{T}/lost.txt')]:
    f = glob.glob(f'{T}/out/runs/{g}/detection_report.yaml')
    o = old.get(g, [])
    if not f: print(f'{g[:30]:30s} {sp[g][:26]:26s} not finished'); continue
    rep = yaml.load(open(f[0]), Loader=L); det = rep.get('detected') or []
    rec += 1; back += bool(det)
    d = det[0] if det else None
    print(f'{g[:30]:30s} {sp[g][:26]:26s} v0.6.0 {o[0]["idiomorph"] if o else "none":6s} -> cap off: ' + (f'{d["idiomorph"]} {d["confidence"]} {d["locus_class"]} {d["contig"]}:{d["start"]}-{d["end"]}' if d else f'no call (suppressed {len(rep.get("suppressed_loci") or [])})'))
print(f'\nfinished {rec} of 15; call recovered in {back}')
