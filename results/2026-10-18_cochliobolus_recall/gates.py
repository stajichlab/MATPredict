"""Where does each clade genome's true MAT gene (located by the audit's annotation-free six-frame ORF scan) stop in the pipeline?
usage: gates.py CLADE_DIR  -> clade_gates.tsv
Per genome and per gene (alpha-box ORF >= 100 bits on the MAT1-1-1 profile, HMG-box ORF >= 80 bits on the MAT1-2-1 profile; the gap in the
distribution of HMG-box scores is 70 (paralogs) vs 86 (MAT1-2-1), see the note)."""
import csv, glob, json, os, re, sys
import yaml
D = sys.argv[1]
scan = {}
for f in glob.glob('results/2026-10-18_dothideomycetes_audit/scan/scan_*.tsv'):
    for r in csv.DictReader(open(f), delimiter='\t'):
        scan[r['genome']] = r
rows = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_genomes.tsv'), delimiter='\t'))
def fl(x):
    try: return float(x)
    except: return None
def loc(w):
    m = re.match(r'(.+):(\d+)-(\d+)([+-])', (w or '').split(';')[0])
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None
def ov(a0, a1, b0, b1, pad=0):
    return a1 >= b0 - pad and a0 <= b1 + pad
out = []
for r in rows:
    g = r['genome']
    rep = f'{D}/{g}/detection_report.yaml'
    if not os.path.exists(rep): continue
    y = yaml.safe_load(open(rep))
    dg = f'{D}/{g}/evidence_diagnostics.jsonl'
    ev = [json.loads(l) for l in open(dg) if '"kind": "evidence"' in l] if os.path.exists(dg) else []
    s = scan.get(g, {})
    called = y.get('detected') or []
    for gene, key, thr in (('MAT1-1-1', 'alpha_prof', 100.0), ('MAT1-2-1', 'hmg_prof', 80.0)):
        sc = fl(s.get(key + '_best'))
        if sc is None or sc < thr: continue
        l = loc(s.get(key + '_where'))
        rec = dict(genome=g, genus=r['genus'], species=r['species'], strain=r['strain'], clade=r['clade'], status=r['status'], gene=gene,
                   scan_score=sc, scan_where=s.get(key + '_where', '').split(';')[0], n50=y['run']['genome'].get('n50'))
        cl = [e for e in ev if l and e['contig'] == l[0] and ov(e['cluster_start'], e['cluster_end'], l[1], l[2], 300)]
        if not cl:
            rec.update(stage='no_cluster', cluster_identity='', cluster_genes='', admitted='', withheld='')
        else:
            e = max(cl, key=lambda e: e['best_identity'] or 0)
            sup = [x for x in (y.get('suppressed_loci') or []) if x['contig'] == e['contig'] and ov(x['start'], x['end'], e['cluster_start'], e['cluster_end'])]
            cal = [x for x in called if x.get('contig') == e['contig'] and ov(x.get('start', 0), x.get('end', 0), e['cluster_start'], e['cluster_end'])] if called else []
            rec.update(cluster_identity=e['best_identity'], cluster_genes=e['gene_count'], cluster_len=e['cluster_end'] - e['cluster_start'] + 1, admitted=e['admitted'])
            if cal:
                rec.update(stage='called', withheld='')
            elif not e['admitted']:
                rec.update(stage='evidence_floor', withheld='')
            elif sup:
                rec.update(stage='withheld', withheld=sup[0].get('withheld_reason'), sup_genes=','.join(sup[0].get('genes_found') or []), sup_polished=sup[0].get('polished_genes'))
            else:
                rec.update(stage='admitted_not_reported', withheld='')
        out.append(rec)
cols = ['genome', 'genus', 'species', 'strain', 'clade', 'status', 'gene', 'scan_score', 'scan_where', 'n50', 'stage', 'cluster_identity', 'cluster_genes', 'cluster_len', 'admitted', 'withheld', 'sup_genes', 'sup_polished']
with open('results/2026-10-18_cochliobolus_recall/clade_gates.tsv', 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore'); w.writeheader(); w.writerows(out)
from collections import Counter
c = Counter((r['clade'], r['genus'], r['gene'], r['stage'], r.get('withheld', '')) for r in out)
for k in sorted(c): print(k, c[k])
