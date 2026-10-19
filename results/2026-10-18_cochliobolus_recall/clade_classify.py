"""Classify each Cochliobolus-clade genome from the audit's detection_report.yaml + evidence_diagnostics.jsonl.
usage: clade_classify.py CLADE_DIR clade_genomes.tsv out.tsv
"""
import json, sys, csv, os
import yaml

d, lst, out = sys.argv[1:4]
rows = list(csv.DictReader(open(lst), delimiter='\t'))
res = []
for r in rows:
    g = r['genome']
    rep = f'{d}/{g}/detection_report.yaml'
    dg = f'{d}/{g}/evidence_diagnostics.jsonl'
    if not os.path.exists(rep):
        continue
    y = yaml.safe_load(open(rep))
    ev = []
    if os.path.exists(dg):
        ev = [json.loads(l) for l in open(dg) if '"kind": "evidence"' in l]
    core = [e for e in ev if 'core_MAT' in e['roles']]
    best = max(core, key=lambda e: e['best_identity'] or 0) if core else None
    gn = y['run']['genome']
    called = y.get('detected') or []
    row = dict(genome=g, genus=r['genus'], species=r['species'], strain=r['strain'], clade=r['clade'], status=r['status'],
               contigs=gn.get('contigs'), n50=gn.get('n50'), length=gn.get('length_bp'),
               n_called=len(called), n_suppressed=y.get('suppressed_unpolished'),
               n_clusters=len(ev), n_admitted=sum(1 for e in ev if e['admitted']))
    if best:
        row.update(best_core_id=best['best_identity'], best_core_contig=best['contig'], best_core_start=best['cluster_start'],
                   best_core_end=best['cluster_end'], best_core_len=best['cluster_end'] - best['cluster_start'] + 1,
                   best_core_genes=best['gene_count'], best_core_admitted=best['admitted'])
        # strongest flank on the same contig
        fl = [e for e in ev if e['contig'] == best['contig'] and 'flanking_conserved' in e['roles'] and e is not best]
        if fl:
            f = max(fl, key=lambda e: (e['best_identity'] or 0))
            dist = min(abs(f['cluster_start'] - best['cluster_end']), abs(best['cluster_start'] - f['cluster_end']))
            row.update(flank_id=f['best_identity'], flank_dist=dist)
            near = [e for e in fl if min(abs(e['cluster_start'] - best['cluster_end']), abs(best['cluster_start'] - e['cluster_end'])) <= 300000 and (e['best_identity'] or 0) >= 70]
            row['n_flank70_within300kb'] = len(near)
        # sup. loci
    nd = (y.get('not_detected') or [{}])[0]
    row['not_detected_reason'] = (nd.get('reason') or '')[:120]
    row['best_fraction'] = nd.get('best_fraction_found')
    row['withheld_reasons'] = ','.join(sorted({s.get('withheld_reason') or '' for s in (y.get('suppressed_loci') or [])}))
    res.append(row)
cols = ['genome', 'genus', 'species', 'strain', 'clade', 'status', 'contigs', 'n50', 'length', 'n_called', 'n_suppressed', 'n_clusters', 'n_admitted',
        'best_core_id', 'best_core_contig', 'best_core_start', 'best_core_end', 'best_core_len', 'best_core_genes', 'best_core_admitted',
        'flank_id', 'flank_dist', 'n_flank70_within300kb', 'best_fraction', 'withheld_reasons', 'not_detected_reason']
with open(out, 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t', extrasaction='ignore')
    w.writeheader(); w.writerows(res)
print(len(res), 'genomes')
