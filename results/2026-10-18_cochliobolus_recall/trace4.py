"""Trace the four C4/C5 genomes through the audit run's diagnostics: the true cluster, the admitted clusters, the suppressed loci.
usage: trace4.py DIAG_DIR -> trace4.tsv (stdout summary)"""
import json, sys, yaml, csv
D = sys.argv[1]
G = ['GCA_028891585.1_CocheC5_4m', 'GCA_000338975.1_CocheC5_3', 'GCF_028858645.1_CocheC4_6', 'GCF_000354255.1_CocheC4_1']
rows = []
for g in G:
    ev = [json.loads(l) for l in open(f'{D}/{g}/evidence_diagnostics.jsonl') if '"kind": "evidence"' in l]
    y = yaml.safe_load(open(f'{D}/{g}/detection_report.yaml'))
    core = [e for e in ev if 'core_MAT' in e['roles']]
    true = max(core, key=lambda e: e['best_identity'] or 0)
    adm = [e for e in ev if e['admitted']]
    adm_core_max = max((e['best_identity'] or 0) for e in adm)
    sup = y.get('suppressed_loci') or []
    gm = y['run']['genome']
    flank = sorted([e for e in ev if e['contig'] == true['contig'] and 'flanking_conserved' in e['roles'] and (e['best_identity'] or 0) >= 70], key=lambda e: -e['best_identity'])
    row = dict(genome=g, contigs=gm['contigs'], n50=gm['n50'], clusters=len(ev), admitted=len(adm), suppressed=len(sup),
               withheld_reasons=','.join(sorted({s['withheld_reason'] for s in sup})), best_admitted_identity=adm_core_max,
               true_contig=true['contig'], true_start=true['cluster_start'], true_end=true['cluster_end'], true_identity=true['best_identity'],
               true_genes=true['gene_count'], true_hits=true['hit_count'], true_admitted=true['admitted'],
               flank_clusters_same_contig=';'.join(f"{e['cluster_start']}-{e['cluster_end']}:{e['best_identity']}" for e in flank),
               true_in_suppressed=any(s['contig'] == true['contig'] and s['start'] <= true['cluster_end'] and s['end'] >= true['cluster_start'] for s in sup),
               not_detected_reason=(y.get('not_detected') or [{}])[0].get('reason', '')[:100])
    rows.append(row)
    print(row)
with open('results/2026-10-18_cochliobolus_recall/trace4.tsv', 'w') as f:
    w = csv.DictWriter(f, list(rows[0].keys()), delimiter='\t'); w.writeheader(); w.writerows(rows)
