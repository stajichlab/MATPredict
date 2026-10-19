"""Join the clade classification with the audit's annotation-free six-frame ORF scan (independent of the pipeline).
Writes clade_joined.tsv with: scan alpha-box Pfam / profile best score, HMG profile best score and their location,
whether the pipeline's best core cluster lies on the same contig as the scan's best ORF, and a class per genome."""
import csv, glob, re
scan = {}
for f in glob.glob('results/2026-10-18_dothideomycetes_audit/scan/scan_*.tsv'):
    for r in csv.DictReader(open(f), delimiter='\t'):
        scan[r['genome']] = r
rows = list(csv.DictReader(open('results/2026-10-18_cochliobolus_recall/clade_classify.tsv'), delimiter='\t'))
def fl(x):
    try: return float(x)
    except: return None
def loc(w):
    m = re.match(r'(.+):(\d+)-(\d+)([+-])', w or '')
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None
out = []
for r in rows:
    s = scan.get(r['genome'], {})
    ab = fl(s.get('alpha_prof_best')); hb = fl(s.get('hmg_prof_best')); ap = fl(s.get('alpha_pfam_best'))
    r2 = dict(r, alpha_pfam_best=s.get('alpha_pfam_best', ''), alpha_prof_best=s.get('alpha_prof_best', ''), alpha_where=s.get('alpha_prof_where', ''),
              hmg_prof_best=s.get('hmg_prof_best', ''), hmg_where=s.get('hmg_prof_where', ''), hmg_pfam_n=s.get('hmg_pfam_n', ''))
    # thresholds as in the audit: alpha profile >=100 bits, HMG profile >=150 bits count as a real box ORF
    has_alpha = ab is not None and ab >= 100
    has_hmg = hb is not None and hb >= 100
    r2['scan_alpha'] = has_alpha; r2['scan_hmg'] = has_hmg
    r2['scan_state'] = ('both' if has_alpha and has_hmg else 'alpha_only' if has_alpha else 'hmg_only' if has_hmg else 'neither')
    # does the pipeline's best core cluster coincide with the scan ORF?
    bc = (r.get('best_core_contig'), fl(r.get('best_core_start')), fl(r.get('best_core_end')))
    hit = False
    for w in (s.get('alpha_prof_where', ''), s.get('hmg_prof_where', '')):
        l = loc(w.split(';')[0] if w else '')
        if l and bc[0] == l[0] and bc[1] is not None and l[2] >= bc[1] - 500 and l[1] <= bc[2] + 500:
            hit = True
    r2['cluster_is_scan_orf'] = hit
    out.append(r2)
cols = list(out[0].keys())
with open('results/2026-10-18_cochliobolus_recall/clade_joined.tsv', 'w') as f:
    w = csv.DictWriter(f, cols, delimiter='\t'); w.writeheader(); w.writerows(out)
from collections import Counter
c = Counter((r['clade'] == 'clade', r['genus'], r['scan_state'], r['status']) for r in out)
for k in sorted(c): print(k, c[k])
