"""LCG held-out test: per-genome table, Zygo 23 subset score, and a comparison
of detect's sexM/sexP models against the old funannotate annotation.

Read-only: uses reports in runs2/, genomes.tsv, leakage.tsv, annot_hmm_hits.tsv.
Writes per_genome.tsv, zygo23_score.txt, model_vs_annotation.tsv, summary.txt.
"""
import collections, csv, os, re
import yaml
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align

D = os.path.dirname(os.path.abspath(__file__))
G = {r['org']: r for r in csv.DictReader(open(f'{D}/genomes.tsv'), delimiter='\t')}
LK = {r['org']: r for r in csv.DictReader(open(f'{D}/leakage.tsv'), delimiter='\t')}
TRUTH = {}
for l in open('/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-09-23_zygo_bar/zygo_truth.tsv'):
    f = l.rstrip('\n').split('\t')
    TRUTH[f[0]] = (f[1], int(f[2]), int(f[3]), f[4])
HMM = collections.defaultdict(list)
if os.path.exists(f'{D}/annot_hmm_hits.tsv'):
    for r in csv.DictReader(open(f'{D}/annot_hmm_hits.tsv'), delimiter='\t'):
        HMM[r['org']].append(r)

aligner = Align.PairwiseAligner(mode='local', open_gap_score=-10, extend_gap_score=-0.5)
aligner.substitution_matrix = Align.substitution_matrices.load('BLOSUM62')


def pid(a, b):
    """Identity over the local alignment and coverage of the shorter protein."""
    if not a or not b:
        return None, None
    aln = aligner.align(a, b)[0]
    same = alen = 0
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for i in range(e1 - s1):
            alen += 1
            same += a[s1 + i] == b[s2 + i]
    return round(100 * same / alen, 1) if alen else 0.0, round(100 * alen / min(len(a), len(b)), 1)


def gff_genes(path):
    """mRNA spans and their CDS (for mapping detect models to annotated genes)."""
    mrna, cds = {}, collections.defaultdict(list)
    for l in open(path):
        if l.startswith('#'):
            continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9:
            continue
        a = dict(kv.split('=', 1) for kv in f[8].split(';') if '=' in kv)
        if f[2] == 'mRNA':
            mrna[a.get('ID')] = (f[0], int(f[3]), int(f[4]), f[6])
        elif f[2] == 'CDS':
            cds[a.get('Parent')].append((int(f[3]), int(f[4])))
    return mrna, cds


def model_protein(genome_seqs, ev):
    exons = sorted((e['start'], e['end']) for e in (ev.get('exons') or []))
    if not exons or ev.get('contig') not in genome_seqs:
        return ''
    s = genome_seqs[ev['contig']].seq
    nt = ''.join(str(s[a - 1:b]) for a, b in exons)
    if ev.get('strand') == '-':
        nt = str(Seq(nt).reverse_complement())
    best = ''
    for fr in range(3):
        p = str(Seq(nt[fr:len(nt) - ((len(nt) - fr) % 3)]).translate(table=1)).rstrip('*')
        if p.count('*') < (best.count('*') if best else 10**9):
            best = p
    return best


rows, cmp_rows = [], []
for org, g in G.items():
    if g['group'] not in ('Mucorales', 'Umbelopsidales'):
        continue
    rp = f'{D}/runs2/{org}/detection_report.yaml'
    rec = dict(org=org, genus=g['genus'], group=g['group'], leakage=LK.get(org, {}).get('status', ''))
    if not os.path.exists(rp) or os.path.getsize(rp) == 0:
        rec['status'] = 'no_report'
        rows.append(rec)
        continue
    r = yaml.safe_load(open(rp))
    calls = r.get('detected') or []
    rec['status'] = 'called' if calls else 'not_called'
    rec['n_calls'] = len(calls)
    wall = f'{D}/runs2/{org}/wall_seconds'
    rec['wall_s'] = open(wall).read().strip() if os.path.exists(wall) else ''
    parts = []
    for c in calls:
        clf = c.get('idiomorph_classifier') or {}
        parts.append('|'.join(str(v) for v in [
            c.get('contig'), c.get('start'), c.get('end'), c.get('idiomorph'), c.get('confidence'),
            c.get('locus_class'), clf.get('classifier_input', ''),
            clf.get('margin', ''),
            'split' if c.get('split_locus') else '', ','.join(c.get('genes_found') or [])]))
    rec['calls'] = ' ;; '.join(parts)
    rec['gate_withheld'] = r.get('suppressed_mat_gene_gate', 0)
    rec['suppressed_reasons'] = ','.join(sorted({str(x.get('reason') or x.get('withheld_reason') or '') for x in (r.get('suppressed_loci') or [])}))
    rec['idiomorphs'] = '+'.join(sorted({str(c.get('idiomorph')) for c in calls}))
    # Zygo 23 subset
    if org in TRUTH:
        sc, ts, te, tid = TRUTH[org]
        hit = [c for c in calls if c.get('contig') == sc and c.get('start') <= te and c.get('end') >= ts]
        rec['zygo_locus_ok'] = bool(hit)
        rec['zygo_idiomorph_ok'] = bool(hit) and any(c.get('idiomorph') == tid for c in hit)
        rec['zygo_truth'] = tid
    # model vs old annotation
    if calls and g.get('gff3') and os.path.exists(g['gff3']):
        try:
            gseq = SeqIO.to_dict(SeqIO.parse(g['genome'], 'fasta'))
            prot = SeqIO.to_dict(SeqIO.parse(g['proteins'], 'fasta'))
            mrna, cds = gff_genes(g['gff3'])
        except Exception as exc:
            gseq, prot, mrna, cds = {}, {}, {}, {}
            rec['annot_error'] = str(exc)[:80]
        hmmscore = {h['protein']: h for h in HMM.get(org, [])}
        for c in calls:
            for ev in c.get('gene_evidence') or []:
                if ev.get('gene') not in ('sexM', 'sexP') or ev.get('status', '').startswith('not_'):
                    continue
                mp = model_protein(gseq, ev)
                ov = [(mid, m) for mid, m in mrna.items()
                      if m[0] == ev.get('contig') and m[1] <= ev['end'] and m[2] >= ev['start']]
                best = None
                for mid, m in ov:
                    ap = str(prot[mid].seq).rstrip('*') if mid in prot else ''
                    i, cv = pid(mp, ap)
                    if best is None or (i or 0) > (best[1] or 0):
                        best = (mid, i, cv, len(ap), len(cds.get(mid, [])), hmmscore.get(mid))
                h = best[5] if best else None
                cmp_rows.append(dict(
                    org=org, gene=ev['gene'], status=ev.get('status'), contig=ev.get('contig'),
                    start=ev.get('start'), end=ev.get('end'), model_len=len(mp),
                    model_exons=len(ev.get('exons') or []), call_idiomorph=c.get('idiomorph'),
                    annotated=bool(ov), annot_id=best[0] if best else '',
                    annot_len=best[3] if best else '', annot_exons=best[4] if best else '',
                    identity=best[1] if best else '', coverage=best[2] if best else '',
                    annot_classifier=(h['best'] if h else ''), annot_margin=(h['margin_P_minus_M'] if h else '')))
    rows.append(rec)

keys = sorted({k for r in rows for k in r}, key=lambda k: ['org', 'genus', 'group', 'leakage', 'status', 'n_calls', 'idiomorphs', 'calls'].index(k) if k in ['org', 'genus', 'group', 'leakage', 'status', 'n_calls', 'idiomorphs', 'calls'] else 99)
with open(f'{D}/per_genome.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys, delimiter='\t'); w.writeheader(); w.writerows(rows)
if cmp_rows:
    with open(f'{D}/model_vs_annotation.tsv', 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(cmp_rows[0]), delimiter='\t'); w.writeheader(); w.writerows(cmp_rows)

# summary
S = []
c = collections.Counter(r['status'] for r in rows); S.append(f'status: {dict(c)}')
for sub in ('clean', 'in_BFD_not_training', 'training_leak', 'zygo23'):
    rs = [r for r in rows if r['leakage'] == sub]
    S.append(f"{sub}: n={len(rs)} called={sum(r['status']=='called' for r in rs)} "
             f"idiomorphs={dict(collections.Counter(r.get('idiomorphs','') for r in rs if r['status']=='called'))}")
z = [r for r in rows if 'zygo_locus_ok' in r]
S.append(f"Zygo23: n={len(z)} locus_ok={sum(r['zygo_locus_ok'] for r in z)} idiomorph_ok={sum(r['zygo_idiomorph_ok'] for r in z)}")
for r in z:
    if not r['zygo_idiomorph_ok']:
        S.append(f"   zygo miss: {r['org']} truth={r['zygo_truth']} calls={r.get('calls','')}")
if cmp_rows:
    ann = [x for x in cmp_rows if x['annotated']]
    S.append(f"models: {len(cmp_rows)}; annotated at locus: {len(ann)}; identity>=95 & cov>=90: "
             f"{sum(1 for x in ann if (x['identity'] or 0)>=95 and (x['coverage'] or 0)>=90)}; "
             f"annotated protein classifier agrees with call: "
             f"{sum(1 for x in ann if x['annot_classifier'] and ((x['annot_classifier']=='sexP')==(x['call_idiomorph']=='Plus')) and x['call_idiomorph'] in ('Plus','Minus'))}"
             f" of {sum(1 for x in ann if x['annot_classifier'] and x['call_idiomorph'] in ('Plus','Minus'))}")
    S.append(f"  unannotated models (annotation gap): {len(cmp_rows)-len(ann)}")
# annotated sexM/sexP-like proteins in genomes with no call
nocall = {r['org'] for r in rows if r['status'] == 'not_called'}
strong = [h for o in nocall for h in HMM.get(o, []) if max(float(h['sexM_bits']), float(h['sexP_bits'])) >= 100]
S.append(f"not-called genomes with an annotated protein scoring >=100 bits to sexM/sexP: "
         f"{len({h['org'] for h in strong})} ({len(strong)} proteins)")
wall = [int(r['wall_s']) for r in rows if r.get('wall_s', '').isdigit()]
if wall:
    wall.sort(); S.append(f"wall s: median {wall[len(wall)//2]} max {wall[-1]}")
open(f'{D}/summary.txt', 'w').write('\n'.join(S) + '\n')
print('\n'.join(S))
