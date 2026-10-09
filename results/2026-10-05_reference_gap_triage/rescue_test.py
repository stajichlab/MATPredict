"""Would deposited MAT loci rescue genomes the v0.6.0 campaign did not call? tblastn of the deposited core and flank proteins against uncalled good-assembly genomes
of Pichiales, Dipodascales, Chaetothyriales and Microascales, plus called controls. Core hit: identity >= 35% over >= 50% of the query. Rescue candidate: an uncalled genome with a core
hit and a flank hit (identity >= 35%, >= 50% of the query) within 30 kb on the same contig. Controls: a called genome whose core hit lies within 10 kb of the call (concordance).
Usage: rescue_test.py <worktree> <out_prefix>  (python3.12 with pyarrow; blast+ on PATH)"""
import collections, csv, gzip, os, random, shutil, statistics as st, subprocess, sys, tempfile
from multiprocessing import Pool
import pyarrow.parquet as pq
W, PFX = sys.argv[1:3]
D = f'{W}/results/2026-10-05_reference_gap_triage'; LIB = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes'
T = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/tables/'
busco = {r['ASMID']: r['complete_pct'] for r in pq.read_table(T + 'busco_genome.parquet', columns=['ASMID', 'complete_pct']).to_pylist()}
asm = {r['ASMID']: r for r in pq.read_table(T + 'asm_stats.parquet', columns=['ASMID', 'contig_count', 'N50_bp']).to_pylist()}
def poor(g):
    b = busco.get(g); a = asm.get(g)
    if b is None and a is None: return True
    return (b is not None and b < 70) or (a is not None and (a['N50_bp'] < 20000 or a['contig_count'] > 5000))
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
G = rd(f'{W}/results/2026-10-03_ascomycota_v060/genomes.tsv'); loci = collections.defaultdict(list)
for r in rd(f'{W}/results/2026-10-03_ascomycota_v060/loci.tsv'): loci[r['genome']].append((r['contig'], int(r['start']), int(r['end'])))
ORDERS = ('Pichiales', 'Dipodascales', 'Chaetothyriales', 'Microascales')
random.seed(20261005); jobs = []
for o in ORDERS:
    v = [r for r in G if r['order'] == o and not poor(r['genome'])]
    unc = [r for r in v if r['status'] != 'called']; cal = [r for r in v if r['status'] == 'called']; random.shuffle(cal)
    jobs += [(o, r['genome'], r['family'], 'uncalled') for r in unc] + [(o, r['genome'], r['family'], 'called') for r in cal[:30]]
print('genomes', len(jobs), collections.Counter((j[0], j[3]) for j in jobs), flush=True)
def run(job):
    o, g, fam, st_ = job; tmp = tempfile.mkdtemp(prefix='rs_', dir=os.environ.get('SCRATCH', '/tmp'))
    try:
        fa = f'{tmp}/g.fa'
        with gzip.open(f'{LIB}/{g}.fa.gz', 'rb') as i, open(fa, 'wb') as out: shutil.copyfileobj(i, out)
        subprocess.run(['makeblastdb', '-in', fa, '-dbtype', 'nucl', '-out', f'{tmp}/db'], capture_output=True, check=True)
        p = subprocess.run(['tblastn', '-query', f'{D}/rescue_queries.faa', '-db', f'{tmp}/db', '-evalue', '1e-5', '-num_threads', '2', '-seg', 'no',
                            '-outfmt', '6 qseqid sseqid pident length qlen bitscore sstart send'], capture_output=True, text=True)
        core = []; flank = []
        for l in p.stdout.splitlines():
            q, s, pid, ln, ql, bs, a, b = l.split('\t'); pid = float(pid); cov = int(ln) / int(ql)
            if pid < 35 or cov < 0.5: continue
            h = (s, min(int(a), int(b)), max(int(a), int(b)), pid, float(bs), q.split('|')[0], q.split('|')[1])
            (core if q.endswith('|core') else flank).append(h)
        best = None; fl = False
        for h in sorted(core, key=lambda h: -h[4]):
            if best is None: best = h
            if any(f[0] == h[0] and min(f[2], h[2] + 30000) - max(f[1], h[1] - 30000) > 0 for f in flank): best = h; fl = True; break
        conc = None
        if st_ == 'called' and core:
            conc = any(h[0] == c and min(b + 10000, h[2]) - max(a - 10000, h[1]) > 0 for h in core for c, a, b in loci[g])
        return (o, g, fam, st_, len(core), best[3] if best else None, best[5] if best else '', best[6] if best else '', fl, conc)
    except Exception as e:
        return (o, g, fam, st_, -1, None, 'error', str(e)[:40], False, None)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
if __name__ == '__main__':
    with Pool(8) as pool: res = pool.map(run, jobs, chunksize=2)
    with open(PFX + '.tsv', 'w') as o:
        o.write('order\tgenome\tfamily\tstatus\tn_core_hits\tbest_core_pident\tdeposit\tgene\tflank_within_30kb\tconcordant_with_call\n')
        for r in res: o.write('\t'.join('' if x is None else str(x) for x in r) + '\n')
    out = open(PFX + '_summary.txt', 'w')
    def P(s=''): out.write(s + '\n'); print(s)
    P('Rescue test: deposited MAT proteins vs uncalled good-assembly genomes (core hit >= 35% identity over >= 50% of the query; flank within 30 kb)')
    P('| order | uncalled genomes | with a core hit | core hit and a flank (rescue candidates) | median best identity of candidates | called controls | control core hit | control concordant with the call |'); P('|---|---|---|---|---|---|---|---|')
    for o in ORDERS:
        u = [r for r in res if r[0] == o and r[3] == 'uncalled']; c = [r for r in res if r[0] == o and r[3] == 'called']
        uc = [r for r in u if r[4] > 0]; cand = [r for r in u if r[4] > 0 and r[8]]; cc = [r for r in c if r[4] > 0]; conc = [r for r in cc if r[9]]
        P(f'| {o} | {len(u)} | {len(uc)} | {len(cand)} | {st.median([r[5] for r in cand]):.1f} | {len(c)} | {len(cc)} | {len(conc)} |' if cand else f'| {o} | {len(u)} | {len(uc)} | 0 | n/a | {len(c)} | {len(cc)} | {len(conc)} |')
    P(''); P('Candidates by family and deposit:')
    for o in ORDERS:
        cand = [r for r in res if r[0] == o and r[3] == 'uncalled' and r[4] > 0 and r[8]]
        P(f'  {o}: ' + '; '.join(f'{f} {n} ({", ".join(sorted({r[6] for r in cand if r[2] == f}))})' for f, n in collections.Counter(r[2] for r in cand).most_common(6)))
    P(f'errors: {sum(1 for r in res if r[4] < 0)}')
