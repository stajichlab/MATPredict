"""Would deposited Botryosphaeriales MAT proteins (Diplodia/Sphaeropsis sapinea KF551229.1 and KF551228.1, Neofusicoccum parvum KX505932.1, Phyllosticta KT708823/24)
match the true loci in Botryosphaeriales genomes better than the current references do? tblastn of the deposited proteins; HSPs within 10 kb of the v0.6.0 / full-run call.
Usage: botry_deposit_test.py <worktree> <out_prefix>"""
import collections, csv, glob, gzip, json, os, random, shutil, statistics as st, subprocess, sys, tempfile
W, PFX = sys.argv[1:3]
D = f'{W}/results/2026-10-05_botryosphaeriales_deposits'; O = f'{W}/results/2026-10-05_dothideomycetes_full'
LIB = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes'
random.seed(20261005)
meta = {r['genome']: r for r in csv.DictReader(open(f'{W}/results/2026-10-03_ascomycota_v060/genomes.tsv'), delimiter='\t')}
loci = collections.defaultdict(list)
for r in csv.DictReader(open(f'{W}/results/2026-10-03_ascomycota_v060/loci.tsv'), delimiter='\t'): loci[r['genome']].append((r['contig'], int(r['start']), int(r['end'])))
bot = [g for g, r in meta.items() if r['order'] == 'Botryosphaeriales' and loci.get(g)]
bd = [g for g in bot if meta[g]['species'].startswith('Botryosphaeria dothidea')]
others = collections.defaultdict(list)
for g in bot:
    if g not in bd: others[meta[g]['species']].append(g)
ctl = [random.choice(v) for sp, v in sorted(others.items())]
random.shuffle(ctl); ctl = ctl[:20]
print('B. dothidea genomes', len(bd), '| controls', len(ctl), flush=True)
MAT = ('MAT1-1-1', 'MAT1-1-4', 'MAT1-2-1', 'MAT1-2-5')
def cluster_identity(g, contig, a, b):
    f = glob.glob(f'{O}/small_*/runs/{g}/evidence_diagnostics.jsonl')
    if not f: return None
    best = None
    for l in open(f[0]):
        x = json.loads(l)
        if x['kind'] == 'evidence' and x['family'] == 'Ascomycota:MAT' and x['contig'] == contig and min(b, x['cluster_end']) - max(a, x['cluster_start']) > 0:
            best = max(best or 0, x['best_identity'])
    return best
rows = []
for grp, gl in (('B. dothidea', bd), ('control', ctl)):
    for g in gl:
        tmp = tempfile.mkdtemp(prefix='bt_', dir=os.environ.get('SCRATCH', '/tmp'))
        try:
            fa = f'{tmp}/g.fa'
            with gzip.open(f'{LIB}/{g}.fa.gz', 'rb') as i, open(fa, 'wb') as o: shutil.copyfileobj(i, o)
            subprocess.run(['makeblastdb', '-in', fa, '-dbtype', 'nucl', '-out', f'{tmp}/db'], capture_output=True, check=True)
            p = subprocess.run(['tblastn', '-query', f'{D}/deposited_mat_queries.faa', '-db', f'{tmp}/db', '-evalue', '1e-3', '-num_threads', '4', '-seg', 'no',
                                '-outfmt', '6 qseqid sseqid pident length qlen evalue bitscore sstart send'], capture_output=True, text=True)
            hsps = [l.split('\t') for l in p.stdout.splitlines()]
            for contig, a, b in loci[g]:
                near = [h for h in hsps if h[1] == contig and min(b + 10000, max(int(h[7]), int(h[8]))) - max(a - 10000, min(int(h[7]), int(h[8]))) > 0]
                best = {}
                for h in near:
                    name = h[0].split('|')[1]; src = h[0].split('|')[0]
                    if name in MAT or name.startswith('MAT'):
                        key = f'{src}:{name}'; sc = float(h[6])
                        if key not in best or sc > best[key][0]: best[key] = (sc, float(h[2]), int(h[3]) / int(h[4]))
                top = max(best.items(), key=lambda kv: kv[1][0]) if best else None
                dsap = [v for k, v in best.items() if k.startswith('Dsap')]
                rows.append((grp, g, meta[g]['species'][:30], contig, a, b, cluster_identity(g, contig, a, b), max((v[1] for v in dsap), default=None), top[0] if top else '', round(top[1][1], 1) if top else '', round(top[1][2], 2) if top else ''))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
with open(PFX + '.tsv', 'w') as o:
    o.write('group\tgenome\tspecies\tcontig\tstart\tend\tcurrent_cluster_best_identity\tbest_Dsapinea_pident_at_locus\tbest_deposit_gene\tbest_deposit_pident\tbest_deposit_query_coverage\n')
    for r in rows: o.write('\t'.join('' if x is None else str(x) for x in r) + '\n')
out = open(PFX + '_summary.txt', 'w')
def P(s=''): out.write(s + '\n'); print(s)
for grp in ('B. dothidea', 'control'):
    R = [r for r in rows if r[0] == grp]
    cur = [r[6] for r in R if r[6] is not None]; ds = [r[7] for r in R if r[7] is not None]; top = [r[9] for r in R if r[9] != '']
    P(f'{grp}: {len(R)} loci in {len({r[1] for r in R})} genomes | current references: best cluster identity median {st.median(cur):.1f} (min {min(cur):.1f}) | '
      f'D. sapinea proteins at the locus: found at {len(ds)} loci, best identity median {st.median(ds) if ds else float("nan"):.1f} (min {min(ds) if ds else float("nan"):.1f}) | any deposit: median {st.median(top) if top else float("nan"):.1f}')
P(''); P('B. dothidea loci (current vs deposited):')
for r in [x for x in rows if x[0] == 'B. dothidea'][:14]: P(f'  {r[1][:28]:28s} {r[3]}:{r[4]}-{r[5]}  current {r[6]}  D. sapinea {r[7]}  best deposit {r[8]} {r[9]}% (cov {r[10]})')
