"""Is SLA2 near the MAT locus in Dothideomycetes? Genome-wide miniprot search of SLA2/APN2/COX13 (NC1011 proteins) in genomes with a called
full locus; distance from each best hit to the called locus. Sordariomycetes (SLA2 beside MAT in most) are the control.
Usage: sla2_distance.py <results_dir_of_ascomycota_v060> <queries.faa> <out.tsv> [threads]
"""
import csv, collections, gzip, os, random, shutil, subprocess, sys, tempfile
from multiprocessing import Pool
RES, QF, OUT = sys.argv[1:4]; NP = int(sys.argv[4]) if len(sys.argv) > 4 else 8
MP = '/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/miniprot'
LIB = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes'
NAMES = {'KAI0195599.1': 'APC5', 'KAI0195601.1': 'SLA2', 'KAI0195602.1': 'COX13', 'KAI0195603.1': 'APN2'}
random.seed(20261005)
loci = list(csv.DictReader(open(f'{RES}/loci.tsv'), delimiter='\t'))
sel = []
for cls, per_order, cap in (('Dothideomycetes', 25, 200), ('Sordariomycetes', 6, 60)):
    by = collections.defaultdict(list); seen = set()
    for r in loci:
        if r['class_'] == cls and r['family_called'] == 'Ascomycota:MAT' and r['locus_class'] == 'mat_locus' and r['genome'] not in seen:
            seen.add(r['genome']); by[r['order']].append(r)
    pick = []
    for o, v in by.items():
        random.shuffle(v); pick += v[:per_order]
    random.shuffle(pick); sel += pick[:cap]
sel = [r for r in sel if os.path.exists(f'{LIB}/{r["genome"]}.fa.gz')]
print('genomes', len(sel), collections.Counter(r['class_'] for r in sel), flush=True)

def run(r):
    tmp = tempfile.mkdtemp(prefix='sla2_', dir=os.environ.get('SCRATCH', '/tmp'))
    try:
        fa = f'{tmp}/g.fa'
        with gzip.open(f'{LIB}/{r["genome"]}.fa.gz', 'rb') as i, open(fa, 'wb') as o: shutil.copyfileobj(i, o)
        p = subprocess.run([MP, '-t2', '-N', '10', '--outs', '0.5', '--outc', '0.3', fa, QF], capture_output=True, text=True)
        hits = collections.defaultdict(list)
        for l in p.stdout.splitlines():
            f = l.split('\t'); qlen = int(f[1]); cov = (int(f[3]) - int(f[2])) / qlen
            tag = {t.split(':')[0]: t.split(':')[2] for t in f[12:] if t.count(':') >= 2}
            pos = int(tag['np']) / qlen
            if cov >= 0.5 and pos >= 0.30: hits[NAMES[f[0]]].append((int(tag['AS']), f[5], int(f[7]), int(f[8]), pos))
        out = []
        ls, le, lc = int(r['start']), int(r['end']), r['contig']
        for g in ('SLA2', 'APN2', 'COX13', 'APC5'):
            h = sorted(hits.get(g, []), reverse=True)
            if not h: out.append((g, 'no_hit', '', '', 0, '')); continue
            sc, ct, a, b, pos = h[0]
            if ct != lc: where, d = 'other_contig', ''
            else:
                d = max(a - le, ls - b, 0); where = 'within_locus' if d == 0 else 'same_contig'
            near = ''
            same = [x for x in h if x[1] == lc]
            if same: near = min(max(x[2] - le, ls - x[3], 0) for x in same)
            out.append((g, where, d, near, len(h), round(pos, 2)))
        return r, out
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

if __name__ == '__main__':
    with Pool(NP) as pool: res = pool.map(run, sel, chunksize=2)
    with open(OUT, 'w') as o:
        o.write('genome\tclass\torder\tlocus_contig\tlocus_start\tlocus_end\tgenes_in_call\tgene\twhere\tdist_best_bp\tdist_nearest_any_hit_bp\tn_hits\tpositives\n')
        for r, out in res:
            for g, where, d, near, n, pos in out:
                o.write(f'{r["genome"]}\t{r["class_"]}\t{r["order"]}\t{r["contig"]}\t{r["start"]}\t{r["end"]}\t{r["genes_found"]}\t{g}\t{where}\t{d}\t{near}\t{n}\t{pos}\n')
    print('done', len(res))
