"""Neighbourhood synteny of the SLA2-APN2 block with miniprot, across Xylariales and outgroup orders.

Queries: the 14 NC1011 proteins KAI0195597.1-KAI0195610.1 (JAJLYR010000004.1:225-272 kb), in NC1011 order:
 1 ARM  2 FRE  3 APC5  4 CIA30  5 SLA2  6 COX13  7 APN2  8 H604  9 CPN10  10 H606  11 GPR1  12 H608  13 H609  14 H610
For each genome the 150-kb window with most distinct queries is taken; per query the best hit in it.
Writes: genomes_miniprot.tsv (per genome), genes_miniprot.tsv (per gene hit), adjacency_miniprot.tsv.
"""
import csv, os, subprocess, sys, tempfile, gzip, shutil, collections
from multiprocessing import Pool
D = sys.argv[1]; NP = int(sys.argv[2]) if len(sys.argv) > 2 else 8
MP = '/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/miniprot'
LIB = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes'
SAMPLES = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv'
NAMES = ['ARM', 'FRE', 'APC5', 'CIA30', 'SLA2', 'COX13', 'APN2', 'H604', 'CPN10', 'H606', 'GPR1', 'H608', 'H609', 'H610']
QIDS = [f'KAI01955{n}.1' for n in range(97, 100)] + [f'KAI0195{n}.1' for n in range(600, 611)]
QMAP = dict(zip(QIDS, NAMES)); IDX = {n: i + 1 for i, n in enumerate(NAMES)}
WIN = 150000
OUT_ORDERS = ['Sordariales', 'Hypocreales', 'Glomerellales', 'Magnaporthales', 'Ophiostomatales', 'Diaporthales', 'Coniochaetales']
W = f'{D}/miniprot'; os.makedirs(W, exist_ok=True)

def genome_list():
    rows = list(csv.reader(open(SAMPLES)))[1:]
    sel = []; seen = collections.defaultdict(set)
    for r in sorted(rows):
        asm, order, fam, genus, sp = r[0], r[10], r[11], r[12], r[13]
        if order in ('Xylariales', 'Amphisphaeriales'):
            sel.append((asm, order, fam, genus, sp, 'ingroup' if order == 'Xylariales' else 'sister'))
        elif order in OUT_ORDERS and genus and genus not in seen[order] and len(seen[order]) < 10:
            seen[order].add(genus); sel.append((asm, order, fam, genus, sp, 'outgroup'))
    return [s for s in sel if os.path.exists(f'{LIB}/{s[0]}.fa.gz')]

def flip(genes):
    """genes: list of (name,start,end,strand) sorted by start -> reversed orientation."""
    return [(n, s, e, '+' if st == '-' else '-') for n, s, e, st in reversed(genes)]

def run(g):
    asm, order, fam, genus, sp, grp = g
    tmp = tempfile.mkdtemp(prefix='mp_', dir=os.environ.get('SCRATCH', '/tmp'))
    try:
        fa = f'{tmp}/g.fa'
        with gzip.open(f'{LIB}/{asm}.fa.gz', 'rb') as i, open(fa, 'wb') as o: shutil.copyfileobj(i, o)
        p = subprocess.run([MP, '-t2', '-N', '10', '--outs', '0.5', '--outc', '0.3', fa, f'{D}/miniprot/q14.faa'], capture_output=True, text=True)
        hits = []
        for l in p.stdout.splitlines():
            f = l.split('\t'); q = f[0]; qlen = int(f[1]); cov = (int(f[3]) - int(f[2])) / qlen
            tag = {t.split(':')[0]: t.split(':')[2] for t in f[12:] if t.count(':') >= 2}
            pos = int(tag['np']) / qlen; asc = int(tag['AS'])
            if cov >= 0.5 and pos >= 0.30: hits.append((QMAP[q], f[5], int(f[7]), int(f[8]), f[4], pos, cov, asc))
        if not hits: return (g, None, [], 'no_hits')
        # best window
        best = None
        by = collections.defaultdict(list)
        for h in hits: by[h[1]].append(h)
        for ct, hs in by.items():
            hs.sort(key=lambda h: h[2])
            for h0 in hs:
                win = [h for h in hs if h0[2] <= h[2] <= h0[2] + WIN]
                nq = len({h[0] for h in win}); sc = sum(h[7] for h in win)
                if best is None or (nq, sc) > best[0]: best = ((nq, sc), ct, win)
        _, ct, win = best
        pick = {}
        for h in win:
            if h[0] not in pick or h[7] > pick[h[0]][7]: pick[h[0]] = h
        genes = sorted([(n, h[2], h[3], h[4]) for n, h in pick.items()], key=lambda x: x[1])
        # normalise orientation: anchor gene strand '+' (first of APC5, CIA30, SLA2, APN2 present)
        anchor = next((a for a in ['APC5', 'CIA30', 'SLA2', 'APN2'] if a in pick), None)
        if anchor and pick[anchor][4] == '-': genes = flip(genes)
        info = {n: pick[n] for n in pick}
        # HMG-like search in the SLA2-APN2 vicinity
        hmg = ''
        if 'SLA2' in pick and 'APN2' in pick:
            a = min(pick['SLA2'][2], pick['APN2'][2]) - 15000; b = max(pick['SLA2'][3], pick['APN2'][3]) + 15000
            a = max(a, 1)
            r = subprocess.run(['samtools', 'faidx', fa]); 
            reg = subprocess.run(['samtools', 'faidx', fa, f'{ct}:{a}-{b}'], capture_output=True, text=True).stdout
            open(f'{tmp}/reg.fa', 'w').write(reg)
            p2 = subprocess.run([MP, '-t2', '--outs', '0.5', '--outc', '0.5', f'{tmp}/reg.fa', f'{D}/miniprot/hmg_refs.faa'], capture_output=True, text=True)
            hs = []
            for l in p2.stdout.splitlines():
                f = l.split('\t'); qlen = int(f[1])
                tag = {t.split(':')[0]: t.split(':')[2] for t in f[12:] if t.count(':') >= 2}
                hs.append((int(tag['AS']), f[0], int(tag['np']) / qlen, int(f[7]) + a, int(f[8]) + a))
            if hs:
                hs.sort(reverse=True); hmg = '%s;np=%.2f;%d-%d' % (hs[0][1][:40], hs[0][2], hs[0][3], hs[0][4])
            else: hmg = 'none'
        return (g, ct, genes, info, hmg)
    except Exception as e:
        return (g, None, [], 'error:%s' % e)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

def core_class(genes):
    d = {n: st for n, s, e, st in genes}; o = [n for n, s, e, st in genes if n in ('APC5', 'SLA2', 'COX13', 'APN2')]
    if len(o) < 4: return 'incomplete(%d/4)' % len(o)
    s = ' '.join(f'{n}({d[n]})' for n in o)
    return {'APC5(+) COX13(-) APN2(+) SLA2(-)': 'canonical', 'APC5(+) SLA2(+) APN2(-) COX13(+)': 'A', 'APC5(+) SLA2(+) COX13(-) APN2(+)': 'B'}.get(s, 'other: ' + s)

if __name__ == '__main__':
    gl = genome_list(); print('genomes', len(gl), flush=True)
    with Pool(NP) as pool: res = pool.map(run, gl, chunksize=2)
    go = open(f'{D}/genomes_miniprot.tsv', 'w'); gd = open(f'{D}/genes_miniprot.tsv', 'w'); ad = open(f'{D}/adjacency_miniprot.tsv', 'w')
    go.write('asm\torder\tfamily\tgenus\tspecies\tgroup\tcontig\tn_genes\tcore_class\torder_string\tsla2_apn2_gap_bp\thmg_in_vicinity\n')
    gd.write('asm\tgene\tcontig\tstart\tend\tstrand_norm\tpositives\tcoverage\n')
    adj_hdr = [f'{NAMES[i]}-{NAMES[i+1]}' for i in range(13)]
    ad.write('asm\torder\tgroup\t' + '\t'.join(adj_hdr) + '\n')
    for r in res:
        g = r[0]
        if r[1] is None:
            go.write('\t'.join(g) + f'\t\t0\t{r[3] if len(r) > 3 else "fail"}\t\t\t\n'); continue
        _, ct, genes, info, hmg = r
        cls = core_class(genes); os_ = ' '.join(f'{n}({st})' for n, s, e, st in genes)
        gp = ''
        nm = {n: (s, e, st) for n, s, e, st in genes}
        if 'SLA2' in nm and 'APN2' in nm: gp = str(max(nm['SLA2'][0], nm['APN2'][0]) - min(nm['SLA2'][1], nm['APN2'][1]))
        go.write('\t'.join(g) + f'\t{ct}\t{len(genes)}\t{cls}\t{os_}\t{gp}\t{hmg}\n')
        for n, s, e, st in genes:
            h = info[n]; gd.write(f'{g[0]}\t{n}\t{ct}\t{s}\t{e}\t{st}\t{h[5]:.2f}\t{h[6]:.2f}\n')
        # adjacency conservation vs NC1011 order (consecutive detected genes; orientation relative)
        pos = {n: i for i, (n, s, e, st) in enumerate(genes)}; stt = {n: st for n, s, e, st in genes}
        # reference orientation of each query in NC1011 (normalised: APC5 '+')
        ref = {'ARM': '-', 'FRE': '-', 'APC5': '+', 'CIA30': '-', 'SLA2': '+', 'COX13': '-', 'APN2': '+', 'H604': '+', 'CPN10': '+', 'H606': '-', 'GPR1': '-', 'H608': '+', 'H609': '-', 'H610': '-'}
        row = []
        for i in range(13):
            a, b = NAMES[i], NAMES[i + 1]
            if a not in pos or b not in pos: row.append('NA'); continue
            adj = abs(pos[a] - pos[b]) == 1
            same = (stt[a] == ref[a]) == (stt[b] == ref[b])
            row.append('1' if adj and same else '0')
        ad.write(f'{g[0]}\t{g[1]}\t{g[5]}\t' + '\t'.join(row) + '\n')
    print('done')
