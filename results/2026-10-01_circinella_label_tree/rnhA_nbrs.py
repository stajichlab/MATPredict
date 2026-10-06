"""Annotated genes within 15 kb of every rnhA miniprot hit (qcov>=0.5), scored with PF00505 and the classifier."""
import csv, sys
gids = sys.argv[1:]
sys.argv = [sys.argv[0]]
exec(open('extract.py').read().split('rows = list(')[0])
with pyhmmer.plan7.HMMFile(str(R / '2026-09-26_sexMP_phylogeny/PF00505.hmm')) as fh:
    hmms['PF00505'] = fh.read()
G = {r['gid']: r for r in csv.DictReader(open('genomes.tsv'), delimiter='\t')}
for gid in gids:
    g = G[gid]
    rn = [m for m in mp_rows(gid) if m['gene'] == 'rnhA' and float(m['qcov']) >= 0.5]
    if not g['gff']:
        print(gid, 'no gff'); continue
    ps = read_fasta(g['proteins'])
    mr = [l.split('\t') for l in open(g['gff']) if '\tmRNA\t' in l]
    for m in rn:
        s, e = int(m['start']), int(m['end'])
        print(f"## {gid} rnhA {m['contig']}:{s}-{e} {m['strand']} qcov={m['qcov']} id={m['identity']}")
        for f in mr:
            if f[0] == m['contig'] and int(f[4]) > s - 15000 and int(f[3]) < e + 15000:
                pid = f[8].split('ID=')[1].split(';')[0]
                sq = ps.get(pid, '')
                sc = score(sq)
                d = max(0, max(int(f[3]), s) - min(int(f[4]), e))
                print(f"   {pid}\t{f[3]}-{f[4]}{f[6]}\td={d/1000:.1f}kb\tlen={len(sq)}\tPF00505={sc['PF00505']}\tsexM={sc['sexM']}\tsexP={sc['sexP']}\tP1={sc['P1']}")
