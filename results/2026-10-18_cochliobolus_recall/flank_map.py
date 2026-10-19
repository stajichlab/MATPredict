"""tblastn the MATPredict roster proteins (MAT core + SLA2/APN2/COX13 flanks) against a genome and report the best hit per roster gene
(identity, query coverage, position, contig length), so that the distance from the MAT gene to its flanks can be read off.
usage: flank_map.py GENOME.fna REFERENCE.faa OUT.tsv GENOME_ID
"""
import subprocess, sys, os, csv
from collections import defaultdict
genome, ref, out, gid = sys.argv[1:5]
os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
clen = {}
k = None
for l in open(genome):
    if l.startswith('>'):
        k = l[1:].split()[0]; clen[k] = 0
    else:
        clen[k] += len(l.strip())
db = out + '.db'
subprocess.run(f'makeblastdb -in {genome} -dbtype nucl -out {db} > /dev/null', shell=True, check=True)
r = subprocess.run(f"tblastn -query {ref} -db {db} -evalue 1e-5 -num_threads 2 -max_target_seqs 500 -outfmt '6 qseqid sseqid sstart send pident length qlen qstart qend bitscore evalue'",
                   shell=True, capture_output=True, text=True, check=True)
best = {}
for l in r.stdout.splitlines():
    f = l.split('\t')
    q = f[0].split('|')
    name = [x for x in q if x.startswith('name=')]
    name = name[0][5:] if name else q[2] if len(q) > 2 else f[0]
    bs = float(f[9])
    row = dict(gene=name, ref=f[0].split('|')[0], contig=f[1], s=min(int(f[2]), int(f[3])), e=max(int(f[2]), int(f[3])), pident=float(f[4]), alen=int(f[5]),
               qlen=int(f[6]), qcov=round(100 * (int(f[8]) - int(f[7]) + 1) / int(f[6]), 1), bits=bs, clen=clen.get(f[1], 0))
    if name not in best or bs > best[name]['bits']:
        best[name] = row
with open(out, 'w') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['genome', 'gene', 'ref', 'contig', 'contig_len', 'start', 'end', 'pident', 'aln_len', 'qlen', 'qcov', 'bits'])
    for n, b in sorted(best.items()):
        w.writerow([gid, n, b['ref'], b['contig'], b['clen'], b['s'], b['e'], b['pident'], b['alen'], b['qlen'], b['qcov'], b['bits']])
for ext in ('.ndb', '.nhr', '.nin', '.njs', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto', '.nsd', '.nsi'):
    try: os.remove(db + ext)
    except OSError: pass
