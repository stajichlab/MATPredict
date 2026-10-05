"""Nucleotide-level breakpoints between congeneric genomes that differ in SLA2/COX13/APN2 order.

Pairs: query = SCA-like genome, target = SAC genome. Windows come from the miniprot gene hits of the 14
neighbourhood genes. minimap2 -x asm20 --cs; blocks and strand changes are written to breakpoints/<pair>.*
"""
import csv, gzip, os, shutil, subprocess, sys, collections
D = sys.argv[1]
LIB = '/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin'
GL = '/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes'
PAIRS = [('GCA_022453505.1_Xylcub1', 'GCA_000966885.1_ASM96688v1'), ('GCA_022495175.1_Hyparg1', 'GCF_022984875.1_Hypfra2')]
PRESETS = [('asm20', []), ('sensitive', ['-k', '13', '-w', '5', '-s', '100', '-n', '1'])]
O = f'{D}/breakpoints'; os.makedirs(O, exist_ok=True)
genes = collections.defaultdict(dict)
for l in csv.DictReader(open(f'{D}/genes_miniprot.tsv'), delimiter='\t'): genes[l['asm']][l['gene']] = (l['contig'], int(l['start']), int(l['end']))

def window(asm, pad=5000):
    g = genes[asm]; ct = g['SLA2'][0]
    xs = [v for v in g.values() if v[0] == ct]
    return ct, max(1, min(x[1] for x in xs) - pad), max(x[2] for x in xs) + pad

def extract(asm, tmp):
    fa = f'{tmp}/{asm}.fa'
    if not os.path.exists(fa):
        with gzip.open(f'{GL}/{asm}.fa.gz', 'rb') as i, open(fa, 'wb') as o: shutil.copyfileobj(i, o)
        subprocess.run(['samtools', 'faidx', fa], check=True)
    ct, a, b = window(asm)
    seq = subprocess.run(['samtools', 'faidx', fa, f'{ct}:{a}-{b}'], capture_output=True, text=True).stdout
    p = f'{O}/{asm}.window.fa'; open(p, 'w').write(seq); return p, ct, a, b

tmp = os.environ.get('SCRATCH', '/tmp')
summ = open(f'{O}/summary.md', 'w')
for q, t in PAIRS:
    qp, qct, qa, qb = extract(q, tmp); tp, tct, ta, tb = extract(t, tmp)
    for pname, extra in PRESETS:
      paf = f'{O}/{q}_vs_{t}.{pname}.paf'
      r = subprocess.run(['minimap2', '-x', 'asm20'] + extra + ['-c', '--cs', '-N', '5', tp, qp], capture_output=True, text=True)
      open(paf, 'w').write(r.stdout)
      blocks = []
      for l in r.stdout.splitlines():
          f = l.split('\t'); 
          if int(f[10]) < 500: continue
          blocks.append((int(f[7]) + ta, int(f[8]) + ta, f[4], int(f[2]) + qa, int(f[3]) + qa, int(f[9]), int(f[10]), int(f[11])))
      blocks.sort()
      summ.write(f'## [{pname}] query {q} ({qct}:{qa}-{qb}) vs target {t} ({tct}:{ta}-{tb})\n\n')
      summ.write('Target = SAC genome, query = SCA-like genome. Blocks >= 500 bp, ordered along the target.\n\n')
      summ.write('| target start | target end | strand | query start | query end | matches | aln len | mapq |\n|---|---|---|---|---|---|---|---|\n')
      for b in blocks: summ.write('| ' + ' | '.join(str(x) for x in b) + ' |\n')
      summ.write('\nGene hits (miniprot) in each genome:\n\n')
      for name, asm in (('target', t), ('query', q)):
          gg = sorted(genes[asm].items(), key=lambda kv: kv[1][1])
          summ.write(f'- {name} {asm}: ' + ', '.join(f'{k}:{v[1]}-{v[2]}' for k, v in gg) + '\n')
      summ.write('\n')
      summ.flush()
print('done')
