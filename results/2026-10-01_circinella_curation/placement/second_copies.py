"""sexP proteins of gained Circinella calls that have no tip in the tree:
C. minor CBS 143.56 (its Plus locus; its tree tip is the sexM) and the second
sexP copy of C. muscae NRRL 1355/1363/2403. Detect's polished model."""
import yaml
from Bio.Seq import Seq
L = "/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/genomes"
R = "/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-01_circinella_curation/lcg/cand"
def contig(p, n):
    on = False; b = []
    for l in open(p):
        if l[0] == '>':
            if on: break
            on = l[1:].split()[0] == n
        elif on: b.append(l.strip())
    return ''.join(b)
todo = [('Circinella_minor_CBS_143.56', 'scaffold_1115'), ('Circinella_muscae_NRRL_1355', 'scaffold_2644'),
        ('Circinella_muscae_NRRL_1363', 'scaffold_1448'), ('Circinella_muscae_NRRL_2403', 'scaffold_1271')]
with open('second_copies.faa', 'w') as o:
    for g, c in todo:
        d = [d for d in yaml.safe_load(open(f'{R}/{g}/detection_report.yaml'))['detected'] if d['contig'] == c][0]
        e = [e for e in d['gene_evidence'] if e['gene'] == 'sexP'][0]
        cs = contig(f'{L}/{g}.sorted.fasta', c)
        nt = ''.join(cs[s - 1:t] for s, t in sorted((x['start'], x['end']) for x in e['exons']))
        if e['strand'] == '-': nt = str(Seq(nt).reverse_complement())
        p = str(Seq(nt[:len(nt) // 3 * 3]).translate()).rstrip('*')
        assert '*' not in p, g
        o.write(f'>COPY__{g}__{c}\n{p}\n')
