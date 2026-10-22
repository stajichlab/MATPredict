"""Independent post-check of the regenerated files: proteins.faa == NCBI protein; locus.gbk CDS (exons) translate to it."""
import sys, yaml, glob
from io import StringIO
from pathlib import Path
from Bio import SeqIO
cache = Path(sys.argv[2])
prot = {}
for p in cache.iterdir():
    t = p.read_text()
    if t.startswith('>') and not p.name.startswith('zz_'):
        r = SeqIO.read(StringIO(t), 'fasta-blast'); prot[r.id] = str(r.seq)
for meta in sorted(glob.glob(sys.argv[1] + '/db/Ascomycota/*/*/metadata.yaml')):
    m = yaml.safe_load(open(meta))
    if m['curation'].get('proposed_by') != 'deposit-to-record-generator': continue
    d = Path(meta).parent; g = m['genes'][0]; acc = g['protein_accession'].split(':')[1]
    faa = SeqIO.read(d / 'proteins.faa', 'fasta')
    gb = SeqIO.read(d / 'locus.gbk', 'genbank')
    cds = [f for f in gb.features if f.type == 'CDS']
    off = m['locus']['core']['segments'][0]['start'] - 1
    ex = [(int(p.start) + 1 + off, int(p.end) + off) for p in cds[0].location.parts]
    tr = str(cds[0].location.extract(gb.seq).translate(table=g.get('transl_table', 1)))
    print(m['record_id'], 'faa==ncbi', str(faa.seq) == prot[acc], 'gbk_cds_exons==meta', ex == [(e['start'], e['end']) for e in g['exons']],
          'gbk_translate==ncbi', tr.rstrip('*') == prot[acc], 'ends_stop', tr.endswith('*'), 'gbk_len', len(gb.seq), 'n_cds', len(cds),
          'gbk_translation_q', cds[0].qualifiers.get('translation', [''])[0] == prot[acc])
