"""Translate GAP1/BGL1 CDS joins from the C4/C5 deposits and compare with the NCBI proteins."""
import sys
from Bio import SeqIO
gb_dir, faa = sys.argv[1], sys.argv[2]
prot = {r.id: str(r.seq) for r in SeqIO.parse(faa, 'fasta')}
for acc in ('AF027687.1', 'AF029913.1'):
    rec = SeqIO.read(f'{gb_dir}/{acc}.gb', 'genbank')
    for f in rec.features:
        if f.type != 'CDS': continue
        g = f.qualifiers.get('gene', ['?'])[0]; pid = f.qualifiers['protein_id'][0]
        tr = str(f.location.extract(rec.seq).translate(table=1)).rstrip('*')
        exons = [(int(p.start)+1, int(p.end), p.strand) for p in f.location.parts]
        print(acc, g, pid, len(tr), 'internal_stop' if '*' in tr else 'ok', 'match' if prot.get(pid) == tr else ('n/a' if pid not in prot else 'MISMATCH'), exons, f.qualifiers.get('translation',[''])[0]==tr)
