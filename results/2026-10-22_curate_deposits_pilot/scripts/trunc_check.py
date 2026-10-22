"""Is the deposit truncated? Local protein alignment of relatives against ORF segments of the deposit's nucleotide."""
import json, sys, re
from io import StringIO
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
cache = Path('ncbi_cache')
def find(prefix):
    for p in cache.iterdir():
        t = p.read_text()
        if t.startswith(prefix): return t
def prot(acc):
    t = find('>'+acc)
    return str(SeqIO.read(StringIO(t),'fasta-blast').seq) if t else None
al = PairwiseAligner(); al.mode='local'; al.substitution_matrix=substitution_matrices.load('BLOSUM62'); al.open_gap_score=-11; al.extend_gap_score=-1
def segs(nuc_acc):
    gb = SeqIO.read(StringIO(find('LOCUS       '+nuc_acc.split('.')[0])),'genbank')
    s = gb.seq
    out=[]
    for strand,seq in ((1,s),(-1,s.reverse_complement())):
        for fr in range(3):
            aa = str(seq[fr:len(seq)-((len(seq)-fr)%3)].translate())
            pos=0
            for piece in aa.split('*'):
                if len(piece)>=15:
                    a0=pos; nt0 = fr+3*a0
                    if strand==1: lo,hi = nt0+1, nt0+3*len(piece)
                    else: hi = len(seq)-nt0; lo = hi-3*len(piece)+1
                    out.append((strand,fr,lo,hi,piece))
                pos += len(piece)+1
    return gb, out
def report(nuc, queries):
    gb, ss = segs(nuc)
    print('#', nuc, len(gb.seq), 'bp')
    for qname,q in queries:
        hits=[]
        for strand,fr,lo,hi,piece in ss:
            sc = al.score(q, piece)
            if sc>=40:
                a = al.align(q, piece)[0]
                qa = a.aligned[0]; qs,qe = qa[0][0]+1, qa[-1][1]
                hits.append((qs,qe,sc,strand,lo,hi))
        for h in sorted(hits): print(f'  {qname} aa {h[0]}-{h[1]} score {h[2]:.0f} on strand {h[3]} ORF-segment {h[4]}..{h[5]}')
args = json.loads(sys.argv[1])
for nuc, qs in args.items():
    report(nuc, [(q, prot(q) or '') for q in qs])
