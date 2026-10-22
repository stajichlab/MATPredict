from io import StringIO
from pathlib import Path
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
cache=Path('ncbi_cache')
def find(prefix):
    for p in cache.iterdir():
        t=p.read_text()
        if t.startswith(prefix): return t
def prot(a): return str(SeqIO.read(StringIO(find('>'+a)),'fasta-blast').seq)
gb=SeqIO.read(StringIO(find('LOCUS       MT819950')),'genbank'); s=gb.seq
al=PairwiseAligner(); al.mode='local'; al.substitution_matrix=substitution_matrices.load('BLOSUM62'); al.open_gap_score=-11; al.extend_gap_score=-1
for name,(lo,hi) in {'orfA 4848..5237':(4848,5237),'orfB 5033..6277':(5033,6277),'orfC 5865..6065':(5865,6065)}.items():
    seg=s[lo-1:hi]; aa=str(seg.translate())
    for q in ['QRL06102.1','UTN00774.1']:
        a=al.align(prot(q),aa)[0]
        (qs,qe),(ss,se)=(a.aligned[0][0][0]+1,a.aligned[0][-1][1]),(a.aligned[1][0][0],a.aligned[1][-1][1])
        print(name,q,'q',qs,qe,'of',len(prot(q)),'nt',lo+3*ss,lo+3*se-1,'score',a.score)
# what follows the deposit exon1 end (5972) in frame: the stop
print('stop codon after orfB end:', s[6277:6280], 'Drepanopeziza C-term', prot('UTN00774.1')[-15:], 'deposit C-term', prot('QRL06102.1')[-15:])
print('orfB translation tail', str(s[5033-1:6277].translate())[-20:])
# Calonectria: relative C-terminal exon
gb2=SeqIO.read(StringIO(find('LOCUS       BK013123')),'genbank'); s2=gb2.seq
for q in ['DBB25572.1','QKY89071.1','QKY89080.1']:
    print(q, prot(q)[-25:])
print('BK013123 exon3 6889..6953 translation', s2[6888:6953].translate(), ' unspliced read-through after exon2 end 6835:', s2[6835:6835+120].translate())
