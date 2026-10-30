import pyhmmer,sys
from pyhmmer.easel import TextSequence,Alphabet,SequenceFile
H="/Users/jstajich/projects/MATPredict/.claude/worktrees/yeast-routing/db/Ascomycota/classifiers/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm"
A=Alphabet.amino()
seqs=[]
with SequenceFile(sys.argv[1],digital=True,alphabet=A) as f: seqs=list(f)
with pyhmmer.plan7.HMMFile(H) as hf: hmm=hf.read()
for hits in pyhmmer.hmmsearch(hmm,seqs,E=10):
    for h in hits:
        d=h.best_domain
        print(h.name,"score",round(h.score,1),"E",h.evalue,"ali",d.alignment.target_from,d.alignment.target_to,"len",len(h.name))
