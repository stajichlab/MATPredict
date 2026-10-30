import sys,collections
from Bio.Seq import Seq
S=sys.argv[1]
seq="".join(l.strip() for l in open(S+"/win.fa") if not l.startswith(">")).upper()
OFF=294999
g=lambda a,b: seq[a-OFF-1:b-OFF]          # 1-based inclusive genome coords
dep={}
for l in open(S+"/depth.tsv"):
    c,p,d=l.split(); dep[int(p)]=int(d)
print("depth per 50bp 310500-313500")
for a in range(310500,313500,50):
    print(a, round(sum(dep.get(i,0) for i in range(a,a+50))/50))
I1=(311842,311903)
# spliced CDS candidates: exon1 ends at 311841; try each frame, walk back to find ATG upstream with no stop
best=[]
for st in range(310800,311841):
    if g(st,st+2)!="ATG": continue
    cds=g(st,311841)+g(311904,313400)
    p=str(Seq(cds[:len(cds)//3*3]).translate())
    stop=p.find("*")
    # require no stop within exon1 part
    e1=(311841-st+1)
    if stop*3<e1: continue
    best.append((stop,st,p[:stop]))
best.sort(reverse=True)
for stop,st,p in best[:4]:
    end=311903+(stop*3+3-(311841-st+1))
    print("start",st,"aa",stop,"stop codon ends",end); print(p)
open(S+"/mata2_candidates.faa","w").write("".join(f">cand_{st}_{311841}_{311904}\n{p}\n" for stop,st,p in best[:4]))
