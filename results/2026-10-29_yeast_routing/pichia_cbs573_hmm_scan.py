import sys,pyhmmer
from pyhmmer.easel import TextSequence,Alphabet
D=sys.argv[1]; S=sys.argv[2]
hm={}
for n,p in [("MAT1-1-1",D+"/MAT_core/MAT1-1-1.hmm"),("MAT1-2-1",D+"/MAT_core/MAT1-2-1.hmm"),("MATa2_sacch",D+"/MATa2_saccharomycotina/MATa2_saccharomycotina.hmm")]:
    with pyhmmer.plan7.HMMFile(p) as f: hm[n]=f.read()
comp=str.maketrans("ACGT","TGCA")
b="TCAG"; aa="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
code={};i=0
for x in b:
    for y in b:
        for z in b: code[x+y+z]=aa[i]; i+=1
def tr(s): return "".join(code.get(s[i:i+3],"X") for i in range(0,len(s)-2,3))
def six(seq):
    for strand,s in (("+",seq),("-",seq.translate(comp)[::-1])):
        for f in range(3):
            p=tr(s[f:]); pos=0
            for seg in p.split("*"):
                if len(seg)>=30: yield strand,f,pos,seg
                pos+=len(seg)+1
for fn in sys.argv[3:]:
    seq="".join(l.strip() for l in open(S+"/"+fn) if not l.startswith(">")).upper()
    print("==",fn,len(seq))
    seqs={}
    orfs=[]
    for st,f,pos,s in six(seq):
        nm=f"{st}f{f}_aa{pos}"; seqs[nm]=len(s)
        orfs.append(TextSequence(name=nm.encode(),sequence=s).digitize(Alphabet.amino()))
    for n,h in hm.items():
        for hits in pyhmmer.hmmsearch(h,orfs,E=1e-1):
            for hit in hits:
                nm=hit.name if isinstance(hit.name,str) else hit.name.decode()
                print(f"  {n:12} {nm:16} orf_aa={seqs[nm]} E={hit.evalue:.1e} score={hit.score:.1f}")
