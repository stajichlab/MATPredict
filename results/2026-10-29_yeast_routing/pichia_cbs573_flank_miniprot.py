import gzip,subprocess,re
W="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-29_yeast_routing/pichia"
L="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/GCF_003054445.1_ASM305444v1.fa.gz"
MP="/bigdata/stajichlab/jstajich/projects/MATPredict/.pixi/envs/default/bin/miniprot"
A,B=295000,325000
n=None;s=[]
for l in gzip.open(L,"rt"):
    if l.startswith(">"):
        if n=="NC_042507.1": break
        n=l[1:].split()[0];s=[]
    else:
        if n=="NC_042507.1": s.append(l.strip())
seq="".join(s); print("contig",n,len(seq))
open(W+"/win.fa","w").write(">NC_042507.1_%d_%d\n%s\n"%(A,B,seq[A-1:B]))
out=subprocess.run([MP,"--gff","-t2",W+"/win.fa",W+"/ref_pk.faa"],capture_output=True,text=True)
best={}
for ln in out.stdout.splitlines():
    f=ln.split("\t")
    if len(f)>8 and f[2]=="mRNA":
        at=dict(x.split("=",1) for x in f[8].split(";") if "=" in x)
        t=at["Target"].split()[0]
        rec=(float(at.get("Identity",0))*100,int(f[3])+A-1,int(f[4])+A-1,f[6],int(at.get("Rank",1)),at.get("Positive",""))
        best.setdefault(t,[]).append(rec)
for t,v in sorted(best.items(), key=lambda x:min(r[1] for r in x[1])):
    for r in sorted(v,key=lambda r:r[1])[:2]: print(f"{t:55} id={r[0]:.0f}% {r[1]}-{r[2]} {r[3]} rank{r[4]}")
print(out.stderr[-300:])
