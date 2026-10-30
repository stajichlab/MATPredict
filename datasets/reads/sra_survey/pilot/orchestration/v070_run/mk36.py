import glob
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"
D=f"{R}/2026-10-26_v070_submit/compare"
tax={}
for f in glob.glob(f"{R}/2026-10-26_v070_ascomycota/lists/wave_*.tsv")+glob.glob(f"{R}/2026-10-26_v070_basidiomycota/lists/*.tsv"):
    for l in open(f):
        p=l.rstrip("\n").split("\t")
        if len(p)>=2: tax[p[0]]=p[1]
def ids(*names):
    out=[]
    for n in names: out+=open(f"{D}/{n}").read().split()
    return out
a5=ids("diff_asco_model.txt"); b31=ids("diff_basidio_call.txt","diff_basidio_model.txt")
A=f"{R}/2026-10-26_v070_submit/attr"
import os; os.makedirs(A,exist_ok=True)
for name,g in (("list_asco5",a5),("list_basidio31",b31),("list_all36",a5+b31)):
    open(f"{A}/{name}.tsv","w").write("".join(f"{x}\t{tax[x]}\n" for x in g))
    print(name,len(g))
