import os
D=os.path.expanduser("~/.taxonkit")
nodes={}; 
for l in open(f"{D}/nodes.dmp"):
    p=[x.strip() for x in l.split("|")]; nodes[p[0]]=(p[1],p[2])
sci={}; byname={}
for l in open(f"{D}/names.dmp"):
    p=[x.strip() for x in l.split("|")]
    if p[3]=="scientific name": sci[p[0]]=p[1]; byname.setdefault(p[1],[]).append(p[0])
gen=[x.strip() for x in open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/genera.txt") if x.strip()]
def lin(t):
    out={}
    while True:
        par,rank=nodes[t]
        out[rank]=sci[t]
        if par==t: break
        t=par
    return out
rows=[]
for g in gen:
    c=[t for t in byname.get(g,[]) if nodes[t][1]=="genus" and lin(t).get("kingdom")=="Fungi"]
    if not c: rows.append((g,"","","","","","NOT FOUND")); continue
    L=lin(c[0]); rows.append((g,c[0],L.get("phylum",""),L.get("class",""),L.get("order",""),L.get("family",""),"" if len(c)==1 else "multiple"))
open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/genus_lineage.tsv","w").write("genus\ttaxid\tphylum\tclass\torder\tfamily\tnote\n"+"".join("\t".join(r)+"\n" for r in rows))
for r in rows: print("\t".join(r))
