import csv,statistics,collections,os
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results/"
def load(p):
    d={}
    for r in csv.DictReader(open(p),delimiter="\t"):
        try: d[r["genome"]]=(float(r["wall_s"]),r["routing"],r["order"])
        except: pass
    return d
v8=load(R+"2026-10-31_v080_ascomycota/genomes.tsv")
for cand in ("2026-10-26_v070_ascomycota/genomes.tsv","2026-10-26_v070_ascomycota/tables/genomes.tsv"):
    if os.path.exists(R+cand): v7=load(R+cand); print("v070 from",cand,len(v7)); break
else: print([x for x in os.listdir(R+"2026-10-26_v070_ascomycota")][:20]); raise SystemExit
common=[g for g in v8 if g in v7]
print("common",len(common))
q=lambda xs,p: sorted(xs)[int(len(xs)*p)]
a=[v7[g][0] for g in common]; b=[v8[g][0] for g in common]
print("median v070 %.0f v080 %.0f  p90 %.0f %.0f  total hours %.0f %.0f"%(statistics.median(a),statistics.median(b),q(a,.9),q(b,.9),sum(a)/3600,sum(b)/3600))
ratio=[v8[g][0]/max(v7[g][0],1) for g in common]; print("median ratio %.2f p90 %.2f"%(statistics.median(ratio),q(ratio,.9)))
by=collections.defaultdict(list)
for g in common: by[v8[g][1]].append(v8[g][0]/max(v7[g][0],1))
for k,v in by.items(): print(k,len(v),"median ratio %.2f"%statistics.median(v))
