import re
O="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates"
for l in open(O+"/references.faa"):
    if l.startswith(">"): pass
n=0;cur=None;ln={}
for l in open(O+"/references.faa"):
    if l.startswith(">"): cur=l[1:].strip(); ln[cur]=0
    else: ln[cur]+=len(l.strip())
for k,v in ln.items():
    if k.endswith("|PM"): print(k,v)
