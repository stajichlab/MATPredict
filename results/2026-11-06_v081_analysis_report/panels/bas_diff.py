import yaml,glob,os,collections
R="/bigdata/stajichlab/jstajich/projects/MATPredict/results"; F=R+"/2026-11-06_seg_followups"
strata={l.split("\t")[0]:l.rstrip().split("\t")[2] for l in open(F+"/bas_panel.tsv")}
def rd(p):
    try: return yaml.safe_load(open(p+"/detection_report.yaml"))
    except Exception: return None
lost=collections.Counter(); gain=collections.Counter(); n=0
for g,s in strata.items():
    bp=(glob.glob(f"{R}/2026-11-05_v081_basidiomycota/runs/*/runs/{g}")+[None])[0]; a=rd(bp); b=rd(f"{F}/bas/runs/{g}")
    if not a or not b: continue
    def L(r): return {(x["family"],x["contig"],x["start"]//20000):x for x in (r.get("detected") or [])}
    la,lb=L(a),L(b)
    out=[]
    for k in set(la)|set(lb):
        if k in la and k in lb:
            ga={e["gene"]:(e.get("status"),e.get("method"),e.get("bitscore")) for e in la[k]["gene_evidence"] if not e.get("lost_to")}
            gb={e["gene"]:(e.get("status"),e.get("method"),e.get("bitscore")) for e in lb[k]["gene_evidence"] if not e.get("lost_to")}
            for x in set(ga)-set(gb): lost[x]+=1; out.append(("lost",x,ga[x]))
            for x in set(gb)-set(ga): gain[x]+=1; out.append(("gained",x,gb[x]))
            if la[k]["idiomorph"]!=lb[k]["idiomorph"] or la[k]["confidence"]!=lb[k]["confidence"]: out.append(("tier/idio",(la[k]["idiomorph"],la[k]["confidence"]),(lb[k]["idiomorph"],lb[k]["confidence"])))
        elif k in la: out.append(("LOCUS LOST",k[0],la[k]["idiomorph"],la[k]["confidence"],sorted(e["gene"] for e in la[k]["gene_evidence"] if not e.get("lost_to"))))
        else: out.append(("LOCUS GAINED",k[0],lb[k]["idiomorph"],lb[k]["confidence"],sorted(e["gene"] for e in lb[k]["gene_evidence"] if not e.get("lost_to"))))
    if out: print(g[:36],s,out[:5])
print("genes lost from loci:",lost.most_common(10)); print("genes gained:",gain.most_common(10))
