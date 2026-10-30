import csv,sys,collections,re
sys.path.insert(0,"/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit")
exec(open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-26_v070_submit/taph_synteny1.py").read().split("os.makedirs(O+\"/dm\",exist_ok=True)")[0])
orth=collections.defaultdict(dict)       # label -> target protein id -> (pjir start, side, pident)
for r in csv.DictReader(open(O+"/pjir_flank_orthologs.tsv"),delimiter="\t"): orth[r["target"]][r["protein"]]=(int(r["pjir_start"])//1000,r["side"],float(r["pident"]))
cand=[]
for r in csv.DictReader(open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/loose/candidates.tsv"),delimiter="\t"): cand.append((r["contig"],int(r["start"]),int(r["end"]),r["ref_gene"],r["identity"],r["species"]))
dom=[]
for r in csv.DictReader(open("/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/domain/domain_hits.tsv"),delimiter="\t"): dom.append((r["contig"],int(r["start"]),int(r["end"]),r["ref_gene"],r["bitscore"]))
def show(label,contig,s,e):
    G=gff_genes(f"{A}/{label}.genomic.gff"); prot,desc=fasta(f"{A}/{label}.protein.faa")
    sel=sorted([g for g in G.values() if g["contig"]==contig and g["end"]>=s and g["start"]<=e],key=lambda g:g["start"])
    print(f"\n== {label} {contig}:{s}-{e}: {len(sel)} genes")
    for g in sel:
        o=orth[label].get(g["pid"]); ot=f"<- Pjir {o[0]}kb ({o[1]}) {o[2]:.0f}%" if o else ""
        c=[f"MAT-cand {x[3]} id{x[4]}" for x in cand if x[0]==contig and min(x[2],g["end"])-max(x[1],g["start"])>0]
        d=[f"domain {x[3]} {x[4]}b" for x in dom if x[0]==contig and min(x[2],g["end"])-max(x[1],g["start"])>0]
        dd=re.sub(r"^\S+\s+","",desc.get(g["pid"],""))[:44]
        print(f"  {g['start']:>8}-{g['end']:<8}{g['strand']} {dd:<44} {ot} {' '.join(c[:1]+d[:1])}")
    # candidates / domain hits not in any annotated gene
    for x in cand+dom:
        if x[0]==contig and x[2]>=s and x[1]<=e and not any(min(x[2],g["end"])-max(x[1],g["start"])>0 for g in sel): print("   (outside annotated genes)",x)
show("Plactucae","MCFI01000029.1",20000,70000)
