import sys; sys.path.insert(0,'.')
import modelkit as mk, yaml
from pathlib import Path
g=sys.argv[1]
d=yaml.load(open(f"runs/runs190r/{g}/detection_report.yaml"),Loader=mk.L)
for l in d["detected"]:
  for e in l["gene_evidence"]:
    if e["gene"]=="MAT1-1-1" and e.get("exons"):
        sq=mk.contigs(Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes")/f"{g}.fa.gz",{e["contig"]})[e["contig"]]
        m=mk.gene_models(e,sq,1)
        for k,v in m.items():
            print(k,len(v["prot"]),"stops at",[i for i,c in enumerate(v["prot"]) if c=="*"],"startM",v["starts_M"],"up_to_stop",v["up_to_stop"],"up_M",v["up_M"],"ext",v["ext_aa"])
            print(v["prot"][:40],"...",v["prot"][-30:], "EXT:",v["ext_seq"][:30])
        print("exons",e["strand"],[(a["start"],a["end"]) for a in e["exons"]],"alt",[(a["start"],a["end"]) for a in (e.get("alternate_model") or {"exons":[]})["exons"]])
