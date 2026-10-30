import re,collections,sys
A="/bigdata/stajichlab/jstajich/projects/MATPredict/results/2026-10-28_taphrinomycotina_mat_candidates/synteny/annot"
def genes(gff):
    out=[]
    for l in open(gff):
        if l.startswith("#"): continue
        f=l.rstrip("\n").split("\t")
        if len(f)<9 or f[2] not in ("gene",): continue
        at=dict(kv.split("=",1) for kv in f[8].split(";") if "=" in kv)
        out.append(dict(contig=f[0],start=int(f[3]),end=int(f[4]),strand=f[6],id=at.get("ID",""),name=at.get("Name",""),gene=at.get("gene",""),biotype=at.get("gene_biotype",""),desc=at.get("description") or at.get("product","")))
    return out
def region(label,contig,s,e,pad=70000):
    G=genes(f"{A}/{label}.genomic.gff"); sel=[g for g in G if g["contig"]==contig and g["end"]>=s-pad and g["start"]<=e+pad]
    sel.sort(key=lambda g:g["start"]); print(f"== {label} {contig}:{s}-{e} (+/-{pad//1000} kb): {len(sel)} genes")
    for g in sel:
        tag="  <== PM" if g["end"]>=s and g["start"]<=e else ""
        print(f"  {g['start']:>8}-{g['end']:<8} {g['strand']} {g['name'] or g['gene']:<14} {g['biotype'][:6]:6} {g['desc'][:70]}{tag}")
region("Pjir_RU7","NW_017264783.1",211240,221394,60000)
