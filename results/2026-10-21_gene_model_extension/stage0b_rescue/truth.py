"""Truth table of deposited MAT1-1-1 / MAT1-2-1 genes from the curated records. python truth.py DB_ROOT OUT_TSV"""
import sys, glob, yaml, os
db, out = sys.argv[1], sys.argv[2]
rows = []
for md in sorted(glob.glob(f"{db}/*/*/*/metadata.yaml")):
    if "/candidates/" in md or "/classifiers/" in md: continue
    d = yaml.safe_load(open(md))
    rid = d["record_id"]; order = md.split("/")[-3]
    prot = {}
    fa = os.path.join(os.path.dirname(md), "proteins.faa"); cur = None
    for l in open(fa):
        if l[0] == ">":
            cur = [x.split("=")[1] for x in l[1:].split("|") if x.startswith("name=")][0]; prot[cur] = ""
        elif cur: prot[cur] += l.strip()
    segs = {i: s for i, s in enumerate((((d.get("locus") or {}).get("core")) or {}).get("segments") or [])}
    for g in d.get("genes") or []:
        if g["name"] not in ("MAT1-1-1", "MAT1-2-1") or g.get("role") != "core_MAT" or not g.get("present", True): continue
        s = segs.get(g.get("segment_index", 0)) or {}
        src = s.get("sequence_source") or {}
        if src.get("type") != "insdc_nucleotide" or not g.get("exons"): continue
        p = prot.get(g["name"], "")
        rows.append(dict(record=rid, order=order, species=(d.get("organism") or {}).get("species"), taxid=(d.get("taxonomy") or {}).get("taxid"),
                         gene=g["name"], accession=src["accession"], strand=g["strand"], gstart=g["start"], gend=g["end"],
                         exons=";".join(f"{e['start']}-{e['end']}" for e in g["exons"]), nexons=len(g["exons"]), prot_aa=len(p), status=(d.get("validation") or {}).get("status"),
                         tier=((d.get("evidence") or {}).get("locus_existence") or {}).get("tier")))
import csv
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, rows[0].keys(), delimiter="\t"); w.writeheader(); w.writerows(rows)
print(len(rows), "genes;", {g: sum(r['gene']==g for r in rows) for g in ('MAT1-1-1','MAT1-2-1')})
