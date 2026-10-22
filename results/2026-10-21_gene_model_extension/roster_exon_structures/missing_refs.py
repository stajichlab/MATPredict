import glob, yaml, os, sys, csv
db=sys.argv[1]; rows=[]
have=set(r['record'] for r in csv.DictReader(open('ref_introns.tsv'),delimiter='\t'))
for md in sorted(glob.glob(f"{db}/*/*/*/metadata.yaml")):
    if "/candidates/" in md or "/classifiers/" in md: continue
    d=yaml.safe_load(open(md))
    for g in d.get("genes") or []:
        if g["name"] in ("MAT1-1-1","MAT1-2-1") and g.get("role")=="core_MAT" and d["record_id"] not in have:
            seg=((d.get("locus") or {}).get("core") or {}).get("segments") or [{}]
            src=(seg[g.get("segment_index",0)] if len(seg)>g.get("segment_index",0) else {}).get("sequence_source") or {}
            rows.append(dict(record=d["record_id"],order=md.split("/")[-3],species=(d.get("organism") or {}).get("species"),gene=g["name"],protein=g.get("protein_accession"),
                             seq_type=src.get("type"),accession=src.get("accession"),strand=g.get("strand"),gstart=g.get("start"),gend=g.get("end"),has_exons=bool(g.get("exons")),
                             tier=((d.get("evidence") or {}).get("locus_existence") or {}).get("tier"),coord_prov=(d.get("locus") or {}).get("coordinate_provenance")))
w=csv.DictWriter(open('missing_refs.tsv','w',newline=''),rows[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(rows)
print(len(rows),'MAT1-1-1/MAT1-2-1 references without exons in the metadata')
import collections
print(collections.Counter((r['seq_type'],r['has_exons'],bool(r['protein'])) for r in rows))
for r in rows: print(r['record'][:40].ljust(40),r['gene'],r['order'][:14].ljust(14),str(r['protein'])[:24].ljust(24),r['seq_type'],r['accession'],r['coord_prov'])
