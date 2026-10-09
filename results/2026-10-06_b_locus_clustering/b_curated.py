#!/usr/bin/env python3
"""Part 1/2: curated B-locus structure per record, and the receptor array that holds it."""
import os, re, glob, sys, yaml
import pandas as pd, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); DB = HERE + "/../../db/Basidiomycota"
REC = re.compile(r"receptor|^bar|^bbr|ste3", re.I); PHE = re.compile(r"pheromone|^bap|^bbp|^mf|^phb", re.I)
recs = {"Agaricales": ["5334_h4-8_Balpha_3", "5334_h4-8_Bbeta_2", "5346_a43-b43-okayama-7_PR_B43"],
        "Polyporales": ["5325_fp-101664-ss1_PR_B1", "5627_9006-11_PR_B1"], "Russulales": ["2830151_kdtol00553_PR_B1", "984962_tc-32-1_PR_B1"]}
out = []
for o, l in recs.items():
    for r in l:
        gl = [x.split("\t") for x in open(f"{DB}/{o}/{r}/locus.gff3") if not x.startswith("#")]
        g = pd.DataFrame([(x[0], int(x[3]), int(x[4]), x[6], re.search("Name=([^;]+)", x[8]).group(1)) for x in gl], columns="contig start end strand name".split())
        g["cls"] = ["R" if REC.search(n) else "P" if PHE.search(n) else "?" for n in g.name]
        m = yaml.safe_load(open(f"{DB}/{o}/{r}/metadata.yaml"))
        rc = g[g.cls == "R"]; pc = g[g.cls == "P"]
        # cassettes: receptor with >=2 curated precursors within 5 kb (edge to edge)
        cas = sum(1 for _, x in rc.iterrows() if ((pc.start - x.end).clip(lower=0).where(pc.start > x.end, (x.start - pc.end).clip(lower=0)) <= 5000).sum() >= 2)
        gaps = np.diff(np.sort(g.start.values)) if len(g) > 1 else []
        out.append(dict(record=r, order=o, species=m["organism"]["species"], contig=g.contig.iloc[0], start=g.start.min(), end=g.end.max(), span_kb=round((g.end.max() - g.start.min()) / 1000, 1),
                        n_recept=len(rc), n_pherom=len(pc), n_other=int((g.cls == "?").sum()), rec_positions=";".join(f"{a}-{b}" for a, b in zip(rc.start, rc.end)),
                        rec_gap_median_kb=round(float(np.median(np.diff(rc.start.values)))/1000, 1) if len(rc) > 1 else None,
                        cassettes_5kb_ge2prec=cas, flank_genes=len(m["locus"].get("extended_flank") or []), completeness=m["locus"]["core"]["completeness"]))
C = pd.DataFrame(out); C.to_csv(HERE + "/curated_agaricomycete_B_records.tsv", sep="\t", index=False); print(C.drop(columns="rec_positions").to_string())
# outgroup PR-type records (receptor + pheromone) for contrast
og = []
for f in glob.glob(DB + "/*/*/locus.gff3"):
    rid = f.split("/")[-2]; o = f.split("/")[-3]
    if o in recs or not re.search(r"redPR|aLocus|MAT_a|MAT_alpha|wallMAT|bLocus", rid): continue
    gl = [x.split("\t") for x in open(f) if not x.startswith("#")]
    nm = [re.search("Name=([^;]+)", x[8]).group(1) for x in gl]
    og.append(dict(order=o, record=rid, n_genes=len(nm), n_recept=sum(bool(REC.search(n)) for n in nm), n_pherom=sum(bool(PHE.search(n)) for n in nm),
                   span_kb=round((max(int(x[4]) for x in gl) - min(int(x[3]) for x in gl)) / 1000, 1), names=",".join(nm)))
O = pd.DataFrame(og); O.to_csv(HERE + "/curated_outgroup_records.tsv", sep="\t", index=False)
print(O.groupby("order").agg(n=("record", "size"), rec=("n_recept", "median"), phe=("n_pherom", "median"), span=("span_kb", "median")))
