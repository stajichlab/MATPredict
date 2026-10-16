import csv, glob, os, sys, collections, yaml
L = getattr(yaml, "CSafeLoader", yaml.SafeLoader)
src, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
meta = {}
for r in csv.DictReader(open("/Users/jstajich/projects/MATPredict/results/2026-10-03_ascomycota_v060/loci.tsv"), delimiter="\t"):
    meta[(r["genome"], r["contig"], int(r["start"]))] = r
TA, TH = 31.6, 35.3          # 5th percentile identity of single-gene loci (alpha, HMG)
def modelled(g): return g["status"] in ("polished_agree", "polished_single", "polished_disagree")
rows = []
for path in glob.glob(src + "/**/detection_report.yaml", recursive=True):
    genome = os.path.basename(os.path.dirname(path))
    doc = yaml.load(open(path), Loader=L)
    for c in doc.get("detected") or []:
        if c.get("family") != "Ascomycota:MAT": continue
        ev = {}
        for g in c.get("gene_evidence") or []:
            if g["gene"] in ("MAT1-1-1", "MAT1-2-1") and not g.get("lost_to"):
                if g["gene"] not in ev or (g.get("bitscore") or 0) > (ev[g["gene"]].get("bitscore") or 0):
                    ev[g["gene"]] = g
        if len(ev) < 2: continue
        a, h = ev["MAT1-1-1"], ev["MAT1-2-1"]
        gap = (a.get("bitscore") or 0) - (h.get("bitscore") or 0)
        higher = "MAT1-1" if gap > 0 else "MAT1-2"
        ma, mh = modelled(a), modelled(h)
        groups = []
        if not ma and not mh and c["idiomorph"] != higher: groups.append("neither_modelled_label_disagrees")
        if ma and mh and abs(gap) < 20: groups.append("both_modelled_gap_lt20")
        if a["identity"] < TA and h["identity"] < TH: groups.append("both_weak")
        m = meta.get((genome, c["contig"], int(c["start"])), {})
        rows.append({
            "genome": genome, "class": m.get("class_", ""), "order": m.get("order", ""), "contig": c["contig"],
            "start": c["start"], "end": c["end"], "label": c["idiomorph"], "confidence": c["confidence"],
            "locus_class": c["locus_class"], "higher_bitscore_gene": higher, "label_disagrees": c["idiomorph"] != higher,
            "bitscore_gap": round(gap, 1),
            "alpha_identity": round(a["identity"], 1), "alpha_bitscore": round(a.get("bitscore") or 0, 1),
            "alpha_status": a["status"], "alpha_method": a["method"],
            "hmg_identity": round(h["identity"], 1), "hmg_bitscore": round(h.get("bitscore") or 0, 1),
            "hmg_status": h["status"], "hmg_method": h["method"],
            "genes_found": "|".join(c.get("genes_found") or []), "groups": ",".join(groups)})
def write(name, sel):
    sel = sorted(sel, key=lambda r: (r["class"], r["genome"], r["contig"], r["start"]))
    with open(os.path.join(out, name), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t"); w.writeheader(); w.writerows(sel)
    return len(sel)
print("ambiguous_loci.tsv", write("ambiguous_loci.tsv", [r for r in rows if r["groups"]]))
dis = [r for r in rows if r["label_disagrees"]]
print("label_disagrees.tsv", write("label_disagrees.tsv", dis))
c = collections.Counter(r["groups"] or "(none)" for r in rows if r["groups"] or r["label_disagrees"])
print(c.most_common())
print("label disagrees by class:", collections.Counter(r["class"] for r in dis).most_common())
print("label disagrees by label:", collections.Counter((r["label"], r["higher_bitscore_gene"]) for r in dis))
print("label disagrees by modelled:", collections.Counter((r["alpha_status"] != "unpolished" and "polished" in r["alpha_status"], "polished" in r["hmg_status"]) for r in dis))
print("label disagrees by confidence:", collections.Counter(r["confidence"] for r in dis), collections.Counter(r["locus_class"] for r in dis))
print("gap quantiles", sorted(abs(r["bitscore_gap"]) for r in dis)[::20])
