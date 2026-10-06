import csv, collections, statistics as st, sys
meta = {r["ASMID"]: r for r in csv.DictReader(open("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/samples.csv"))}
rows = list(csv.DictReader(open("besthits.tsv"), delimiter="\t"))
by = collections.defaultdict(dict)
for r in rows:
    g, (src, gene) = r["genome"], r["query"].split("|")
    r["ev"] = float(r["evalue"]); r["nx"] = float(r["next_evalue"]) if r["next_evalue"] else None
    key = (gene, src)
    by[g][key] = r
out = []
def fmt(e): return f"{e:.0e}" if e is not None else "-"
print(f"{'genome':34s} {'family':15s} {'SLA2 Mimp E / %id / next / RBH':38s} {'APN2 Mimp E / %id / next / RBH':38s} {'MAT gene':22s} {'SLA2-MAT':>9s} {'APN2-MAT':>9s}")
stats = collections.Counter(); sla_ev=[]; apn_ev=[]; sla_ev_afu=[]; apn_ev_afu=[]; dist_s=[]; dist_a=[]
for g in sorted(by):
    d = by[g]; fam = meta.get(g, {}).get("FAMILY", "?")
    def cell(gene):
        r = d.get((gene, "Mimp"))
        if not r: return "none", None
        return f"{fmt(r['ev'])} / {float(r['pident']):.0f} / {fmt(r['nx'])} / {'Y' if r['rbh_afu']=='True' else 'n'}", r
    cs, rs = cell("SLA2"); ca, ra = cell("APN2")
    mats = [d[k] for k in d if k[0] in ("MAT1-1-1", "MAT1-2-1") and d[k]["ev"] < 1e-5]
    mat = min(mats, key=lambda r: r["ev"]) if mats else None
    def dist(r):
        if not (r and mat and r["contig"] and r["contig"] == mat["contig"]): return None
        a1, b1, a2, b2 = int(r["start"]), int(r["end"]), int(mat["start"]), int(mat["end"])
        return max(0, max(a1, a2) - min(b1, b2))
    ds, da = dist(rs), dist(ra)
    mcell = f"{mat['query'].split('|')[1]} {fmt(mat['ev'])}" if mat else "none <1e-5"
    print(f"{g:34s} {fam[:15]:15s} {cs:38s} {ca:38s} {mcell:22s} {('' if ds is None else ds):>9} {('' if da is None else da):>9}")
    stats["genomes"] += 1
    if rs: stats["sla2_found"] += 1; sla_ev.append(rs["ev"]); stats["sla2_rbh"] += rs["rbh_afu"] == "True"
    if ra: stats["apn2_found"] += 1; apn_ev.append(ra["ev"]); stats["apn2_rbh"] += ra["rbh_afu"] == "True"
    for gene, lst in (("SLA2", sla_ev_afu), ("APN2", apn_ev_afu)):
        r = d.get((gene, "Afum"));  lst.append(r["ev"]) if r else None
    if mat: stats["mat_found"] += 1
    if ds is not None: dist_s.append(ds); stats["sla2_same_contig_as_mat"] += 1
    if da is not None: dist_a.append(da); stats["apn2_same_contig_as_mat"] += 1
print("\nSUMMARY", dict(stats))
for lab, v in (("SLA2 Mimp", sla_ev), ("APN2 Mimp", apn_ev), ("SLA2 Afum", sla_ev_afu), ("APN2 Afum", apn_ev_afu)):
    if v: print(f"  {lab}: n={len(v)} median E={st.median(v):.0e} worst E={max(v):.0e}")
for lab, v in (("SLA2-MAT gap bp", dist_s), ("APN2-MAT gap bp", dist_a)):
    if v: print(f"  {lab}: n={len(v)} median={st.median(v):.0f} max={max(v)} <=20kb={sum(x<=20000 for x in v)}")
