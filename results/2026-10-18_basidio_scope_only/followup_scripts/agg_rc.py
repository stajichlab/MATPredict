import json, csv, glob, collections, os
W = "/bigdata/stajichlab/jstajich/projects/mat-fam18/"
grp = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(W + "rc/groups.tsv")}
rows = []
for f in sorted(glob.glob(W + "rc/out/*.json")):
    d = json.load(open(f)); acc = d["acc"]
    gen, status, kind = grp[acc][1:4]
    mb = d["miniprot_best"]
    def top(g):
        v = mb.get(g) or []
        return v[0] if v else None
    e, w = top("bE"), top("bW")
    tp = d["tblastn_pairs"]; mpp = d["miniprot_pairs"]
    full = [p for p in mpp if p["bE_qcov"] >= 0.8 and p["bW_qcov"] >= 0.8]
    if full: cat = "A_pair_complete_models"
    elif mpp or tp: cat = "B_pair_on_one_contig_fragmentary"
    elif (e or d["tblastn"]["bE"]) and (w or d["tblastn"]["bW"]): cat = "C_both_genes_hit_not_on_one_contig"
    elif e or w or d["tblastn"]["bE"] or d["tblastn"]["bW"]: cat = "D_one_gene_only"
    else: cat = "E_no_hit"
    best_t = tp[0] if tp else None
    hd = d["pfam_homeodomain"]
    rows.append(dict(acc=acc, genus=gen, status=status, group=kind, genome_mb=round(d["genome_bp"] / 1e6, 1), contigs=d["n_contigs"], n50=d["n50"],
                     category=cat,
                     mp_bE=f"{e['ident']:.2f}/{e['qcov']}/fs{e['frameshift']}/{e['bp']}bp" if e else "", mp_bW=f"{w['ident']:.2f}/{w['qcov']}/fs{w['frameshift']}/{w['bp']}bp" if w else "",
                     mp_pair=(f"gap{mpp[0]['gap']}/{mpp[0]['orient']}/clen{mpp[0]['contig_len']}" if mpp else ""),
                     tb_bE=(f"{d['tblastn']['bE'][0]['best_bits']}/{d['tblastn']['bE'][0]['qcov']}" if d["tblastn"]["bE"] else ""),
                     tb_bW=(f"{d['tblastn']['bW'][0]['best_bits']}/{d['tblastn']['bW'][0]['qcov']}" if d["tblastn"]["bW"] else ""),
                     tb_pair=(f"gap{best_t['gap']}/{best_t['orient']}/clen{best_t['contig_len']}" if best_t else ""),
                     pfam_domains=len({(x['contig'], x['domain']) for x in hd}), pfam_contigs=len({x['contig'] for x in hd})))
with open(W + "rc/summary.tsv", "w") as fh:
    w_ = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w_.writeheader(); w_.writerows(rows)
print(len(rows), "genomes")
for k in ("target", "uncalled_other", "control"):
    sel = [r for r in rows if r["group"] == k]
    print(k, len(sel), dict(collections.Counter(r["category"] for r in sel)))
print("by genus/status (category):")
c = collections.defaultdict(collections.Counter)
for r in rows: c[(r["genus"], r["status"])][r["category"][0]] += 1
for k, v in sorted(c.items()): print("  ", k, dict(v))
n50 = collections.defaultdict(list)
for r in rows: n50[(r["status"])].append(r["n50"])
import statistics
for k, v in n50.items(): print("N50", k, "n", len(v), "median", statistics.median(v), "min", min(v), "max", max(v))
for g in ("Hemileia", "Endocronartium", "Uromyces", "Melampsora", "Puccinia", "Phakopsora", "Cronartium", "Austropuccinia"):
    sel = [r for r in rows if r["genus"] == g]
    if sel: print(g, "n", len(sel), "N50 median", statistics.median(r["n50"] for r in sel), "genome Mb median", statistics.median(r["genome_mb"] for r in sel), "contigs median", statistics.median(r["contigs"] for r in sel))
