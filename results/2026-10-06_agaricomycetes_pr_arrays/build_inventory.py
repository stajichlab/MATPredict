#!/usr/bin/env python3
"""Assemble inventory/ from loci_all.tsv.gz, genome_table.tsv, BFD taxonomy and the per-genome outputs of
extract_loci_proteins.py (prot_out/*.prot.tsv, *.prot.faa).
Usage: build_inventory.py PROT_OUT_DIR TAXONOMY_TSV   (TAXONOMY_TSV: ASMID taxid phylum subphylum class order family genus species)
qpass = BUSCO complete >= 70 and N50 >= 20 kb and contigs <= 5000."""
import glob, gzip, os, sys
import numpy as np, pandas as pd

PROT, TAX = sys.argv[1], sys.argv[2]
OUT = "inventory"; os.makedirs(OUT, exist_ok=True)
tax = pd.read_csv(TAX, sep="\t", header=None, dtype=str, keep_default_na=False,
                  names=["genome", "taxid", "phylum", "subphylum", "class", "order", "family", "genus", "species"])
G = pd.read_csv("genome_table.tsv", sep="\t")
G["qpass_busco_n50"] = G.qpass
G["qpass"] = G.qpass & (G.contigs <= 5000)
G = G[G.scanned].merge(tax, on="genome", how="left", suffixes=("_g", ""))
for r in ("order", "family"):
    G[r] = G[r].replace("", np.nan).fillna(G[r + "_g"]).fillna("unclassified")
G["species_bfd"] = G.species.replace("", np.nan).fillna(G.sp)   # BFD SPECIES (may carry strain/variant text)
G["species"] = G.sp                                              # binomial used in the analysis (705 qpass species)
G["genus"] = G.genus.replace("", np.nan).fillna(G.sp.str.split().str[0])
G["class"] = G["class"].replace("", np.nan).fillna("Agaricomycetes"); G["phylum"] = G.phylum.replace("", np.nan).fillna("Basidiomycota")
TC = ["phylum", "class", "order", "family", "genus", "species", "species_bfd", "taxid"]

L = pd.read_csv("loci_all.tsv.gz", sep="\t")
L["locus_key"] = L.genome + "|" + L.contig + ":" + L.start.astype(str) + "-" + L.end.astype(str) + L.strand
P = pd.concat([pd.read_csv(f, sep="\t", dtype={"complete": str}) for f in sorted(glob.glob(f"{PROT}/*.prot.tsv"))])
assert P.locus_id.is_unique and L.locus_key.is_unique
miss = set(L.locus_key) - set(P.locus_id); extra = set(P.locus_id) - set(L.locus_key)
print("loci_all", len(L), "protein-table", len(P), "unmatched in table", len(miss), "extra", len(extra))
T = L.merge(P.drop(columns=["genome", "contig", "start", "end", "strand"]), left_on="locus_key", right_on="locus_id", how="left")
T = T.merge(G[["genome", "qpass"] + TC], on="genome", how="left")
T["caax_flag"] = T["flag"]; T["pipeline_call"] = T["in_call"]; T["withheld_cluster"] = T["in_withheld"]
T["array_id"] = T.genome + "|" + T.array.astype(str)
T["array_size"] = T.groupby("array_id").locus_id.transform("size")
T["complete"] = T.complete.map({"True": True, "False": False}).fillna(False)
T["partial_reason"] = T.partial_reason.fillna("")
T["call_status"] = np.where(T.pipeline_call, "in_pipeline_call", np.where(T.withheld_cluster, "in_withheld_cluster", "not_called"))
cols = ["locus_id", "genome", "qpass"] + TC + ["contig", "contig_len", "start", "end", "strand", "best_query", "best_qcov",
        "best_ident", "best_ref_ident", "n_queries", "model_query", "model_start", "model_end", "model_identity", "model_positive",
        "model_qcov", "n_cds", "prot_len", "complete", "partial_reason", "start_met", "has_stop", "n_internal_stop", "frameshift",
        "array_id", "array_size", "caax_flag", "hx", "d_T", "nT_10kb", "d_Hx", "nHx_10kb", "pipeline_call", "withheld_cluster", "call_status"]
T = T[cols].rename(columns={"best_ident": "locus_ident", "best_ref_ident": "locus_ref_ident", "hx": "precursor_homology_near"})
T.to_csv(f"{OUT}/ste3_loci_table.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
G[["genome", "qpass", "qpass_busco_n50"] + TC + ["size_bp", "contigs", "n50_bp", "busco_complete_pct", "status", "n_ste3"]].rename(
    columns={"n_ste3": "n_loci"}).to_csv(f"{OUT}/genomes_scanned.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0})

# FASTA
qp = set(T[T.qpass == True].locus_id)
fa, fq = gzip.GzipFile(f"{OUT}/ste3_loci_proteins.faa.gz", "wb", mtime=0), gzip.GzipFile(f"{OUT}/ste3_loci_proteins_qpass.faa.gz", "wb", mtime=0)
nf = nq = 0
for f in sorted(glob.glob(f"{PROT}/*.prot.faa")):
    txt = open(f).read()
    for rec in txt.split(">")[1:]:
        lid = rec.split(" ", 1)[0]
        fa.write((">" + rec).encode()); nf += 1
        if lid in qp:
            fq.write((">" + rec).encode()); nq += 1
fa.close(); fq.close()
print("fasta records all", nf, "qpass", nq)

# ---- distribution tables
def agg(df, rank, tag):
    g = df.groupby(rank)
    arr = df.groupby(rank).array_id.nunique()
    out = pd.DataFrame({"loci": g.size(), "arrays": arr, "genomes_with_loci": g.genome.nunique(), "species_with_loci": g.species.nunique(),
                        "complete_models": g.complete.sum(), "caax_flagged_loci": g.caax_flag.sum(),
                        "in_pipeline_call": g.pipeline_call.sum(), "in_withheld_cluster": g.withheld_cluster.sum()})
    return out.add_prefix(tag + "_")
for rank in ("class", "order", "family"):
    gs = G.groupby(rank).agg(genomes_scanned=("genome", "size"), species_scanned=("species", "nunique"))
    gq = G[G.qpass].groupby(rank).agg(qpass_genomes=("genome", "size"), qpass_species=("species", "nunique"))
    d = gs.join(gq).join(agg(T, rank, "all")).join(agg(T[T.qpass == True], rank, "qpass")).fillna(0)
    d = d.sort_values("all_loci", ascending=False)
    if rank != "class":
        d.insert(0, "order_of_family" if rank == "family" else "class", G.groupby(rank)["order" if rank == "family" else "class"].agg(lambda s: s.mode().iat[0]))
    d.astype({c: int for c in d.columns if c not in ('class', 'order_of_family')}).to_csv(f"{OUT}/tax_loci_arrays_genomes_species_by_{rank}.tsv", sep="\t")
BINS = [0, 1, 2, 3, 4, 5, 7, 11, 10**6]; LAB = ["0", "1", "2", "3", "4", "5-6", "7-10", "11+"]
for rank in ("order", "family"):
    rows = []
    for setname, S in (("qpass", G[G.qpass]), ("all_scanned", G)):
        for k, s in S.groupby(rank):
            n = s.n_ste3.astype(int); h = pd.cut(n, BINS, right=False, labels=LAB).value_counts().reindex(LAB)
            rows.append(dict(set=setname, **{rank: k}, n_genomes=len(s), n_species=s.species.nunique(), mean=round(n.mean(), 2),
                             median=n.median(), q1=n.quantile(.25), q3=n.quantile(.75), min=n.min(), max=n.max(),
                             **{f"genomes_with_{l}": int(h[l]) for l in LAB}))
    pd.DataFrame(rows).sort_values(["set", "n_genomes"], ascending=[False, False]).to_csv(f"{OUT}/per_genome_receptor_count_by_{rank}.tsv", sep="\t", index=False)
rows = []
for setname, S, TT in (("qpass", G[G.qpass], T[T.qpass == True]), ("all_scanned", G, T)):
    for k, s in S.groupby("family"):
        t = TT[TT.family == k]; st = s.status.value_counts()
        rows.append(dict(set=setname, family=k, order=s.order.mode().iat[0], n_genomes=len(s), genomes_with_loci=int((s.n_ste3 > 0).sum()),
                         **{f"genomes_status_{x}": int(st.get(x, 0)) for x in sorted(G.status.dropna().unique())},
                         loci=len(t), loci_in_pipeline_call=int(t.pipeline_call.sum()), loci_in_withheld_cluster_only=int((t.call_status == "in_withheld_cluster").sum()),
                         loci_not_called=int((t.call_status == "not_called").sum()), caax_flagged_loci=int(t.caax_flag.sum()),
                         flagged_not_called=int((t.caax_flag & ~t.pipeline_call).sum())))
pd.DataFrame(rows).sort_values(["set", "loci"], ascending=[False, False]).to_csv(f"{OUT}/call_status_by_family.tsv", sep="\t", index=False)

# reconciliation + model summary
rec = [("loci_all.tsv.gz rows (all scanned genomes with >=1 locus)", len(L), L.genome.nunique()),
       ("genomes scanned (agari_genomes.tsv scan==True)", "", len(G)),
       ("qpass BUSCO>=70 & N50>=20kb only (no contig limit)", int(L.genome.isin(G[G.qpass_busco_n50].genome).sum()), int(L[L.genome.isin(G[G.qpass_busco_n50].genome)].genome.nunique())),
       ("qpass final (+ contigs<=5000): loci", int(T.qpass.sum()), T[T.qpass == True].genome.nunique()),
       ("qpass final: genomes in scan / species", G.qpass.sum(), G[G.qpass].species.nunique()),
       ("not qpass", int((T.qpass != True).sum()), T[T.qpass != True].genome.nunique())]
pd.DataFrame(rec, columns=["set", "loci_or_count", "genomes"]).to_csv(f"{OUT}/set_reconciliation.tsv", sep="\t", index=False)
for nm, S in (("all", T), ("qpass", T[T.qpass == True])):
    print(nm, "loci", len(S), "complete", int(S.complete.sum()), "no_model", int((S.partial_reason == "no_model").sum()))
    print(S.partial_reason.value_counts().head(12).to_string())
