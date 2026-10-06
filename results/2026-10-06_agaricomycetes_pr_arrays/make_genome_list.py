#!/usr/bin/env python3
"""Agaricomycetes genome list of the v0.6.0 Basidiomycota run + quality filter.
Quality pass = BUSCO complete >= 70 and N50 >= 20 kb (BFD asm_stats / busco_genome)."""
import pandas as pd
V = "../2026-10-03_basidiomycota_v060/"
g = pd.read_csv(V + "genomes.tsv", sep="\t")
g = g[g.class_ == "Agaricomycetes"].copy()
q = pd.read_csv("bfd_quality.tsv", sep="\t")
g = g.merge(q, on="genome", how="left")
g["qpass"] = (g.busco_complete_pct >= 70) & (g.n50_bp >= 20000)
g["scan"] = g.size_bp.fillna(0) < 500_000_000
g[["genome", "species", "order", "family", "size_bp", "contigs", "n50", "busco_complete_pct", "n50_bp", "qpass", "scan",
   "status", "n_loci", "families_called"]].to_csv("agari_genomes.tsv", sep="\t", index=False)
print(len(g), "genomes;", int(g.qpass.sum()), "pass quality;", int((~g.scan).sum()), ">=500 Mb;",
      int(g.busco_complete_pct.isna().sum()), "no BUSCO")
print(g.groupby("order").agg(n=("genome", "size"), qpass=("qpass", "sum")).sort_values("n", ascending=False))
