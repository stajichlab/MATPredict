# Non-MAT HMG genes in the sexP/sexM gene trees (2026-10-04)

Report: analysis/2026-10-04_unclassified-hmg-in-gene-trees.md.

- `place_hmg.py` -> `hmg_outgroup_placement.tsv`: for each outgroup tip in the full-length
  and HMG-box trees, the first clade with MAT tips (types, UFBoot, sizes), nearest MAT
  tip, genome, contig/coords from results/2026-09-27_sexMP_fasttree/candidates.tsv (blank
  for non-locus copies there), the genome's v0.6.0 calls.
- `sexP_sister_block.txt`: walk up from the sexP MRCA and the sister groups added at each
  node (orders, the genomes' MAT calls).
- `GCF_025528875.1_Mycafr1_h1.faa`, `GCA_025716815.1_Dicele1_h5.faa`: the two genes checked
  by tblastn against their genomes (results in the report).
