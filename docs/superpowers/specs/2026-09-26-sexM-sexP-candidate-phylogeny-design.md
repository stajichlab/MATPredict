# sexM/sexP candidate phylogeny for the early-diverging scan

Status: spec, not started. Written 2026-09-26 for a later agent.

## Question

The 2026-09-26 early-diverging scan found HMG-box hits that match sexM or sexP
in most Mortierellomycota (93/100) and Kickxellomycota (179/190) genomes, and
in Lichtheimiaceae, where the evidence bar holds most loci back. Best
identities are 37-46%. For each hit, we want to know if it is:

1. a sexM or sexP ortholog (a true MAT idiomorph gene), or
2. another HMG-box paralog that detection matched to sexM/sexP.

A gene tree that includes the curated sexM/sexP proteins and known non-MAT
HMG-box proteins answers this. It also tests each idiomorph label: a copy
that detection labelled Minus should group with the sexM clade.

## Inputs

- Scan results: `results/2026-09-26_early_diverging/{Mucoromycota,Mortierellomycota,Kickxellomycota}/runs/*/detection_report.yaml`.
  The `detected` and `suppressed_loci` blocks give contig, start, end,
  idiomorph and `genes_found`.
- Code and database: worktree `run-634dda4` (the code that made the scan).
- Genomes: `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes/<asmid>.fa.gz`.
- Curated references: the sexM/sexP proteins in the database (41 Mucoromycota
  proteins in each report's `_reference.faa`).
- Prior work: `testset/Zygo/sexM_sexP_hits.fa`, its alignment, the ClipKIT
  trim and the FastTree tree. Reuse these as the Mucorales backbone.

## Gap: the reports do not keep the hit sequences

`detection_report.yaml` and `evidence_diagnostics.jsonl` keep cluster spans
and identities only. They keep no protein sequence and no per-hit coordinates.
The genome FASTAs were on node-local scratch and are gone. So the first step
must re-extract the candidates.

## Steps

1. **Collect loci.** For each genome, list every locus in `detected` and in
   `suppressed_loci` where `genes_found` includes sexM or sexP. Record
   genome, taxonomy (from `samples.csv`), contig, span, idiomorph, called or
   withheld, and the genes found with it.
2. **Re-extract proteins.** Decompress the genome into `$SCRATCH`. Take the
   locus span with 5 kb padding. Align the curated sexM and sexP proteins to it
   with miniprot, and keep the best model for each HMG hit. Record the model's
   coverage of its reference. Use the genetic code from the report.
   Do not assume frame 1 or table 1 (see the genetic-code memory).
3. **Outgroups.** Add non-MAT HMG-box proteins from the same genomes, and from
   Mucorales genomes where the MAT locus is known: for example the HMG-box
   genes that the Zygo 23 regression reported as `partial_locus` (the
   sexP|sexM|glrA HMG-paralog region). Add one or two Ascomycota MAT1-2-1
   HMG proteins as a distant outgroup.
4. **Deduplicate.** Collapse identical proteins, and keep a table that maps
   each kept sequence to all the genomes it stands for. Many Coemansia
   isolates are near-clonal: RSA1085 and RSA1250 give the same locus at the
   same coordinates.
5. **Align and trim.** Use MAFFT (L-INS-i if n < 1000) and ClipKIT. Align the
   HMG box. Report the length of the trimmed alignment.
6. **Tree.** IQ-TREE with ModelFinder and 1000 UFBoot, plus SH-aLRT. Run it
   as a SLURM job on the `epyc` partition (see the nextflow-hpcc skill for
   partitions; a single job is enough).
7. **Read the tree.** For each candidate, record which clade it falls in:
   sexM, sexP, or outgroup/paralog. Compare that with the detection label
   (Plus/Minus/undetermined) and with called vs withheld.

## Outputs

`results/<date>_sexMP_phylogeny/`:
- `candidates.tsv`: one row per locus with the fields from step 1, the
  model coverage and the tree clade.
- `candidates.faa.gz`, `dedup_map.tsv`, the alignment, the trimmed alignment
  and the tree file.
- A note in `docs/notes/` that gives, per phylum/family, the counts:
  sexM clade, sexP clade, paralog, and label agreement. State the support
  values. Do not claim orthology from branches with UFBoot < 95.

## Curator rulings (2026-09-26)

- Include Lichtheimiaceae and Syncephalastraceae with the two non-Mucoromycota
  phyla: **yes**.
- Outgroup: **same-genome non-MAT HMG-box paralogs**. Add an Ascomycota
  MAT1-2-1 only if it aligns across the HMG box; report whether it did.
- **One tree**, if the sequences allow it. First measure how well the
  candidates align: pairwise identity spread, and the trimmed alignment length
  over the HMG box. Try an HMM-guided alignment (hmmalign to the Pfam HMG-box
  model, PF00505, or to a profile built from the curated sexM/sexP set) and
  compare it with MAFFT on the same input. Use the alignment that keeps more
  informative HMG-box columns. If no single alignment is usable, report why
  and fall back to one tree per phylum.

## Out of scope

- Changing detection or the evidence bar.
- The flank-ortholog synteny check (reciprocal best hits of tptA/rnhA). That
  check is separate, and the scan note lists it as a precondition. Run it in
  parallel, and use both results before the curator rules on these phyla.
