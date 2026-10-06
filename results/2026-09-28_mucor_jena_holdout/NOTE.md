# Mucor_Jena held-out test (2026-09-28)

Blind leave-out set: /bigdata/stajichlab/shared/projects/ZyGoLife/Mucor_Jena/annotation.
Nothing from it entered training, curation or the classifier. Keyed by strain
folder only; species names in file names were not used.

Code: frozen worktree .claude/worktrees/run-076afe4 (PR #9 head 076afe4).
Run: `matpredict detect --genome <file> --phylum Mucoromycota`, no taxid.
Routing was `explicit_phylum` for all 64 genomes. Genome phylum is unknown,
so no call carries the out-of-phylum `unverified` label (verification null in
all 66 calls).

Inputs: 65 folders. CBS169_25 has no assembly (dropped). CBS334_71 has only
predict_results (scaffolds used; no contigs file). 64 genomes run on
scaffolds, 63 on contigs.

## Leakage
- CBS293_63: its BFD copy GCA_052058895.1 supplied a sexP protein to the
  classifier training_extra set. Exclude from scoring.
- CBS210_80: its BFD copy GCA_060309335.1 is in BFD and was one of the seven
  label-change cases used to set min_margin 25. Treat as not independent.
- No strain matches a curated db record, Zygo 23, or other training entries.

## Results (see summary.txt)
- Scaffolds: 61/64 called; 66 calls; 5 genomes with 2 calls.
- Calls: Minus 36 (29 high, 7 medium), Plus 29 (25 high, 4 medium),
  undetermined 1 (low). Classifier input: model 65, hsp_fragment 1.
- MAT-gene gate withheld 3 loci. Split-locus fired 0 times.
- Contigs: identical calls to scaffolds for all 63.
- Not called: CBS186_87, CBS538_80, CBS564_66.

## Gene models vs funannotate (gene_models.tsv)
- 65 rows for the called idiomorph's gene: 15 not annotated at the locus;
  49 match the annotated protein at >=95% identity; 6 annotations are
  shorter than 70% of detect's model.
- Classifier on the annotated protein agrees with detect in 98 of 99
  comparable rows (1 disagreement: CBS210_80 second call).

## Flags for the curator (annotated_sexP_vs_call.tsv)
- Five strains called Minus carry a strong annotated sexP-like protein
  (219-307 bits) on another scaffold: CBS336_62, CBS421_70, CBS608_78,
  CBS763_74 (sexP 300+), and weaker CBS156_58/CBS230_35/CBS893_73 (~87-90).
  Detect drops most of those loci at the fraction floor. Candidate
  homothallic strains or sexP paralogs; needs truth/taxonomy.
- Two strains have a Plus (margin ~250) and a weak Minus (margin 30-34) call
  on different scaffolds: CBS169_57, CBS251_35.
- CBS538_80 is uncalled but carries a 290-bit sexP-like annotated protein on
  a short scaffold (split locus that the rule did not recover).
