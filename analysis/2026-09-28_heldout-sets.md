# Held-out test sets: Zygo 23, Mucor_Jena, ZyGoLife LCG
Status: decided (runs complete); open (curator truth tables for Jena and LCG)

## Question
How does detection perform on genomes that never entered training or curation,
and are the MAT proteins called properly compared with the funannotate
annotation?

## Data and code version
- Code: frozen worktree `run-076afe4` (PR #9 at 076afe4); `detect --phylum
  Mucoromycota`, no taxid.
- Zygo 23: 23 Mucorales LCG genomes (16 Plus, 7 Minus; truth from the curator's
  earlier BLAST-based calls in `testset/Zygo/`).
- Mucor_Jena: `/bigdata/stajichlab/shared/projects/ZyGoLife/Mucor_Jena/annotation`,
  65 strains (newer funannotate); keyed by strain folder only.
- LCG: `/bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/`, 897
  genomes (older funannotate); 621 Mucoromycotina scored, the rest discovery-only.
- Neither set may enter training or curation (curator ruling 2026-09-28).

## Method
1. Leakage check against curated records, classifier training, BFD, Zygo 23
   and each other.
2. Detection on scaffolds (Jena also on contigs).
3. Compare detect's gene models with the annotated proteins; run the classifier
   on the annotated proteins.

## Results

### Zygo 23
Saturated: 23/23 locus and mating type under every rule change so far (Fable
review). It no longer discriminates between rules.

### Mucor_Jena
Source: `results/2026-09-28_mucor_jena_holdout/NOTE.md`.
- Leakage: CBS293_63 excluded (its BFD copy supplied a sexP protein to
  training); CBS210_80 not independent (one of the 7 cases used to set the
  25-bit margin floor). No other overlap.
- 61/64 called; 66 calls; Minus 36 (29 high), Plus 29 (25 high), undetermined 1.
  Contig input gave identical calls. Uncalled: CBS186_87, CBS538_80, CBS564_66.
- Gene models vs annotation: 15 of 65 not annotated at the locus; 49 match at
  >= 95% identity; classifier agrees on 98 of 99 comparable rows.
- Flags: 7 Minus-called strains carry a sexP-like annotated protein elsewhere
  (CBS336_62, CBS421_70, CBS608_78, CBS763_74 at 219-307 bits; CBS156_58,
  CBS230_35, CBS893_73 ~87-90). With the P1 class (see classifier-builds report)
  CBS763_74 becomes Plus (309.6), and CBS221_71 / CBS223_63 are called Minus.

### ZyGoLife LCG
Source: `results/2026-09-28_lcg_holdout/NOTE.md`.
- Leakage tiers: clean 533; identical annotated MAT protein to training 29;
  training segment >= 60 aa 22; in BFD but not training 8; training strains 6;
  Zygo 23. No overlap with Jena, so old vs new funannotate was not compared on
  the same strain.
- 536/621 called (86.3%); clean tier 449 called: Minus 224, Plus 190, both 28,
  undetermined 7. Zygo subset 23/23.
- 36 genomes with Plus and Minus on different contigs (see two-idiomorphs
  report).
- 85 uncalled: 34 Lichtheimiaceae (known gap); A. blakesleeana 0/3 (known gap,
  see paralogs report); 62 of 85 carry an annotated protein >= 100 bits to
  sexM/sexP.
- Old annotation vs models: 322 of 576 (56%) of detect's MAT models are not
  annotated at the locus; classifier agrees 247/247 where both exist.

## What changed in detection
None (evaluation only). Later changes measured on these sets: P1 class,
gate threshold per build.

## Limits
- No known mating types yet for most Jena/LCG genomes; no accuracy claimed.
- LCG leakage by identical sequence covers annotated proteins only (undercounts).
- 11 LCG genomes have no collection number, so the strain check misses them.
- The Jena file names carry species names, so the test is not fully blind.

## Curator decisions
- Made 2026-09-28: both sets are held out; results keyed by strain.
- Open: curator taxonomy and mating types for
  `results/2026-09-28_mucor_jena_holdout/per_strain.tsv` and
  `results/2026-09-28_lcg_holdout/curator_table.tsv`.

## Files
`results/2026-09-28_mucor_jena_holdout/` (per_strain.tsv, gene_models.tsv,
annotated_protein_classifier.tsv, annotated_sexP_vs_call.tsv, leakage.tsv);
`results/2026-09-28_lcg_holdout/` (per_genome.tsv, curator_table.tsv,
leakage.tsv, model_vs_annotation.tsv).

## Update 2026-10-01
Both sets were rerun on PR #9 `52b3ff9` with curator species names applied:
see 2026-10-01_heldout-rerun-and-curator-names.md.
