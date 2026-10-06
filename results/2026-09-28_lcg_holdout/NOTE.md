# LCG held-out test (Mucoromycotina), 2026-09-28

## Question
How does detect (PR #9 code 076afe4, frozen worktree run-076afe4) call MAT loci
in the ZyGoLife LCG genomes, used as a leave-out set, and do the old funannotate
annotations contain the MAT genes detect models?

## Data and code
- Genomes: /bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/genomes/<Org>.sorted.fasta
  (names match the annotation scaffolds and zygo_truth.tsv); annotation from
  annotate/<Org>/{annotate,predict}_results (proteins.fa, gff3).
- 897 genome folders (+3 script files). Scope by NCBI genus lineage
  (genus_taxonomy.tsv; Absidia and Echinosporangium assigned by hand):
  Mucorales 620, Umbelopsidales 1 = 621 in scope. Discovery-only, not scored:
  Kickxellomycota 171, Mortierellomycota 84, Entomophthoromycota 17,
  Endogonales 1, Basidiobolomycota 1, Ascomycota 1 (Alternaria).
- Command: `matpredict detect --genome <sorted.fasta> --phylum Mucoromycota`
  (no taxid), run_chunk.slurm, 2 short jobs (29198158/9), 16-way.
- Nothing entered the db, classifier or curation. No src/db edits.

## Leakage (leakage.tsv, protein_identity_leak.tsv)
| tier | n | called |
|---|---|---|
| clean (no strain match, no identical annotated protein) | 533 | 449 |
| clean strain, annotated MAT protein identical to a training sequence | 29 | 29 |
| clean strain, annotated protein contains a >60-aa training segment | 22 | 21 |
| in BFD, not in training | 8 | 8 |
| training leak (strain in curated record or classifier training) | 6 | 6 |
| Zygo 23 (known answers; excluded from classifier training) | 23 | 23 |
- Training leaks: Mucor circinelloides f. lusitanicus NRRL 3631 and Phycomyces
  NRRL 1555 (curated records); Apophysomyces ossiformis NRRL A-21654, Mycotypha
  africana NRRL 2978, Pilaira anomala RSA 1997, Radiomyces spectabilis NRRL 2753
  (training_extra genomes from BFD).
- Overlap with Mucor_Jena: 0 strains (by CBS/NRRL key). So no old-vs-new
  funannotate comparison was possible.
- 11 genomes have no collection number (XY*, JES_114, SC16, MES_3091, FLAS);
  they count as clean by strain but could not be checked by key.
- The identical-sequence check covers annotated proteins only; 56% of detect's
  MAT models are unannotated, so it undercounts.

## Results (summary.txt, per_genome.tsv, curator_table.tsv)
- 536/621 called (86.3%); 0 failed runs; wall median 162 s, max 378 s.
- Clean tier (533): 449 called: Minus 224, Plus 190, Minus+Plus 27+1,
  undetermined 7.
- Zygo 23: locus 23/23, idiomorph 23/23 (sorted.fasta input).
- 36 genomes carry both a Plus and a Minus call on different contigs, 23 of
  them high/high. They include genera/species described as homothallic
  (Zygorhynchus x2, Syzygites, Mucor genevensis x5, Rhizopus azygosporus) plus
  21 other Mucor and 10 Rhizopus. None is labelled homothallic_candidate.
  Not resolved: homothallism, mixed assemblies, or paralogs.
- 85 not called. 26 are "Absidia" corymbifera (22) and ramosa (4) and 8 are
  Lichtheimia: 34 Lichtheimiaceae, the known gap (no Mucorales flanks).
  Absidia blakesleeana 0/3 called although each has an annotated sexP scoring
  108-146 bits margin; withheld as below_fraction_floor + modelled_gene_bar.
  62/85 uncalled genomes have an annotated protein >=100 bits to sexM/sexP
  (45 with |margin| >= 25); these are candidate misses or paralogs.

## Old annotation vs detect models (model_vs_annotation.tsv)
Winning-allele models (gene matching the call's idiomorph): 576.
- Annotated at the locus: 254 (44%). Not annotated: 322 (56%) -- the
  annotation-gap pattern.
- Of the 254: 193 match at >=95% identity and >=90% coverage; 61 differ
  (mostly truncated annotated or model proteins at 77-88% coverage; one merged
  5-exon annotated gene).
- Classifier on the annotated protein agrees with the call's idiomorph:
  247/247 where the protein scored.

## Limits
- No truth beyond Zygo 23: no accuracy is claimed for the other 598.
- Strain matching is by collection-number keys; the identity check covers only
  annotated proteins.
- Folder names carry species names; taxonomy was used only for scope.

## Housekeeping
runs/ and lists/ hold the output of a cancelled first submission that pointed
at the GFF3 column by mistake (621 empty run folders). They can be deleted;
the removal was blocked by a safety check. Valid runs are in runs2/.

## Files
genomes.tsv, genus_taxonomy.tsv, leakage.tsv, protein_identity_leak.tsv,
annot_hmm_hits.tsv, per_genome.tsv, curator_table.tsv (blank curator columns),
model_vs_annotation.tsv, summary.txt, analyze.py, annot_hmm.py,
run_chunk.slurm, runs2/.
