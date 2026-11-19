# Weak-fragment cap: the 66 loci it would lower, checked against annotation (#341)

Date 2026-11-19. Follow-up to `2026-11-19_weak-fragment-cap-idiomorph.md` (PR #360). The cap stays off and in no preset.

## What was checked
The 66 `Ascomycota:MAT` loci (65 genomes; GCA_043170715.1 Cryan3 has two) that the redefined cap lowers in the v0.9.0 Ascomycota run (46 `high`, 20 `medium`). None is a FungiDB genome, so the FungiDB annotation of #290 does not apply. The independent annotation used is the BFD funannotate gene set (`Fungi_BFD_runs/genome_annotation/<genome>/predict_results/*.gff3.gz`; contig names equal the assembly's). For each locus the verdict is on the strongest core gene OF THE CALLED IDIOMORPH, the gene the cap ranks: `annotated_gene_present` (an annotated gene covers at least 50% of the model), `annotation_partial`, `no_annotated_gene`, `no_annotation` (no annotation folder was found for the genome).
Script `scripts/cap_annotation_review.py` (test `tests/test_cap_annotation_review.py`), job 32225678, commit b2980e6; table `results/2026-11-19_weak_fragment_cap_341/annotation_review.tsv` (summary next to it); source `results/2026-11-09_v090_ascomycota/reports_all.tar.zst`.

| Verdict on the called idiomorph's strongest core gene | Loci |
|---|---|
| `annotated_gene_present` | 55 |
| `no_annotated_gene` | 4 (GCA_018994105.1, GCA_022814695.1, GCA_025594185.1, GCA_964257465.1) |
| `annotation_partial` | 1 (GCA_937863395.1, a MAG) |
| `no_annotation` | 6 (GCA_002003505.1, XylFL0594, XyFL1777_1, XyFL0064_1, GCA_053525395.1, GCA_035927985.1) |

## What it does and does not say
- In 55 of 66 loci an annotated gene lies over the model that carries the call. The model is a 29 to 50 bit hit, so this says only that the region is inside a predicted gene. An HMG-box paralog would be annotated as well, so the check cannot tell a real MAT1-2-1 from a paralog, and an annotated gene is not evidence the thin model is the MAT gene. It does show the cap would not lower these on the ground that the region is empty.
- 56 of the 66 loci carry at least one flank gene (SLA2 or APN2) in `genes_found`; this is the pipeline's own evidence, not independent.
- The four `no_annotated_gene` loci (three `high`, one `medium`; 31 to 50 bits) are the cases where the cap's concern, a fragment with no gene behind it, is most plausible. They have not been examined further.
- Six genomes have no annotation verdict. I matched the folder by species and strain from `Fungi_BFD/samples.csv` and did not investigate why these six are not found (a naming difference is the likely cause); I did not look for the folders by hand.
- Several rows look like repeated assemblies of one species: the ASM2975xx/ASM2976xx assemblies (the first, GCA_002975535.1, is *Corynespora cassiicola*) give the same 45.4 to 47.4 bit MAT1-2-1 model. I did not group the 66 loci by species, so they are likely fewer independent observations than 66.

## Reading for the decision
The annotation does not argue for or against the cap as a tier rule: it has no way to separate a real MAT1-2-1 from an HMG paralog at 40 to 50 bits. What would: the domain content of the translated model (HMG-box for MAT1-2-1, alpha-box for MAT1-1-1) and flank synteny, as in the #261 check, or read typing where a species has a panel. These were not done. Until one of them is, the cap should stay out of every preset.
