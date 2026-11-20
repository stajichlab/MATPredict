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

## Domain content and flank synteny (added 2026-11-19; the check proposed above)
Script `scripts/cap_domain_flank_check.py` (test `tests/test_cap_domain_flank_check.py`), job 32229856, commit cd1a9ec; table `results/2026-11-19_weak_fragment_cap_341/domain_flank.tsv` (the four loci with no annotated gene first). For the called idiomorph's strongest core gene, the genome around its model (+/- 600 bp) was translated in six frames (stop-free stretches of 40 aa or more, no intron model) and searched with the Pfam model of the gene's domain at the gathering cutoff (MAT1-1-1: PF04769 alpha-box; MAT1-2-1 and MAT1-1-3: HMG-box PF00505, PF09011, PF14887). Flank synteny: a flank gene (SLA2, APN2, COX13) of the report on the same contig within 30 kb. Genomes from `Fungi_BFD_runs/input_clean_genomes/clean`.

| Verdict | Loci |
|---|---|
| expected domain found, flank gene within 30 kb | 47 |
| domain not found, flank gene within 30 kb | 19 |
| flank gene absent (either domain result) | 0 |

Every one of the 66 has a flank gene within 30 kb. The domain was found for 41 of 59 MAT1-2-1 genes, 4 of 5 MAT1-1-1 and both MAT1-1-3 genes. The deciding gene's bitscore does not separate the two groups (median 45.1 bits with the domain, 45.4 without; range 29.6 to 49.7).

**The four loci with no annotated gene**
| Genome | Tier | Deciding gene (bits) | Domain | Flank genes (distance) |
|---|---|---|---|---|
| GCA_018994105.1 | high | MAT1-2-1, polished_single (49.7) | HMG-box found (PF00505, 31.8 bits) | APN2 -3,998; COX13 -7,786 |
| GCA_022814695.1 | high | MAT1-2-1, polished_single (47.8) | not found | COX13 +1,467 |
| GCA_025594185.1 | high | MAT1-2-1, polished_single (31.2) | not found | APN2 -954; COX13 -3,565 |
| GCA_964257465.1 (MAG) | medium | MAT1-2-1, polished_single (45.1) | not found | SLA2 +4,039 |

**The 19 with the domain not found.** Eight are the GCA_002975535.1 to GCA_002976175.1 assemblies, which have the same MAT1-2-1 model (45.4 or 47.4 bits) and the same flank offsets (APN2 968 bp, COX13 3,267 bp from it); they are probably one species and count as one observation. The rest are single genomes (GCA_001752565.1, GCA_003788615.2, GCA_005059845.1 (MAT1-1-1, 29.6 bits), GCA_007644065.1, GCA_033117895.1, GCA_036971165.1, GCA_937863395.1, GCF_002105095.1 and the three above).

## What this settles and what it does not
- All 66 loci sit next to a flank gene, and 47 of them also encode the expected HMG-box or alpha-box domain at the model. That is evidence the capped calls are not empty fragments: the cap, as redefined, would lower a call whose gene model lies in a MAT-like locus and mostly encodes the right domain. It is not evidence for the call's idiomorph: an HMG-box paralog has the same domain, and no test here separates MAT1-2-1 from a paralog.
- The 19 without a domain hit are not shown to be defective. The HMG box of MAT1-2-1 can be split by an intron, the search had no intron model, a 40 aa minimum ORF and the gathering cutoff, so a thin or split box can be missed. A 29 to 50 bit model next to APN2 at 1 kb is the pattern of a divergent or partial MAT1-2-1, not of noise, but I did not test that.
- Neither check was run on loci the cap does not lower, so the 47 and 19 have no control group: I cannot say how many ordinary calls would also lack the domain at this window.
- Read typing was not used.

## Reading for the decision
The annotation, the domain search and the flanks together argue against treating the 66 as empty fragments, and nothing here shows that lowering them to `low` is warranted. The cap's premise (a thin fragment with nothing behind it) is best supported for the three `domain not found` loci with no annotated gene, and even those sit within 4 kb of a flank gene. I recommend keeping the cap out of every preset. If a rule on these is still wanted, a control (the same check on a random 66 of the loci the cap leaves alone) and an intron-aware domain search are the next measurements.
