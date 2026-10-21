# `--homothallic-policy both_genes`: measurement (2026-10-21)

Patched run: 3,181 genomes (sets in `results/2026-10-18_homothallic_patch/sets.tsv`), `--partial-core-policy coverage_profile`, otherwise defaults.
Re-run locally from the committed reports with `compare_patch.py` and `analyze_patch.py` (outputs in `tables/`). Not run: the
legacy-cause breakdown (needs a table of per-genome failing gates that is not in the repo; `why_not.py` needs the HPCC).

## Result

- **No side effects on legacy output.** 2,595 Dothideomycetes loci appear in both this run and the 10-18 audit run: 0 differ in idiomorph, confidence,
  legacy class or genes. The patch changes only the locus class (`mat_locus` 63 or `idiomorph_gene_only` 5 -> `homothallic_candidate`; 68 loci, 68 genomes).
- **Controls untouched:** 0 of 1,430 baseline heterothallic Dothideomycetes (*Alternaria*, *Parastagonospora*, *Zymoseptoria*, *P. teres*, ...),
  0 of 40 *A. fumigatus*, 0 of 15 *Neurospora*/*Sordaria*, 0 of 47 two-locus genomes, 0 of *Macrophomina* (466) and *Botryosphaeria dothidea* (141).
- **Recall on the 298 conflicting genomes: 29 (10%).** *Didymella* 14, *Stemphylium* 5, other 10. **Not recovered: *Pyrenophora* 0 of 62 and
  *Aureobasidium* 0 of 119**, the two groups that dominate the conflict set. Why they fail is not measured here.

## Problem: false positives on heterothallic *Botrytis*

16 genomes whose BioSample states MAT1-2 (15 *B. cinerea*, 1 *B. fabae*) are re-classed `homothallic_candidate` (idiomorph MAT1-2, medium):
a `MAT1-1-1` gene scores 116-124 bits on the alpha profile and `MAT1-2-1` 221-223 bits on the HMG profile, 1,579 bp apart. This is the failure tracker's F3:
the `MAT1-1-1` is a 61% remnant (see `docs/notes/2026-10-17_botrytis-mislabels.md`), and the alpha floor (min re-classed value 116 bits) does not separate a remnant from a gene.
Of the 17 truth-set genomes re-classed, 1 is correct (*Sclerotinia minor*, stated homothallic or diploid) and 16 are wrong. Adopting the policy as it stands would
call heterothallic *B. cinerea* homothallic.

Other re-classed sets: flip-review 35 of 48 loci (the fused-locus verdicts: Polytolypa, Monilinia, S. nivalis, S. sativa re-classed; Calogaya not);
holdout 2 of 74 (`KEW_Xsp1`, `KEW_Xsp2`, both genes strong at 360-393 bits, 3,354 bp apart; not judged here).

## Decision for the curator

Default stays off, and it should not join the clean run as is. Options: (a) require the alpha gene to be full-length against the roster (remnant status, `remnant-status-step2`
removes remnant proteins from references) and re-measure; (b) raise the alpha floor above ~125 bits, which also drops the other 28 true candidates' margin
(min 116, median 189) so needs the distribution by set first; (c) leave the policy off and keep `undetermined`.
