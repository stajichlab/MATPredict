# Verticillium and proposed records: the five call changes
Status: open (curator review of the changes)

## Question
Which reference record causes each of the five call changes seen in the 11-record
regression panel, and is each change an improvement, an error, or neither?
Handoff item 1 of `docs/HANDOFF-2026-10-11-verticillium-refs.md`.

## Data and code version
- Code: commit 8e34533fac (frozen worktree `.claude/worktrees/run-8e34533fac-vd3`). `src/` is identical to e9060daf32; only db records and one note differ.
- Database for each variant: the e9060daf32 db (`.claude/worktrees/run-e9060daf32-base/db`, 69 Ascomycota records) plus the named records, built as symlink copies in `results/2026-10-12_item1_ablation/dbs/`.
- Inputs: the 5 genomes from `testset/regression_panel.tsv`.
- Committed-state panel (item 6): `results/2026-10-11_regression_vd3/` (3 records vs base).

## Method
1. Ran `detect` on the 5 genomes with 7 db variants: 3 Verticillium records singly (TV103, NBRC 6119, VdLs.17) and 4 species pairs (P. kewense 1132638 (syn. Eupenicillium crustaceum), P. chrysogenum 5076, F. sacchari 42676, F. fulva 5499). Jobs 29874859 to 29874865, all exit 0, 35/35 reports.
2. Compared called loci (contig, span, idiomorph, confidence, genes found) against base and against the all-11 run (`results/2026-10-11_regression_vd11/cand`).
3. For the changed core gene models, extracted the exon sequences from the genome, translated them, counted internal stops and read the splice sites.
4. For Bertia, translated the MAT1-2-1 and APN2 hits on the new locus and on the existing MAT1-2 locus and aligned them (BLOSUM62, local).

## Results
Taxon names were looked up in NCBI Taxonomy on 2026-10-12: taxid 1132638 is *Penicillium kewense* (synonym *Eupenicillium crustaceum*); 5076 *Penicillium chrysogenum*; 42676 *Fusarium sacchari*; 5499 *Fulvia fulva* (synonym *Cladosporium fulvum*). The ablation directory `G_ecrust` holds the *P. kewense* pair; the name was not changed on disk.

Cause of each change (a variant reproduces the change when it appears in that variant):

| Genome | Change | Reproduced by | Not reproduced by |
|---|---|---|---|
| GCA_001557505.1 Reticulascus | MAT1-1-1 hit changes from 5518_3639 (432 bp, 39.6%, 112 bits) to 27337_nbrc6119 (495 bp, 35.6%, 118 bits); call unchanged | NBRC 6119 alone | TV103, VdLs.17, all 4 pairs |
| GCA_025802115.2 Bertia | withheld -> called MAT1-2 (medium) at JAOYNC020017679.1:1278-5158; genes APN2 + MAT1-2-1 | TV103 alone; VdLs.17 alone | NBRC 6119, all 4 pairs |
| GCA_004918105.1 Aureobasidium | confidence medium -> high; MAT1-2-1 model changes | F. fulva pair (confidence, model); F. sacchari pair (span 32127 -> 12482 only) | Verticillium records, P. kewense, P. chrysogenum |
| GCA_026023515.1 Claussenomyces | gene MAT1-1-2 added to the MAT1-1 call | F. sacchari pair | all others |
| GCA_026023555.1 Chrysothrix | MAT1-2-1 model changes 264 bp/43.2% -> 252 bp/50.0% | P. kewense pair (252 bp, 100 bits); F. fulva pair gives 225 bp/53.3% | Verticillium records, F. sacchari |

Reading of each change (by sequence, from step 3 and 4):

- **Aureobasidium: improvement.** The base MAT1-2-1 model is 454 bp (not a multiple of 3) with 3 internal stop codons. The new model is 2 exons (GT..AG), 459 bp, 0 stops, and contains the HMG box (RPANAFILYR ... KIIGKMW). Identity 47.0% -> 49.3%, bits 85.5 -> 95.9.
- **Reticulascus: neither.** Both models are one exon with 0 stops. Both are partial alpha-box fragments with no Met start. The new model drops 10 residues at the N-terminus and adds 31 at the C-terminus. The new reference is from the same order (Glomerellales). Identity falls because the added tail is less conserved.
- **Chrysothrix: neither.** All models are partial HMG-box fragments, one exon, 0 stops, no Met start. The P. kewense model adds an N-terminal basic region (RQNK) and drops a C-terminal tail (RATRHSAGPRTM).
- **Claussenomyces: weak evidence.** The added MAT1-1-2 hit is 25.0% identity, 33.1 bits, 375 bp, no exon model (`exons: null`). It sits between MAT1-1-1 and MAT1-2-1 in a locus that already lists both.
- **Bertia: not resolved.**
  - First-pass resolution between MAT1-1-3 and MAT1-2-1 flips on a margin of 0.27 identity points (base, MAT1-1-3) to 1.66 (new, MAT1-2-1). The losing MAT1-1-3 hit has the higher coverage (51 against 35).
  - The MAT1-2-1 hit is a 93-aa HMG-box fragment, 35.9% identity, 65.5 bits, from reference 06cwm-g22 (no Verticillium reference is listed at this locus). The contig is 15,590 bp, and the locus starts 1,277 bp from the edge.
  - The genome already has a MAT1-1 call (JAOYNC020110414.1) and a high-confidence MAT1-2 call (JAOYNC020059261.1). The new call is a third locus.
  - The new locus is not a duplicate of the MAT1-2 locus. Local alignment of the two MAT1-2-1 proteins gives 25 identical residues (score 129). The two APN2 proteins give 187 identical residues (score 922) over a 375-aa and 569-aa pair.
  - This analysis cannot say whether the HMG-box fragment is a MAT gene or another HMG-box gene.

Side findings:
- G_pchrys gives bits 98.2 for the Chrysothrix MAT1-2-1 hit with the same reference, coordinates and identity that give 88.2 bits in base. I did not find the cause.
- The base APN2 on JAOYNC020059261.1 has 1 internal stop in its translation. This is outside item 1.

## What changed in detection
None. No code or db change in this study.

## Limits
- One variant per record group, 5 genomes. Interactions between records (for example TV103 with VdLs.17) were not tested. The two reproduce Bertia separately, and the all-11 run also has it.
- The ablation gives the cause of a change. It does not say whether the change is correct.
- Translations use genetic code 1. The genomes' codes were not looked up.
- No expression or RNA data. "Partial fragment" is judged from the translation only (no Met start, extent of the HMG or alpha box).
- Reticulascus, Chrysothrix and Aureobasidium were not checked against a protein database or domain model.

## Curator decisions
Open.
- Whether to accept the Bertia third-locus call (withheld before). Needs a judgement on HMG-box evidence at 36%.
- Eight proposed records (handoff item 8): the Aureobasidium, Claussenomyces and Chrysothrix changes come from them, so these results bear on that decision.

## Files
- `results/2026-10-12_item1_ablation/` (`runs/<variant>/runs/<asm>/detection_report.yaml`, `jobs.tsv`, `dbs/`)
- `results/2026-10-11_regression_vd3/diff/summary.md`
- `results/2026-10-11_regression_vd11/diff/ascomycota/regression_summary.md`
