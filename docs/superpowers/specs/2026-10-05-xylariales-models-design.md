# Xylariales MAT models: design (proposal)

Status: proposal for curator review. Nothing here is implemented.

## Problem
Xylariales have no known canonical MAT locus (Robinson & Natvig 2019; see
`docs/notes/2026-09-24_mat-reference-gap-literature.md`, Xylariales). Two
observations from our runs sit awkwardly with that:

- *Xylaria flabelliformis* NC1011 (GCA_022453505.1) gets one call, MAT1-2/medium
  (`idiomorph_gene_only`), for a 734-aa HMG protein (KAI0192626.1) matching
  MAT1-2-1 at 44%, among housekeeping genes. Its intact SLA2-APN2 block carries
  no MAT gene and is withheld.
- Other Xylariales (*Microdochium paspali*, sp. XT01, *Didymobotryum rigidum*)
  carry a MAT-like gene between SLA2 and APN2 at 39-47 bits.

The paper's data (`docs/notes/2026-10-05_robinson-natvig-2019-xylariales.tsv`)
shows the same split: 16 of 35 genomes have a MATA_HMG gene linked to
SLA2/APN2, 19 do not, and nearly all linked genes are closest to *N. crassa*
NCU03481, a non-MAT sexual-development regulator.

## Goal
Call an HMG gene as a MAT1-2-1 candidate only when it falls in the MAT1-2-1
clade or is positionally and phylogenetically supported, and report
"MAT-like, idiomorph not determinable" otherwise.

## Components

### 1. MATA_HMG family tree and clade models
- Alignment of the HMG core (residues 122-233 of *N. crassa* Mat a-1,
  AAA33598, plus ~35 downstream residues) for: known MAT1-2-1 and MAT1-1-3
  proteins already in the database; NCU03481/PaHMG8 and fmf-1/PaHMG5 orthologs
  across Pezizomycotina; the Xylariales HMG proteins from the paper's S2/S3.
  The paper's alignment is on TreeBase (S23036) and can seed this.
- Tree with RAxML-NG (tooling from the `raxml-hmg-check` analysis; note that
  analysis found the sexP clade not robust, so robustness of the MAT1-2-1
  clade must be tested first).
- Profile HMMs per clade: MAT1-2-1, MAT1-1-3, NCU03481-like, fmf-1-like.

### 2. Paralog decoys in classification
Each HMG hit is scored against all clade models. It supports a MAT call only
if the MAT1-2-1 (or MAT1-1-3) score beats the paralog scores by a margin set on
the benchmark. Candidates that place with NCU03481/fmf-1 are reported as
`hmg_regulator_like` and not as an idiomorph. Relation to the existing decoy
work: `docs/decoy-feasibility.md`.

### 3. Xylariales flank model
SLA2 and APN2 are linked in 28 of 35 genomes but their order and orientation,
and the HMG gene's position, vary. The Xylariales profile should accept either
order and orientation, and must not require the Sordariomycete layout.
Positional seeds: the 16 linked genes in the TSV (verify each against the
genome before use, since gene models were AUGUSTUS-predicted).

### 4. Order-level outcome
For Xylariales the report gains a category: "no canonical MAT; MATA_HMG
candidate at SLA2-APN2, clade placement X". The reason is recorded, with the
alternatives from the paper (MAT1-1 lost with MAT1-2-1 diverged; upstream
PaHMG5/PaHMG8 control; unisexual reproduction).

## Checks before any model is trusted
1. **NC1011 interval.** Translate six frames and run tBLASTn on the SLA2-APN2
   interval, using the Xylariales MAT proteins as queries, and map RNA-seq
   (SRR8861595) over it. The spec's gene order places COX13 between SLA2 and
   APN2 (unverified); if so the interval is ~2.1 kb and shared with COX13.
2. **KAI0192626.1 placement.** Place it in the HMG tree. Expectation (not
   tested): NCU03481-like.
3. ***D. rigidum*.** Confirm MAT1-1-2 (PF17043) and the position; the paper
   found no MAT1-1-2 in any member.
4. **Order-wide tally** over the 257 BFD Xylariales: flank block intact with MAT
   gene, intact and empty, flank block absent; and for the HMG genes, clade
   placement.
5. **Seed audit.** Re-run the paper's 35 genomes (accessions in the TSV) with
   the current pipeline and compare to Table S2 positions.

## Decisions for the curator
1. Should HMG-only calls in Xylariales without clade support be downgraded or
   withheld now (NC1011 is the case), ahead of the models?
2. Is `hmg_regulator_like` the right label, or should such hits be silent?
3. Pestalotiopsis (Amphisphaeriales in our taxonomy) is in the paper's set:
   keep it in the Xylariales seed set or hold it out?
