# sexM/sexP candidate phylogeny, early-diverging scan

2026-09-26. Spec: `docs/superpowers/specs/2026-09-26-sexM-sexP-candidate-phylogeny-design.md`
(commit 22d4024). Input: `results/2026-09-26_early_diverging/` reports (code 634dda4).

## Summary

- The sexP clade is well supported: all 6 curated sexP references plus 65 other
  domains, SH-aLRT 98.1, UFBoot 99.
- The sexM references are not monophyletic. The best clade that has no sexP and
  no MAT1-2-1 holds 7 of 9 sexM references, at UFBoot 57. Two references fall
  outside it: 35722 ATCC11077 and 36080 R7B (*Mucor lusitanicus*). So no
  candidate can be called a sexM ortholog from this tree.
- **Mortierellomycota and Kickxellomycota:**
  - No candidate falls in the sexP clade.
  - 4 withheld candidates fall in the weak sexM clade: 2 *Mortierella alpina*,
    1 *Spiromyces*, 1 *Coemansia*.
  - All 13 calls fall in "other HMG" or have no HMG box.
  - This agrees with the flank synteny check: none of these calls is supported.
- **Lichtheimiaceae:** 2 of 79 withheld loci carry a sexP-clade HMG box. None
  falls in the sexM clade.
- **Syncephalastraceae:** 1 of 4 called loci carries a sexP-clade HMG box. None
  of the 25 withheld loci does.
- **A detection defect found on the way (tree-independent):**
  - 33 Mucoromycota calls find sexM and not sexP, yet are labelled Plus. All 33
    are *Rhizopus*: *R. arrhizus* 19, *R. delemar* 9, *R. stolonifer* 3,
    *R. microsporus* 2.
  - In these calls sexM hits at about 96% identity and 99% coverage, and the
    idiomorph resolution names sexM the winner. The idiomorph score still
    favours Plus, by about 870 to 380.
  - This explains a large part of the *R. arrhizus* Plus skew noted in the scan
    note. 19 of 43 *R. arrhizus* "Plus" calls have only sexM. See
    `label_check.txt`.

## Method

1. **Loci:** every detected or withheld locus in the reports whose genes_found
   has sexM or sexP. That gives 1,780 loci in 565 genomes (`loci.tsv`).
   Suppressed ASMIDs were skipped. Every report says genetic code 1.
2. **Extraction** (`extract.py`): each locus span ±5 kb. miniprot `--trans`
   with the 15 curated sexM/sexP proteins.
   - Where miniprot made no model, a tblastn fallback took one protein per
     3 kb HSP cluster (≥ 30 aa).
   - Result: 3,857 candidate proteins, 215 of them from miniprot.
3. **Outgroups:**
   - Same-genome HMG copies that overlap no reported locus: the 3 best per
     genome, in one genome per genus, 68 genera (`outgroup_clusters.tsv`).
   - Plus 2 Ascomycota MAT1-2-1 proteins, from *Fusarium* 5518 and *Tuber*
     55307. Both align across the full HMG box (69/69 match columns), so both
     were kept.
4. **HMG box and deduplication:**
   - hmmsearch against Pfam PF00505 (Pfam 38.2), the best domain per protein at
     i-Evalue ≤ 0.01, ±5 aa.
   - 1,241 proteins have an HMG box. 2,864 do not; they are mostly short
     tblastn fragments.
   - Identical domains were collapsed: 878 unique (`dedup_map.tsv`).
5. **Alignment choice** (`alignment_comparison.txt`):
   - A: hmmalign to PF00505, match columns only. 69 columns, 68 parsimony-informative, 12% gaps.
   - B: MAFFT L-INS-i plus ClipKIT smart-gap. 161 columns, 107 informative, 57% gaps.
     Counting only columns with ≤ 50% gaps, B has 76 columns (75 informative).
   - I kept A for the main tree. It has almost the same number of well-filled
     informative columns, and a quarter of the gap content. The extra
     B columns are mostly gaps outside the box.
   - B was run as a sensitivity tree.
6. **Tree:** IQ-TREE 3.0.1, ModelFinder (Q.INSECT+F+I+R6 by BIC), 1000 UFBoot,
   1000 SH-aLRT. Job 29114152, 5 h 15 min wall.
   Files: `tree_hmm.treefile`, `tree_hmm.contree`. Leaf IDs map to names in
   `taxon_ids.tsv`.
7. **Clade assignment** (`read_tree.py`): for each type, the smallest tree side
   that holds that type's references and no reference of the other type and no
   MAT1-2-1. For sexM, which is not monophyletic, I used the pure side that holds
   the most sexM references. `candidates.tsv` has one row per candidate
   protein. `clade_table.txt` has one row per locus: a locus counts as sexM or
   sexP if any of its proteins falls in that clade.

## Loci by clade (HMM tree)

The label columns count loci in the sexM or sexP clade whose detection label
matches the clade (Minus with sexM, Plus with sexP).

| group / family | status | loci | sexM | sexP | both | other HMG | no HMG box | label agrees | label disagrees |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Kickxellomycota | called | 7 | 0 | 0 | 0 | 6 | 1 | 0 | 0 |
| Kickxellomycota | withheld | 730 | 2 | 0 | 0 | 348 | 380 | 0 | 1 |
| Mortierellomycota | called | 6 | 0 | 0 | 0 | 4 | 2 | 0 | 0 |
| Mortierellomycota | withheld | 383 | 2 | 0 | 0 | 143 | 238 | 2 | 0 |
| Lichtheimiaceae | called | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0 |
| Lichtheimiaceae | withheld | 79 | 0 | 2 | 0 | 45 | 32 | 2 | 0 |
| Syncephalastraceae | called | 4 | 0 | 1 | 0 | 3 | 0 | 1 | 0 |
| Syncephalastraceae | withheld | 25 | 0 | 0 | 0 | 16 | 9 | 0 | 0 |
| Rhizopodaceae | called | 99 | 45 | 52 | 0 | 1 | 1 | 63 | 34 |
| Mucoraceae | called | 75 | 8 | 37 | 0 | 30 | 0 | 41 | 4 |
| Backusellaceae | called | 13 | 4 | 9 | 0 | 0 | 0 | 7 | 6 |
| Cunninghamellaceae | called | 19 | 4 | 9 | 0 | 6 | 0 | 13 | 0 |
| Umbelopsidaceae | called | 15 | 2 | 8 | 0 | 0 | 5 | 8 | 2 |

The other families are in `clade_table.txt`.

- **Label disagreements:**
  - 33 of the 34 disagreements in Rhizopodaceae are the *Rhizopus* labelling
    defect above.
  - The Backusellaceae (6) and Mucoraceae (4) disagreements are Minus calls whose
    HMG box falls in the sexP clade. This includes 4 *Apophysomyces* genomes.
    I did not check them further.
- **Withheld loci that carry a sexP-clade box:** 7 in Mucoromycota.
- **Non-locus outgroup copies:** of the 231, 7 fall in the sexP clade and 7 in the
  weak sexM clade. These may be MAT copies that no reported locus covers.
  They are not checked.

## HMG-box copies not included

These are copies with tblastn E ≤ 1e-3 to a curated sexM/sexP protein that
overlap no reported locus. Per genome (`hmg_copies_per_genome.tsv`), the
median is:

| Phylum | Genomes | Median copies in loci | Median copies not included (range) |
|---|---:|---:|---:|
| Mucoromycota | 293 | 1 | 10 (1-20) |
| Mortierellomycota | 100 | 1 | 6 (2-9) |
| Kickxellomycota | 190 | 1 | 8 (2-20) |
| Chytridiomycota | 24 | 0 | 4 (1-11) |

Only 231 of these copies entered the tree, as outgroups.

## Limits

- The alignment is one 69-column domain. Support is low across most of the
  tree, and sexM itself is not recovered as a clade (UFBoot 57). No placement
  in the sexM clade is evidence of orthology.
- The sexP clade is supported (UFBoot 99). Membership in it is the one strong
  result.
- 2,864 of 4,105 proteins have no detectable HMG box. Most are tblastn
  fragments from withheld loci, so "no HMG box" often means "fragment", not
  "absent".
- miniprot modelled only 215 of 3,857 candidates. The rest are HSP chains, which
  can merge exons wrongly or join paralogs.
- The outgroup sampling is one genome per genus and 3 copies per genome. Other
  paralog families may be missing.
- MAFFT sensitivity tree: not finished when this note was written. SLURM job
  29114153 (`tree_mafft.sh`) was at iteration 360 after 5 h 26 min, and the
  likelihood was still improving. Its time limit is 12 h. To score it, run
  `python3 read_tree.py tree_mafft` (this overwrites `candidates.tsv` and
  `clade_summary.tsv`, so copy them first).

## Files

- Inputs: `loci.tsv`, `clusters.tsv`, `outgroup_clusters.tsv`,
  `refs_sexMP.faa`, `refs_MAT121.faa`.
- Extraction: `extract2/` (per genome), `candidates_raw.tsv`,
  `candidates.faa.gz`.
  `extract/` holds a first run that I superseded, because its tblastn fallback
  kept one protein per window. It can be deleted.
- Alignments: `domains.faa`, `aln_hmm.afa`, `aln_mafft.clipkit.afa`, and the
  `*.ids.fa` inputs.
- Trees: `tree_hmm.*`, `tree_mafft.*`.
- Results: `candidates.tsv`, `clade_summary.tsv`, `clade_table.txt`,
  `tree_hmm.clades.txt`, `label_check.txt`.
