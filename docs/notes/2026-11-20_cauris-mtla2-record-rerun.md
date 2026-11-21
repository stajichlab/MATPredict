# Rerun of the 313 FungiDB genomes with the corrected C. auris MTLA2 record (#364, #349)

Date 2026-11-20. Before/after for PR #364 (record `498019_b11243_MTL_A` version 4: MTLA2 as the 134 bp GT-AG intron, 160 aa), run as the repository rule requires for any change that can move a call. Same frozen code and preset as v0.9.0, so the only difference is the db record (`results/2026-11-20_cauris_mtla2_rerun/RUN.txt`). All 313 reports carry commit 7d73be5, preset v0.9.1 and the new db content hash 87c8afb5.

## Result: 307 genomes identical, 6 changed
Compared per locus: outer span, tier, idiomorph, and per gene status, span, exons, strand, identity and alternate-model tool. No tier, idiomorph or family changed in any genome.

| Genome | Gene or locus | Before | After |
|---|---|---|---|
| *C. auris* 6684 | MTLA2 | `polished_disagree`, exons 179339-179580, 179739-179952, identity 98.68 | `polished_agree`, exons 179339-179575, 179710-179952, identity 98.75 |
| *C. auris* B11245 | MTLA2 | `polished_disagree`, exons 1585220-1585461, 1585620-1585833, identity 100.0 | `polished_agree`, exons 1585220-1585456, 1585591-1585833, identity 100.0 |
| *C. auris* B8441 | MTLA2 | `polished_disagree`, exons 852451-852692, 852851-853064, identity 98.68 | `polished_agree`, exons 852451-852687, 852822-853064, identity 98.75 |
| *C. pseudohaemulonii* B12108 | MTLA2 | `polished_disagree`, exons 149346-149565, 149722-149963 (-), identity 52.0 | `polished_agree`, exons 149346-149594, 149727-149963 (-), identity 52.2 |
| *C. auris* B11243 (the record's own genome) | MTLA2 | `not_polish_candidate`, 14845-15461, identity 100.0 | `not_polish_candidate`, same span, identity 93.8 |
| *D. fabryi* CBS789 | MTL locus (alpha, high; MTLA2, MTLA1, MTLalpha1 unchanged) | outer span LMYN01000014:86090-89485 | outer span 86090-102372; `beyond_core_bp` 0 becomes 12,887 |

- The three *C. auris* genomes behave as predicted: the primary model now has the 134 bp GT-AG intron, equal to the former alternate, so the status is `polished_agree`. The introns are the ones the RNA-seq supports (`2026-11-19_model-disagreement-splice-check.md`).
- *C. pseudohaemulonii* B12108 changes model as well, which was not tested before: its MTLA2 now takes the exonerate model matching the new reference (intron 149595-149726, 132 bp) in place of the former 156 bp intron; it was one of the ten "tools disagree" genes where the motifs favoured miniprot with close identities, and no RNA-seq was checked for it. Identity to the reference is unchanged in practice (52.0 to 52.2).
- B11243's own identity falls from 100.0 to 93.8 because the reference protein is now the 160 aa model and the best protein hit covers the new exon boundary differently; its status and span do not change.
- *D. fabryi* CBS789 is an unexplained change: the outer span of its MTL locus grows by 12.9 kb past the core genes, with the same genes, tier and idiomorph, and `contig_edge_distance` 23,010 to 10,123. I did not trace which evidence produced the extension. It is not an effect of the MTLA2 intron (the genome's MTLA2 is unchanged); it is listed for the curator's approval or tracing, as the repository rule asks.

## What this does and does not show
- Calls are unaffected: 0 changes of tier, idiomorph or family in 313 genomes. The change is confined to the reported gene model of MTLA2 in four genomes and the span of one.
- It covers the 313 FungiDB genomes only; the 510 panel genomes were not rerun. Other genomes that use this record as a reference could move too.
- The *D. fabryi* span change needs a decision (approve as a legitimate extension, or trace and reverse). Nothing was changed in the repository for it.
