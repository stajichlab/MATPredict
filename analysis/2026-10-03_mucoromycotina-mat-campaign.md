# Mucoromycotina MAT campaign: calls, locus size, synteny, sexP/sexM trees (2026-10-03)

Status: running (HMG-box tree and RAxML-NG check not in yet); figures are drafts.

## Question
Across BFD, LCG and Jena Mucoromycotina genomes: which genomes carry a sexP
(Plus) or sexM (Minus) MAT locus, how large is the locus per genus and
idiomorph, how is gene order kept across genera, and do sexP and sexM form
separate clades in a rooted gene tree?

## Ruling
J. Stajich, 2026-10-03: run Mucoromycotina first; locus size = inner ends of
the flanking genes; root the gene tree with a non-MAT HMG-box outgroup; a
synteny figure of about 20-30 pruned loci.

## Data
- 974 genomes: BFD Mucorales + Umbelopsidales 289 (288 after suppress), LCG 621,
  Jena 64. LCG and Jena are held-out sets (no training, curation or classifier
  use). Code: v0.6.0 (`run-7c7ed99`).
- Taxon overrides (`db/taxon_overrides.tsv`): confirmed rows are counted under
  the likely genus; `hold out` and `unconfirmed` rows are dropped from locus
  size and the tree set.

## Method
See `results/2026-10-03_mucoromycotina_mat/NOTE.md` for every step and file.

## Results
- Detection: 865 of 973 genomes called; 911 loci (Plus 482, Minus 425,
  undetermined 4). BFD 256/288, LCG 548/621, Jena 61/64. 45 genomes have more
  than one locus; 34 are flagged `two_idiomorphs`.
- Locus size (579 loci, 23 genera, polished flanks on both sides). Largest:
  Umbelopsis (median Plus 13,063 / Minus 12,256 bp) and Absidia (9,611 / 9,249).
  Intermediate: Phycomyces (6,531 / 4,166), Backusella (6,113 / 6,057), Pilaira
  (4,726 / 2,968). Most Mucoraceae and Rhizopus: about 1.4-2.6 kb. The flank
  pair is not the same in every lineage, so sizes compare one flank interval per
  lineage, not one fixed gene set.
- Synteny: 26 loci (Plus and Minus for 13 genera). Figure:
  `synteny/mucoromycotina_MAT_synteny.png`; interactive:
  `mucoromycotina_MAT_synteny.clinker.html.zst`.
- Full-length gene tree (526 sequences x 526 sites, JTT+F+R8, UFBoot 1000):
  - sexP clade: 163 tips, 0 sexM tips inside, UFBoot 40.
  - sexM clade: 173 tips, 0 sexP tips inside, UFBoot 93.
  - 3 outgroup HMG genes from Mucorales genomes fall inside these clades (2 in
    sexP, 1 in sexM); they may be MAT-gene paralogs away from the locus.
  - The outgroup is not monophyletic; 132 outgroup tips fall between sexP and
    sexM after rooting. This tree does not support or reject sexP + sexM as
    sister groups.

## Decision
None yet. Figures are drafts for the curator.

## Open
- HMG-box tree (job 29389343) and RAxML-NG check.
- Compare locus size within one flank pair, or show all with the pair marked?
- May held-out LCG/Jena genomes appear in publication figures? (They are used
  here for description only, not for any build or score.)
- LCG called 548 genomes on v0.6.0 versus 536 on 52b3ff9; not yet compared
  call by call.
- Identity of the 3 outgroup tips inside the sexP/sexM clades.
