# Genomes with both mating types, and homothallism in Mucorales
Status: decided (neutral report field); open (causes unresolved per genome)

## Question
LCG held-out genomes include 36 with both a Plus and a Minus call on different
contigs, among them species described as homothallic. What does the literature
say about homothallic Mucorales loci, and how should MATPredict report such
genomes without over-claiming homothallism?

## Data and code version
- Calls: LCG held-out reports (code 076afe4); `two_idiomorphs` re-derived from
  those reports with bcd1e0d / 48ee7ff (PR #9).
- Literature: The Mycota I, 3rd ed. (2016), ch. 10 and 16 (local EPUB; book text
  not copied); Idnurm 2011 (PMID 21908600); Gryganskyi 2018 (PMID 29674435);
  Lee & Idnurm 2017 (PMID 28332467); Schulz et al. 2016, Endocytobiosis Cell Res
  27:39-57 (local PDF; paraphrased only).

## Method
1. Literature review of homothallic Mucorales locus structure.
2. A neutral genome-level `two_idiomorphs` statement: arrangement (same_locus
   <= 50 kb, same_contig_distant, unlinked), per-call evidence, contig GC,
   identity of flank genes present at both calls, and `possible_causes`
   (homothallism, hybrid_or_fusion, duplication, mixed_culture_or_heterokaryon,
   assembly_artefact). It never asserts homothallism and changes no call.
3. Test on the 36 LCG genomes against literature positives and heterothallic
   species.

## Results

### Literature
Source: `results/2026-09-29_mucoro_homothallism_literature/NOTE.md` (incl.
the Schulz 2016 addendum).
- Three arrangements are described. (The first literature pass, before Schulz
  2016 was read, found no one-locus case; Schulz 2016 corrects that.)
  - Zygorhynchus heterogamus: one locus, sexM and sexP 5.3 kb apart (Schulz 2016).
  - Mycotypha africana: same scaffold, ~150 kb apart (Schulz 2016).
  - Syzygites megalocarpus: two separate loci, each with rnhA/glrA copies, one
    flank copy pseudogenised (Idnurm 2011).
- The tptA-sex-rnhA order is not kept in the homothallics studied.
- Plus/Minus protoplast fusion in Absidia glauca gives homothallic strains
  (Schulz 2016), so fusion cannot be excluded from sequence alone.
- Ch. 16 of the book has no Mucorales locus data; ch. 10 lists sexM/sexP in
  Syzygites without the arrangement.

### Test on 36 LCG genomes
Source: `results/2026-09-29_two_idiomorphs/NOTE.md`.
- All 36 are `unlinked`.
- The literature signal (divergent intact shared flanks) did not discriminate:
  it supported homothallism in 13/16 heterothallic vs 5/9 reported homothallic
  genomes. It is reported as evidence only.
- 22 of 36 carry a weak second call (mostly Minus near 30 bits); later shown to
  include the P1 paralog (see paralogs report).
- GC difference between the two contigs never exceeded 5 points.
- Literature positives: Syzygites sp. MES 3091 Plus and Minus, both strong;
  S. megalocarpus Minus withheld at the floor; Z. heterogamus called Minus only
  (sexP missed, split assembly); Z. moelleri Plus only; R. azygosporus supports
  hybrid/fusion. Mycotypha africana NRRL 2978 is a training strain and not a
  clean positive.

### Syzygites (notable finding 023)
Both Syzygites genomes carry sexP and sexM on separate contigs, and the sexM
protein is 100% identical between them. Evidence and caveats:
`docs/publication-notable-findings/023-syzygites-both-idiomorphs.md`.

## What changed in detection
bcd1e0d, 48ee7ff (PR #9): `src/MATPredict/detect/two_idiomorphs.py`; roster
field `two_idiomorphs_report: true` for Mucoromycota. Report-only.

## Limits
- Pseudogene detection, genome-wide duplicated single-copy genes and read depth
  are not assessed (listed under `not_assessed`).
- The degraded-flank signal rests on one species (Syzygites).
- The two Syzygites genomes may not be independent isolates.

## Curator decisions
- Made 2026-09-29: two-idiomorph genomes are not confirmed homothallics; test
  first; use a neutral label with possible causes.
- Open: whether a strong-core floor rescue should be combined with these checks.

## Files
`results/2026-09-29_two_idiomorphs/` (per_genome.tsv and variants);
`results/2026-09-29_mucoro_homothallism_literature/NOTE.md`;
`docs/publication-notable-findings/023-syzygites-both-idiomorphs.md`.
