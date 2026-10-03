# Two idiomorphs in one genome: a neutral statement, tested on LCG
Status: open

## Question
When one family calls both a Plus and a Minus idiomorph in a genome, what is
the arrangement, and which causes does the evidence support? Curator's
ruling (2026-09-29): test first; this is NOT a confirmed homothallism signal.
Causes can be homothallism, hybrid/fusion, duplication, mixed culture or
heterokaryon, or assembly artefact.

## Data and code version
- Code: commits bcd1e0d, 48ee7ff on branch two-idiomorphs (from PR #9 076afe4).
  Module `src/MATPredict/detect/two_idiomorphs.py`; roster field
  `two_idiomorphs_report: true` on Mucoromycota:MAT.
- Calls: the LCG held-out reports made by frozen 076afe4
  (results/2026-09-28_lcg_holdout/runs2/). The statement was re-derived from
  those reports with the committed module; detection was not re-run, so no
  call could change.
- Genomes: /bigdata/stajichlab/shared/projects/ZyGoLife/LCG/Annotation/genomes/<Org>.sorted.fasta.
- Panel: the 36 LCG genomes with a determined Plus and Minus call, plus 6
  literature positives with one idiomorph called. Mycotypha africana NRRL
  2978 is excluded (a training strain).
- Literature: results/2026-09-29_mucoro_homothallism_literature/NOTE.md
  (Idnurm 2011; Gryganskyi 2018; Schulz et al. 2016).

## Method
1. Group determined calls by family; a statement is written when two
   idiomorphs are called. The best (highest confidence, then margin) Plus and
   Minus calls are compared.
2. Arrangement: same_locus (homothallic_candidate, or <= 50 kb on one
   contig), same_contig_distant (> 50 kb), unlinked (different contigs).
3. Evidence: confidence, classifier input and margin, flanks per call, contig
   GC, and the protein identity of any flank gene present at both calls
   (translated from the genome, global alignment, BLOSUM62).
4. Support rules (stated thresholds): >1 call of an idiomorph or a weak call
   (margin < 50 or fragment-typed) -> duplication (possible HMG paralog); GC
   difference > 5 points -> mixed_culture_or_heterokaryon; shared flank >= 95%
   identical -> hybrid_or_fusion (>= 99.5% also assembly_artefact);
   same_locus -> homothallism. Every cause is always listed, with its support
   (possibly none). Homothallism is never asserted.
5. Groups for checking: "reported homothallic" (Mucor genevensis, Rhizopus
   azygosporus, Zygorhynchus, Syzygites) vs "conventionally heterothallic"
   (Mucor hiemalis, racemosus, indicus/rouxii, Rhizopus stolonifer/nigricans,
   microsporus group, Backusella, Circinella). The grouping is general
   literature, not checked strain by strain.

## Results
Source: per_genome.tsv, statements.json (this folder).

- All 36 two-idiomorph genomes are `unlinked` (different contigs). None is
  same_locus or same_contig_distant.
- Groups: reported homothallic 9, conventionally heterothallic 16, unknown 11.
- **The rule proposed from the literature did not discriminate.** A shared
  flank whose copies are intact and divergent (55-95%, the Z. heterogamus
  pattern) supported homothallism in 13/16 heterothallic vs 5/9 reported
  homothallic genomes (per_genome_v1_no_weak_guard.tsv). Guarded by call
  strength: 5/16 vs 0/9 (per_genome_v2_weak_guard.tsv). It is now recorded as
  evidence only (commit 48ee7ff).
- **Weak second calls dominate.** 22 of 36 genomes have a call with margin
  < 50 (all but one a Minus call near 30 bits; just above the 25-bit floor):
  all 5 M. genevensis, both Z. exponens, 9 heterothallic, 6 unknown.
  Such a call may be an HMG paralog, not a sexM.
- **Ten genomes from different Mucor names share an identical score profile**
  (Plus margin 248.9, Minus margin 29.7, glrA copies 87.4%): M. hiemalis x2,
  M. indicus x2, M. racemosus NRRL 1427, M. rouxianus, M. rouxii,
  M. subtilissimus, Mucor sp. x2 (and Backusella ctenidia NRRL 6239).
  Identical scores imply identical proteins across nominal species: possible
  strain misidentification or one species complex. Not resolved here.
- Final supported causes: reported homothallic — duplication 7,
  hybrid_or_fusion 1, none 1; heterothallic — duplication 9,
  hybrid_or_fusion 6, none 2; unknown — duplication 6, hybrid_or_fusion 5,
  assembly_artefact 1, none 1. GC difference never exceeded 5 points.
- Literature positives:
  - Syzygites sp. MES 3091: Plus and Minus, unlinked, both margins strong
    (226/160), no shared flank; no cause supported (consistent with Idnurm
    2011, two separate loci).
  - Syzygites megalocarpus SC16: Plus called; a Minus locus on scaffold_38 was
    withheld at the fraction floor. Two loci are expected (Idnurm 2011); the
    second is lost to the floor.
  - Zygorhynchus heterogamus NRRL 1489 (the Schulz genome; sexM-sexP 5.3 kb
    at one locus): only Minus called (margin 71.3); no Plus call and no
    withheld Plus. A likely detection miss of sexP.
  - Z. moelleri x4: Plus only, consistent with Schulz 2016 (sexM not found).
  - Rhizopus azygosporus NRRL 13165: shared rnhA 95.1% -> hybrid_or_fusion,
    consistent with Gryganskyi 2018's fusion caveat.

## What changed in detection
Report-level only: a genome-level `two_idiomorphs` list (bcd1e0d, 48ee7ff).
No call, idiomorph or confidence changes (pipeline helper test; 912 tests).

## Limits
- Species groups are general literature, not strain-verified truth.
- Flank identity uses only flanks detection placed in both calls.
- Not assessed: degraded flank copies, genome-wide duplicated single-copy
  genes, read depth. GC is whole-contig.
- The 36 come from one sequencing project; 10 share one score profile.

## Curator decisions
- Made (2026-09-29): test first; neutral label; homothallism never asserted.
- Open: (1) whether weak Minus calls near 30 bits (22/36) should be treated
  as paralogs by a stronger typing floor; (2) trace the Z. heterogamus sexP
  miss and the Syzygites megalocarpus floor-withheld Minus locus; (3) check
  the 10 identical-profile Mucor genomes for strain identity.

## Files
- results/2026-09-29_two_idiomorphs/per_genome.tsv (final), statements.json,
  per_genome_v1_no_weak_guard.tsv, per_genome_v2_weak_guard.tsv, evaluate.py
