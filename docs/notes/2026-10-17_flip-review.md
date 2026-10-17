# Review of the `--partial-core-policy coverage_profile` flips (2026-10-17)

Branch `flip-review`. Analysis only: nothing in `src/` or `db/` changed. Data and scripts: `results/2026-10-17_flip_review/`
(`flip_review.tsv` is the per-flip table, one row per vote note, 49 rows; `genes.tsv`, `loci.tsv` are the per-gene and per-locus
extractions; `mp_ranges*.tsv` the miniprot reference ranges; `bs.tsv` the BioSample fetch; scripts `extract_loci.py`,
`miniprot_ranges.py`, `genome_wide_boxes.py`, `fetch_biosample.py`, `build_tsv.py`). HPCC scratch: `/bigdata/stajichlab/jstajich/projects/flips/`
(jobs `flips-extract`, `flips-mp`, `flips-mp2`, `flips-gw`).

## Verdicts (46 flips, 2 undetermined, 1 edge-kept)

| group | n | verdict | confidence | what the evidence is |
|---|---|---|---|---|
| Botrytis, truth set (BioSample `mating_type` 1-2) | 17 | flip right | high | 17 of 17 correct on truth; margins 95-98 |
| Botrytis sample of 5 outside truth (T4 2011 assembly, B. hyacinthi, 3 B. cinerea) | 5 | flip right | high | identical locus to the truth genomes |
| Botrytis not individually examined | 9 | flip right | medium-high | identical gene table; no miniprot ranges |
| Botryotinia calthae | 1 | flip right | medium-high | Botrytis architecture |
| Stromatinia cepivora | 1 | flip right | medium | margin 66; 10 internal stops in the remnant |
| Ophiostoma (5 species) | 5 | flip right | medium-high | HMG gene full, "MAT1-1-1" is the C-terminal half with no alpha box; same remnant sits in the db's own O. montium MAT1-2 record |
| Leptographium, Intubia, Sporothrix | 3 | flip right | medium | short non-alpha-box "MAT1-1-1" beside full MAT1-2-1; no alpha-box ORF in the genome (Sporothrix: only the 60-bit locus ORF) |
| Monilinia vaccinii-corymbosi | 1 | BOTH present (fused) | high | all four MAT genes intact; published homothallism potential |
| Sclerotinia nivalis (2), S. sativa | 3 | BOTH present (fused, S. sclerotiorum Inv+ like) | medium-high | MAT1-1-5 intact plus MAT1-1-1 pieces both orientations |
| Sclerotinia spermophila | 1 | unclear | low | Botrytis-like, but Sclerotinia switching goes the other way |
| Polytolypa hystricis (undetermined) | 1 | both present; `undetermined` is right | high | published homothallic, fused |
| Calogaya cf. arnoldii (undetermined) | 1 | both present; `undetermined` is right | medium | metagenome; one colinear 63 kb contig |
| Aspergillus fumigatus (edge-kept) | 1 | not a flip; call itself suspect | low | one contig carries full MAT1-2-1 and a truncated MAT1-1-1 |

So 41 of the 46 flips look right, 4 (Monilinia, three Sclerotinia) are wrong in the sense that the locus holds both idiomorphs and neither label is right,
and 1 (S. spermophila) cannot be settled. No flip was shown to be right by independent evidence outside Botrytis; "right" for the others means the locus has
the same structure as a verified MAT1-2 (a complete HMG gene with a degenerate alpha-box remnant, and the other idiomorph's companion genes absent).

## What was examined and how

Per flip: the `on`-run report (`/bigdata/.../pcv/runs/on/<genome>/detection_report.yaml`, the run that produced the vote note); every gene of the locus
with coverage, identity, reference record, status and bitscore (`genes.tsv`); each protein translated from its exons (primary and alternate model, best
of the two) and scored on the alpha-box and HMG-box profiles with `score_core_profiles.py` functions (`alpha_best_of_models`, `hmg_best_of_models`);
internal stops; distance to the contig end; N in 3 kb; miniprot (`-u -N 3`) of the curated MAT1-1-1/MAT1-2-1 references (B. cinerea SAS56 MAT1-1,
SAS405 MAT1-2, S. sclerotiorum 44Ba1 for Sclerotiniaceae; O. montium UAMH 1363 and 11059 for Ophiostomatales) over the locus +-2 kb, to see which residues of
the reference the "MAT1-1-1" really covers; full BioSample attributes by E-utilities (`bs.tsv`); a genome-wide six-frame search with both profiles
(`gw_boxes.tsv`, section below). Read-level evidence was not used (not cheap; the Botrytis read check of the earlier note stands for Rf1 only).

## BioSample and independent statements of mating type

None of the 18 non-Botrytis genomes (incl. Botryotinia and A. fumigatus) has a mating-type attribute in its BioSample (`bs.tsv`; the `mating` grep hits are `material`/Botrytis BioProject
PRJNA624742). So there is no independent statement of mating type for any non-Botrytis flip. The only independent statements are in the literature,
all opened here through PubMed and cited by DOI; abstracts were read, not full texts.

- Monilinia vaccinii-corymbosi, reference genome (Yow et al. 2021, G3, doi:10.1093/g3journal/jkaa052): "The MAT locus sequences indicate the potential for
  homothallism in the reference genome, but a deletion allele of the MAT locus, characterized in a second isolate, indicates heterothallism." Our assembly
  GCA_017357885.1 (strain RL-1, 9 chromosomes, N50 4,062,641) is that genome (the report gives n50 4062641, 9 contigs, 30.0 Mb).
- Polytolypa hystricis (Munoz et al. 2018, Sci Rep, doi:10.1038/s41598-018-22816-6): "both of the basal branching, non-pathogenic species appear homothallic, with both
  mating type locus idiomorphs fused at a single locus".
- Sclerotinia sclerotiorum is homothallic with the idiomorphs fused end to end, Inv- and Inv+ alleles both self-fertile (Chitrampalam et al. 2013,
  doi:10.1371/journal.pone.0056895); in Inv+ MAT1-2-1/MAT1-2-4 are inverted and MAT1-1-1 is truncated at the 3' end, bounded by a 250-bp partial MAT1-1-1.
  S. trifoliorum switches unidirectionally by deleting the MAT1-2 genes, so a self-sterile strain keeps MAT1-1 (Xu et al. 2016, doi:10.1038/srep27083).
  Amselem et al. 2011 (doi:10.1371/journal.pgen.1002230): the loci differ in organisation and indicate evolution of heterothallism from homothallism.
- Monilinia fructicola/laxa/fructigena are heterothallic with MAT1-1 or MAT1-2 per isolate and small opposite-idiomorph fragments dMAT1-1-1/dMAT1-2-1 at the
  edge of the locus (Abate et al. 2018, doi:10.1094/PHYTO-03-18-0074-R). This is the same remnant structure as Botrytis, in a sister species; M. vaccinii-corymbosi
  differs by carrying both genes complete.

Not found or not verified: nothing for Stromatinia/Sclerotium cepivorum, S. nivalis, S. sativa, S. spermophila, any Ophiostoma, Leptographium, Intubia, Sporothrix
or Calogaya (PubMed searches for mating type returned nothing relevant; two Ophiostoma abstracts opened were not about these strains). The recollection that
Xanthoria parietina (Teloschistaceae) is homothallic was NOT verified (no hit).

## Botrytis: the sample of 5 and the rest

The 17 truth genomes are as in the earlier curation note. The 5 sampled non-truth genomes (GCA_000227075.1 T4 2011 assembly; B. hyacinthi Bh0001; B. cinerea
GCA_019186565.1 (chromosome-level), GCA_040955445.1, GCA_056956565.1) have, in miniprot ranges and gene table: MAT1-2-1 98-100% identical to SAS405 (KF944386) with MAT1-2-10
present, a "MAT1-1-1" that aligns reference residues 127-353 only (alpha profile 106-126 against 311 for the intact B05.10 protein, 1-5 internal stops), no MAT1-1-5, SLA2/APN2/COX13 in
order. Botryotinia calthae and Stromatinia cepivora, Sclerotinia spermophila give the same remnant shape (residues 105-139 to 348-353). The other 9 Botrytis have identical gene tables
(margins 95-98); only their tables were looked at. Verdict flip right for all.

## Sclerotiniaceae that are not Botrytis-like: Monilinia and Sclerotinia

Monilinia vaccinii-corymbosi RL-1, chromosome contig CP063406.1 of 5.2 Mb: SLA2 (1450531) - MAT1-1-1 model (3 exons, 1461302-1461789, 1470368-1470655, 1470705-1470963; 97.7% of the reference, alpha profile 228.5) -
MAT1-2-10 (1462116-1463046, 98.9%) - MAT1-2-1 (1463571-1464496, model 69%, HMG profile 188, miniprot spans reference residues 1-371) - MAT1-1-5 (1471217, 97.6%, 313 bits) - APN2 - COX13. Note that the
"MAT1-1-1" model bridges an 8.6 kb gap that contains MAT1-2-10 and MAT1-2-1: it is an N-terminal piece (1461302-1461789) joined by the aligner to a C-terminal piece (1470368-1470963), the alpha-box score 228 belongs to
the joined protein, and whether this is one functional MAT1-1-1 interrupted by the MAT1-2 block or two pieces (as in S. sclerotiorum Inv+, below) I could not tell. Either way all four MAT genes and both DNA-binding
domains are present in 35 kb. The bitscore vote said MAT1-2 (348 vs 225 bits) and the rule flipped it to MAT1-1; neither is right. The margin of 40.5 bits is the difference between two strong scores.

Sclerotinia nivalis SnTB1 and SM8 and S. sativa (same locus in all three, identical numbers): MAT1-1-5 intact (100%, 547 bits, 95% identity to S. sclerotiorum), MAT1-2-10, MAT1-2-1 (99.7%)
and MAT1-1-1 in two pieces: a C-terminal piece (reference residues 116-335, the 65% model the rule saw, alpha 158) and an N-terminal piece (residues 1-210 at 96% identity to S. sclerotiorum
MAT1-1-1, at 583584-584262 in SnTB1, 14281-14959 in SM8, on the opposite strand to the C-terminal piece in SM8). That is the Inv+ type of S. sclerotiorum (inverted MAT1-2 genes, MAT1-1-1
broken, MAT1-1-5 present). The profile margin of 83 bits compared the C-terminal fragment with a full HMG gene and so said nothing about the N-terminal alpha box 2 kb away. Verdict: fused locus,
neither MAT1-1 nor MAT1-2. Three assemblies from different laboratories agree (SnTB1 China, SM8 Russia, S. sativa). S. nivalis/sativa have never been checked by us for self-fertility.

Sclerotinia spermophila CBS 219.46: SLA2 - MAT1-1-1 (residues 115-348, alpha 175, 1 stop) - MAT1-2-10 - MAT1-2-1 - APN2, no MAT1-1-5 within +-2 kb, no N-terminal piece in miniprot. This
is the Botrytis MAT1-2 shape and the flip is probably right, but the S. trifoliorum result says Sclerotinia lineages derive from the fused locus by losing the MAT1-2 genes, not the MAT1-1 genes,
and a MAT1-2-only Sclerotinia would be unexpected. alpha 175 is also the highest remnant score. Unclear, low confidence; the genome-wide search below may help.

## Ophiostomatales and the others

Five Ophiostoma (O. piceae, O. perfectum, Ophiostoma sp. 15807, O. araucariae, O. pseudominus) have a full MAT1-2-1 (98.5-100% of the O. montium MAT1-2 reference, HMG 160-172, alpha 1-3) beside
a 47-54%-coverage "MAT1-1-1" (365 aa) that aligns reference residues 315-360 to the end of the 711-aa O. montium MAT1-1-1 (42% identity) and has alpha profile 5-18, i.e. no alpha box. The curated
O. montium MAT1-2 record in the db (`230073_uamh-11059_MAT_MAT1-2`) itself holds a 224-aa MAT1-1-1 remnant beside MAT1-2-1 (alpha 71, aligns 59-224 of it at 54%), so the shape is that of a MAT1-2 idiomorph
by the database's own curated example. O. piceae, O. perfectum and sp. 15807 have identical scores and are probably the same lineage. Ophiostoma sp. 15807 is an 8.2 kb contig (locus 1.2 kb from the end,
no SLA2): the least supported of the five. Leptographium lundbergii (112 kb contig), Intubia oerlemansii (40.7 kb contig) and Sporothrix brunneoviolacea have a full MAT1-2-1 (HMG 155-201) and a short "MAT1-1-1"
(17-35% coverage, alpha 0-60). These are low-medium: the evidence is the absence of an alpha box at the locus, and the Sporothrix and Leptographium hits align 500-660 residues of the O. montium reference at 35-36% identity
without any alpha-box signal. For Leptographium and Intubia the genome-wide search finds no alpha-box ORF; Sporothrix has one 60-bit ORF, at the locus, so a divergent MAT1-1-1 homolog there cannot be fully excluded.

## Genome-wide box search

`genome_wide_boxes.py` translates every contig in six frames, scores each ORF of >= 70 aa (unspliced, so intron-split proteins score low) on the two profiles, and reports the best
(`gw_boxes.tsv`, 25 genomes: all non-Botrytis flips except S. sativa, plus 6 Botrytis). Best alpha-box ORF per genome (bits): Ophiostoma (5) 20-24, Leptographium 0, Intubia 0 (no alpha-box protein
anywhere; the HMG ORFs score 80-107); Sporothrix 60 at the locus only; Botrytis 123-126 and Botryotinia 70 and S. spermophila 133 and Stromatinia 134 (the remnant itself, nothing better elsewhere); S. nivalis
182 and 159 (two ORFs on the locus contig, the N-terminal and C-terminal pieces); Monilinia 144 and 62; Polytolypa 204; Calogaya 266. So the MAT1-2 flips of Ophiostoma, Leptographium and Intubia are supported
by an absence of any alpha box in the genome (limits: unspliced ORFs; the assembly may be incomplete, Ophiostoma sp. 15807 is in 8 kb contigs), and S. spermophila has no second alpha-box ORF like the S. nivalis N-terminal piece.

## The four both-strong genera: what the loci look like

| genus | locus | reading |
|---|---|---|
| Monilinia | one contiguous chromosome, SLA2-MAT1-1-1-MAT1-2-10-MAT1-2-1-MAT1-1-5-APN2-COX13 | homothallic-type fused locus, all four genes intact; published as such for the reference genome |
| Sclerotinia | SLA2-(MAT1-1-1 C-piece)-MAT1-2-10-MAT1-2-1-(MAT1-1-1 N-piece)-MAT1-1-5-APN2 | fused locus, S. sclerotiorum Inv+-like |
| Polytolypa | SLA2-MAT1-2-1-MAT1-1-1-MAT1-1-4-APN2-COX13, 285 kb contig | fused, homothallic (published) |
| Calogaya cf. arnoldii | SLA2-MAT1-2-1-MAT1-1-1-APN2-COX13, 63 kb contig; BioSample is a lichen metagenome | most likely a fused homothallic locus in the mycobiont; a chimeric/mixed assembly cannot be excluded; MAT1-2-1 model has 7 internal stops |

None is a mixed assembly in the sense of two idiomorph contigs; none is a flip. Calogaya's profile scores: the pipeline's note has 323.4 vs 323.8; my reconstruction of the HMG protein
gives 139.6 (stops kept as X) so the HMG number is translation-dependent and unverified (both are strong, alpha 323 is robust).

Would the existing `homothallic_candidate` rule catch them? No, for a reason that is not the loci. `classify_locus` requires both genes to be in `full_length_models`, which excludes any gene
that took part in an idiomorph cross-match (every winner and loser). In 47 of the 49 loci the MAT1-1-3 reference hits MAT1-2-1 (MAT1-1-3 is an HMG-box gene) and the resolution ("winner MAT1-2-1,
loser MAT1-1-3") is recorded; MAT1-2-1 is therefore never "full length" for the rule. All four loci have it (Monilinia: MAT1-2-1/MAT1-1-3 identity 61/33; Calogaya, Polytolypa, Sclerotinia the same).
The other conditions are met in all four: same contig, MAT1-1-1 model >= 50% of its reference (97.7, 72.7, 62.7, 65.4 %), different `gene_class`, and for Calogaya, Polytolypa and Sclerotinia a gap of 2.0-2.3 kb (<20 kb) between the two genes. Monilinia would fail a second condition: the MAT1-1-1 model spans the MAT1-2-1 gene (separation <= 0). The
only `homothallic_candidate` among the 49 is the Aspergillus locus (MAT1-1-1 and MAT1-2-1, no MAT1-1-3 hit). `idiomorph_gene_only` is Ophiostoma sp. 15807 (no SLA2). So
the rule is effectively unreachable at Ascomycota loci where MAT1-1-3 cross-hits MAT1-2-1, which is most of them (47/49 here). This is a finding about the rule, not a change; not tested by simulation.

## Margin calibration

Profile margin = winner-side protein profile minus the other protein's, from the vote note (the pipeline's own numbers):

| genus (n) | margin (bits) | verdict |
|---|---|---|
| Botrytis (31) | 95.2-182.9 (17 truth: 95-98) | right |
| Botryotinia calthae | 179.0 | right |
| Ophiostoma (5) | 148.7-158.3 | right |
| Leptographium / Intubia | 185.1 / 191.3 | right (low-medium) |
| Sclerotinia spermophila | 228.6 | unclear |
| Sporothrix brunneoviolacea | 89.3 | right (low-medium) |
| **Sclerotinia nivalis x2, S. sativa** | **82.8, 83.2, 83.2** | **wrong (fused)** |
| **Stromatinia cepivora** | **66.2** | right (medium) |
| **Monilinia vaccinii-corymbosi** | **40.5** | **wrong (fused)** |
| Calogaya / Polytolypa | 0.4 / -14.0 | undetermined, correct |

Loci at 25-70 bits: two, Monilinia 40.5 (wrong) and Stromatinia 66.2 (right). A margin of 50 turns Monilinia into `undetermined` (correct) and keeps Stromatinia; 70 would make Stromatinia
`undetermined` too (a lost correct flip, but a cheap one). A margin cannot be calibrated further from this: a wrong flip sits at 83 (Sclerotinia) above a right one at 66 (Stromatinia) and 89 (Sporothrix), so
no single threshold separates them (n = 4 wrong, one clade each). The margin measures how different two protein scores are, not whether the loser is a remnant: S. nivalis' 83 compares a
C-terminal fragment to a full HMG protein. Proposal: **50 bits** as a floor (confidence ~40%: it is the value that separates the one observed boundary case, from two data points), and rely on two
vetoes that do separate the groups in these data:

1. loser-companion veto: if a gene specific to the set-aside idiomorph and not a remnant is intact in the locus (MAT1-1-5 for the alpha side in Leotiomycetes: Monilinia, S. nivalis, S. sativa, all 100/97.6% coverage; absent in all 31 Botrytis, Botryotinia, Stromatinia, S. spermophila) the call is `undetermined` (or a new homothallic class). 4/4 wrong flips caught, 0 right flips lost (confidence in the rule: moderate, because it rests on one gene and one family; MAT1-1-5 is Leotiomycete-specific).
2. both-strong veto: set-aside protein >= 200 bits and winner >= 150 bits (Monilinia 228/188, Polytolypa 261/247, Calogaya 323/>=140). No right flip has a set-aside protein above 175 (S. spermophila) or 131 (others). It misses S. nivalis (159), which veto 1 catches.

I am least sure of the thresholds in veto 2 (three data points). The curator's proposal ("keep loci where both score strongly as `undetermined`") is supported, with the caveat that the strength has to be measured on the whole
region, not on the model the vote kept (S. nivalis' alpha box is in a piece the model did not take).

## Recommendation on the default

Do not turn it on globally. Evidence: 17/17 on truth for Botrytis and 5/5 sampled Botrytis consistent; it is also probably right for Ophiostomatales, Botryotinia, Stromatinia (41/46 flips), but 4 of the 46 flips are wrong (Monilinia, three
Sclerotinia) and they all fail in the same way, a fused homothallic locus mislabeled as a single idiomorph, while the old call was equally wrong for them. So:

- **Safe to enable now for Botrytis/Botryotinia** (the only clades with truth); a clade allow-list is the only option here that rests on evidence.
- **Ophiostomatales, Stromatinia: probably right (medium-high / medium)**: no alpha-box ORF anywhere in the Ophiostoma, Leptographium and Intubia genomes; no independent statement exists; enable after the companion veto exists.
- **Not for Sclerotinia or Monilinia** until the vetoes are in; the rule is wrong on 4 of 5 loci there (S. spermophila unclear).
- Before default ON campaign-wide: add the loser-companion veto and the both-strong veto, then re-run the 54 triggers; and address the `homothallic_candidate` cross-match exclusion so fused loci get a class rather than `undetermined`.

## Not verified / caveats

- No flip outside Botrytis has independent truth other than the Monilinia and Polytolypa literature statements (abstracts only, same-assembly lineage inferred for Polytolypa, exact for Monilinia by N50).
- Sample of 5 Botrytis only; 9 Botrytis never looked at beyond the gene table. S. sativa was not in the genome-wide search (same locus numbers as S. nivalis).
- Calogaya HMG score depends on the translation (139.6 vs 323.8); the metagenomic origin means the locus may not be one genome.
- Leptographium, Intubia, Sporothrix: the short hit plus the genome-wide ORF search (unspliced ORFs; a divergent alpha box could be missed).
- Ophiostoma: no strain-level statement; the "remnant in the db's own MAT1-2 record" is the curated record's content, taken from its proteins.faa (not independently re-derived from GenBank).
- S. nivalis/sativa structure inferred from miniprot ranges (no full gene model of the N-terminal piece, no read check, no RNA).
- Read-level evidence: none.
- Whether `homothallic_candidate` would fire was reasoned from the code and the report fields, not by running the code on a patched build.
- The 6 loci that no longer trigger and the 5 non-triggered `undetermined`-adjacent loci were not examined.
