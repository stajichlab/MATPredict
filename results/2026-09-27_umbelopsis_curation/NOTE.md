# Tier-2 Umbelopsis MAT records (2026-09-27, branch curation-umbelopsis, pending sign-off)

No Umbelopsidales MAT deposit or paper exists (NCBI, PubMed). Two genome-derived
records, in db/Mucoromycota/Umbelopsidales/; Mucoromycota:MAT scope now adds
Umbelopsidales (1913639).

- 41833_gzumbrama1_MAT_Plus: U. ramanniana gzUmbRama1, GCA_977110945.1
  (23 contigs, N50 1.65 Mb), CDSBDH010000020.1:157,014-174,784. sexP
  CAO3688886.1, 320 aa, single exon, PF00505 aa 108-172; classifier sexP 145.2
  vs sexM 34.2. glrA CAO3688830.1, tptA CAO3688890.1, algA CAO3688894.1.
- 44442_wa0000051536_MAT_Minus: U. vinacea WA0000051536, GCA_016758895.1
  (150 scaffolds, N50 1.33 Mb), JAEPRA010000016.1:141,140-158,574. sexM is not
  annotated: single-exon ORF 154,336-154,788 (-), 150 aa, PF00505 aa 33-94,
  no match in the Plus proteome; classifier sexM 39.8 vs sexP 26.6. 3' end
  uncertain (curated sexM are 185-249 aa). glrA KAG2174489.1, tptA
  KAG2174464.1, algA KAG2174465.1.
- rnhA is not at the Umbelopsis locus (present: false; ortholog ~500 kb away
  or on another scaffold).
- Correction during curation: the first Minus model fused the HMG exon to
  KAG2174494.1, a conserved neighbour present in both idiomorphs (89% to
  CAO3688882.1 beside sexP). Fixed in c0a4b71; error log B1.7/B1.8.
- Accessioned genes validate at 100% identity/coverage.
- Classifier rebuilt by scripts/build_idiomorph_hmms.py (it reads curated
  records): LOO 87/87, worst correct margin 22.5 -> 13.3 bits (the held-out
  Umbelopsis sexM). Roster min_margin unchanged at 25.

## Evaluation (293 Mucoromycota genomes, 4457c3a vs c0a4b71; compare_c0a4b71.txt)

Zygo 23: 23/23 locus and 23/23 idiomorph, scaffold and contig inputs.
Called genomes 231 -> 236; 273 unchanged, 6 gained, 1 lost, 9 label, 4 conf.

Umbelopsidaceae 8 -> 14/14 genomes with a call; routing phylum_fallback -> lineage.
- The 5 previously lost loci: 4 recovered as Minus/high (WA50703 x2, U. vinacea
  x2, core sexM now modelled, classifier margin 37); U. nana Minus/high (32.1).
- U. sp. M5902: Minus/low -> Plus/medium (model-based, margin 77.3).
- U. sp. AD052 (never called): undetermined/medium (margin 23.9).
- U. isabellina x4 and PMI_123: Plus, model-based, margins 76-143.
- gzUmbRama1 high -> medium (an unpolished distant rnhA hit sits in the cluster).

Side effects to review:
- 13 genomes gain a second 'undetermined/medium' call (classifier margin
  19.5-24.8, below the floor): 5 Apophysomyces, 7 Umbelopsis, U. ramanniana AG.
  These are 85-155 kb clusters of an HMG paralog plus flank-paralog hits at
  28-38% (tptA, algA, rnhA, glrA) now modellable against the new Umbelopsis
  flanks. Likely spurious.
- Rhizomucor pusillus GCA_900175165.2: Minus/medium -> uncalled (0 modelled genes).
- Mucor griseocyanus CBS116.08: Minus -> undetermined (classifier retrained).
