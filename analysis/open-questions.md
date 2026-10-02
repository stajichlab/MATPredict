# Open questions and what each waits on

Updated 2026-10-01. Closed rows (deterministic rebuild, Umbelopsis merge) are in decisions.md.

| Item | Waits on | Where |
|---|---|---|
| Absidia blakesleeana missed (flanks off-scaffold) | known gap (option c); revisit with the Lichtheimiaceae exploration | `results/2026-09-29_strain_labels_and_absidia/NOTE.md` |
| Lichtheimiaceae MAT gene order and flanks | exploration/research; hypotheses: count a strongly modelled rnhA next to the core as flank support; admit best-hit core-only clusters when the classifier types a full model above the build threshold (option b) | `results/2026-09-29_circinella_trace/NOTE.md`; `results/2026-09-27_mucoro_curation_guard/NOTE.md` |
| Strong-core floor rescue (e.g. S. megalocarpus Minus locus) | on hold; combine with the two_idiomorphs checks (adds 28 LCG loci, 21 as a second idiomorph) | `results/2026-09-29_sexM_like_paralog/floor_rescue.tsv` |
| Core-only admission for split loci (Z. heterogamus sexP on a 4.8 kb scaffold) | design and test with the gate and P1 class | `results/2026-09-29_sexM_like_paralog/NOTE.md` |
| Two-idiomorph genomes: resolve causes per genome | pseudogene test, duplicated single-copy genes, read depth (not assessed) | `results/2026-09-29_two_idiomorphs/NOTE.md` |
| P1 paralog class trained on one sequence; P1b and the Z. exponens gene not covered | more paralog copies | `results/2026-09-29_r4_paralog/NOTE.md` |
| CAAX-dependent calls: review the unverified label | a labelled set of >= 100 non-mating STE3 loci | `results/2026-09-28_validation_f3_f4/NOTE.md`; handoff review-later item |
| Group distant subloci by conserved flanks (mip/beta-fg) | design and test; S. commune Aα–Aβ ~450–550 kb apart | `results/2026-09-28_subloci_literature/NOTE.md` |
| F9: tests that assert roster state, not intent | rewrite | `results/2026-09-28_fable_review/README.md` |
| Future classifiers (Sporidiobolales A1/A2 first; Ascomycota MAT1-1-1/MAT1-2-1; Serinales MTLa/alpha; Basidiomycota HD later) | build with the chosen aligner; training-diversity check | `results/2026-09-27_receptor_explore/NOTE.md`; `2026-09-29_classifier-builds.md` |
| Receptor curation for Boletales, Hymenochaetales, rusts | locus deposits or precursor data | `results/2026-09-27_pheromone_positional/NOTE.md` |
| Receptor copy number at a locus is not visible in reports | report design | `results/2026-09-27_russulaceae_receptor/NOTE.md` |
| Shorter homeodomain-only redHD queries to cut HD cost | test | `results/2026-09-27_hd_prescreen/NOTE.md` |
| C. auris flank outside the idiomorph | curation | `results/2026-09-26_gap_zygosity_validation/` |
| Homothallism candidates (D. hansenii CBS767, C. kikuchii, Syzygites) | reads or literature | `docs/publication-notable-findings/` |
| Discovery-only lineages (Mortierellomycota, Kickxellomycota, Zoopagomycota, Glomeromycota, Blastocladiomycota, Chytridiomycota, Endogonales) | discovery projects, not curation | `results/2026-09-26_flank_synteny_ED/NOTE.md` |

## Added 2026-10-01
- Runtime check pending (`results/2026-10-01_runtime_check/`): held-out median
  ~260 s vs ~52 s per genome on unmatched nodes. Accuracy first; only
  call-neutral optimisations. Waits on: the same-node timing and profile.

## Future research: ploidy of Circinella strains with two sex-locus calls (tagged 2026-10-01, J. Stajich)

- Circinella muscae NRRL 1355, 1363 and 2403: each carries two Plus (sexP)
  calls on separate scaffolds. Circinella minor CBS 143.56 carries a Plus and
  a Minus locus.
- Candidate causes: two idiomorphs (homothallism or fusion), gene duplication,
  diploidy/heterokaryosis, or contamination (mixed culture). The curator notes
  that the same pattern in three C. muscae strains is strange, which argues
  against chance contamination.
- Proposed tests: genome-wide ploidy from read k-mer spectra and allele
  frequencies where reads exist; BUSCO duplication rate; read depth of each
  sex-locus contig; sequence identity of the two copies and their flanks.
- Evidence: results/2026-10-01_circinella_curation/NOTE.md;
  ANNOTATION_ERRORS_FIXED_REPORT.md entry C5 (branch curation-circinella).
- Waits on: availability of raw reads for these strains.

## Ruled 2026-10-01, work pending

| Item | Next step | Where |
|---|---|---|
| Ascomycota locus synonym map (MAT/MATtub/MATyl/MATsc/MTL) | draft for curator edit | Fable review Q5 |
| LCG name check | marker genes, Mucorales first | `results/2026-09-28_lcg_holdout/` |
| Absidia sp. NRRL 3163 identity | optional marker-gene check (not recorded as misidentified) | `results/2026-10-01_circinella_curation/NOTE.md` |
