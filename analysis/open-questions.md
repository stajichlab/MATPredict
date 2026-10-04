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

## Side note: reference genomes with foreign rDNA (found 2026-10-02, kept at curator's request)

The B12 ITS check (`results/2026-10-01_lcg_name_check/its_check/NOTE.md`)
ran barrnap + ITSx on the reference genomes too. In 4 of them, the rDNA does
not match the genome's name:

| Genome | rDNA best match | Use in MATPredict |
|---|---|---|
| BFD Cunninghamella bertholletiae 175 (GCA_000697215) | ITS 98.9% Rhizopus arrhizus CBS 112.07 (type); LSU 99.6% R. arrhizus | its sexP is a Mucoromycota MAT classifier training sequence (`training_extra.faa`) |
| BFD Lichtheimia ramosa B5399 (GCA_000738555) | ITS 98.5% Mucor circinelloides CBS 195.68 (type); LSU 99.8% M. circinelloides | none found in db/ |
| BFD Lichtheimia ramosa PG115-04A (GCA_037041495) | ITS 100% Nothophoma pruni (type); LSU 99.4% Epicoccum proteae (Ascomycota) | 2 proteins in `paralog_negatives.faa` |
| Jena Pilaira anomala CBS 695.68 (held-out set) | ITS 100% Mucor saturninus CBS 974.68 (type) | held-out; tree reference for the LCG Pilaira call |

- Only rDNA was checked. Nuclear markers were not. rDNA from a minor
  contaminant can assemble in a short-read assembly, so this is not yet
  evidence that the nuclear genome is misnamed.
- Check made 2026-10-02: the C. bertholletiae 175 sexP training protein
  matches the 3 other Cunninghamella sexP training proteins at 74-79%
  identity (blastp), and no Rhizopus sexP is in its top hits. The training
  label looks correct.
- Use: a diagnostic of reference-genome errors (contamination or
  mislabelled cultures). A possible routine check for reference genomes:
  rDNA identity vs genome name.
- Next step, if wanted: marker-protein placement of these 4 genomes; check
  the 2 PG115-04A paralog negatives are fungal Mucorales proteins, not
  ascomycete contaminant proteins.

## Ruled 2026-10-01, work pending

| Item | Next step | Where |
|---|---|---|
| S. pombe P and Yarrowia A/B idiomorph class (unassigned) | protein-domain evidence (Pc vs alpha box; MATA/MATB proteins) | `db/Ascomycota/order.yml` (B9) |

## Opened 2026-10-04 (campaign runs)

| Item | Next step | Where |
|---|---|---|
| R. toruloides HD: 2026-09-26 generic-HD locus (PR contig, ~200 kb from PR, with MIP1) vs v0.6.0 `redHD` locus (other contig, no MIP1); 26 genomes | curator ruling on which is the MAT-linked HD pair | `2026-10-04_basidiomycota-ascomycota-v060-campaign.md` |
| Basidiomycota PR calls: 1,360 unverified (strict-CAAX only); most phylum-fallback calls are PR only | keep the caveat in any summary; review with the receptor-loci study | same |
| Mucoromycotina locus size: flank pair differs by lineage | compare within one pair, or plot all with the pair marked | `2026-10-03_mucoromycotina-mat-campaign.md` |
| Held-out LCG/Jena genomes in publication figures | curator ruling | same |
| 3 outgroup HMG genes inside the sexP/sexM clades (Dicele1 h5, GCA_016758965.1 h6, Mycafr1 h1) | check locus and identity | same |
