# Open questions and what each waits on

Updated 2026-09-29.

| Item | Waits on | Where |
|---|---|---|
| Approve the deterministic full rebuild of the Mucoromycota sexM/sexP (explicit MAFFT L-INS-i) | replay on Mucoromycota, Zygo 23, LCG and Jena; curator approval | `results/2026-09-29_aligner_default/` (running); `2026-09-29_classifier-builds.md` |
| Merge curation-umbelopsis into PR #9 | P1 class and gate threshold added via fast paths; re-measure | `results/2026-09-29_umbelopsis_p1/` (running); `2026-09-27_umbelopsis.md` |
| Mucor_Jena and LCG truth tables | curator taxonomy and known mating types | `results/2026-09-28_mucor_jena_holdout/per_strain.tsv`; `results/2026-09-28_lcg_holdout/curator_table.tsv` |
| Absidia blakesleeana missed (flanks off-scaffold) | known gap (option c); revisit with the Lichtheimiaceae exploration | `results/2026-09-29_strain_labels_and_absidia/NOTE.md` |
| Lichtheimiaceae MAT gene order and flanks | exploration/research; hypotheses: count a strongly modelled rnhA next to the core as flank support; admit best-hit core-only clusters when the classifier types a full model above the build threshold (option b) | `results/2026-09-29_circinella_trace/NOTE.md`; `results/2026-09-27_mucoro_curation_guard/NOTE.md` |
| Strong-core floor rescue (e.g. S. megalocarpus Minus locus) | on hold; combine with the two_idiomorphs checks (adds 28 LCG loci, 21 as a second idiomorph) | `results/2026-09-29_sexM_like_paralog/floor_rescue.tsv` |
| Core-only admission for split loci (Z. heterogamus sexP on a 4.8 kb scaffold) | design and test with the gate and P1 class | `results/2026-09-29_sexM_like_paralog/NOTE.md` |
| Two-idiomorph genomes: resolve causes per genome | pseudogene test, duplicated single-copy genes, read depth (not assessed) | `results/2026-09-29_two_idiomorphs/NOTE.md` |
| Misidentified strains (B. ctenidia NRRL 6239, R. arrhizus NRRL 1470, T. repens NRRL 6240, the M. indicus lineage) | how to record them in BFD metadata | `results/2026-09-29_strain_labels_and_absidia/misidentified_strains.tsv` |
| P1 paralog class trained on one sequence; P1b and the Z. exponens gene not covered | more paralog copies | `results/2026-09-29_r4_paralog/NOTE.md` |
| CAAX-dependent calls: review the unverified label | a labelled set of >= 100 non-mating STE3 loci | `results/2026-09-28_validation_f3_f4/NOTE.md`; handoff review-later item |
| Group distant subloci by conserved flanks (mip/beta-fg) | design and test; S. commune Aα–Aβ ~450–550 kb apart | `results/2026-09-28_subloci_literature/NOTE.md` |
| Ascomycota idiomorph synonym map (MAT vs MATtub vs MATyl vs MATsc vs MTL labels) | curator-supplied mapping | Fable review Q5; `results/2026-09-27_locus_merge/NOTE.md` |
| F9: tests that assert roster state, not intent | rewrite | `results/2026-09-28_fable_review/README.md` |
| Future classifiers (Sporidiobolales A1/A2 first; Ascomycota MAT1-1-1/MAT1-2-1; Serinales MTLa/alpha; Basidiomycota HD later) | build with the chosen aligner; training-diversity check | `results/2026-09-27_receptor_explore/NOTE.md`; `2026-09-29_classifier-builds.md` |
| Receptor curation for Boletales, Hymenochaetales, rusts | locus deposits or precursor data | `results/2026-09-27_pheromone_positional/NOTE.md` |
| Receptor copy number at a locus is not visible in reports | report design | `results/2026-09-27_russulaceae_receptor/NOTE.md` |
| Record assembly accession: 92 of 115 records are locus deposits with no assembly | decide how the self-check handles them | `results/2026-09-28_next_fixes/NOTE.md` |
| Shorter homeodomain-only redHD queries to cut HD cost | test | `results/2026-09-27_hd_prescreen/NOTE.md` |
| C. auris flank outside the idiomorph | curation | `results/2026-09-26_gap_zygosity_validation/` |
| Homothallism candidates (D. hansenii CBS767, C. kikuchii, Syzygites) | reads or literature | `docs/publication-notable-findings/` |
| Discovery-only lineages (Mortierellomycota, Kickxellomycota, Zoopagomycota, Glomeromycota, Blastocladiomycota, Chytridiomycota, Endogonales) | discovery projects, not curation | `results/2026-09-26_flank_synteny_ED/NOTE.md` |
| Pucciniales genome GCA_025617555.3 timed out | longer per-genome limit (GENOME_TIMEOUT) | `results/2026-09-26_basidiomycota_full/ANALYSIS.md` |

## Added 2026-10-01
- Runtime check pending (`results/2026-10-01_runtime_check/`): held-out median
  ~260 s vs ~52 s per genome on unmatched nodes. Accuracy first; only
  call-neutral optimisations. Waits on: the same-node timing and profile.
- Thamnostylum lucknowense RSA_1015_Plus-T is called Minus (medium) against its
  Plus file label. Waits on: a trace (label vs paralog vs real Minus).
- Five Jena strains have no curator name: CBS169_57, CBS334_71, CBS417_77,
  CBS564_66, CBS608_78. Waits on: the curator.
- LCG names beyond the Jena-table matches cannot be verified (no framework).
  Waits on: a verification approach (e.g. marker-gene comparison).
