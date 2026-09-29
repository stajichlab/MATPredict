# Open questions and what each waits on

| Item | Waits on | Where |
|---|---|---|
| Sign-off of Umbelopsis Plus/Minus and S. racemosum NRRL 2496 records | re-measure on 076afe4; whether the guard is still needed under the MAT-gene gate | branch curation-umbelopsis; `results/2026-09-28_umbelopsis_rebased/` (running) |
| Mucor_Jena held-out results (65 strains) | detection runs; curator's taxonomy and mating types | `results/2026-09-28_mucor_jena_holdout/` (running) |
| ZyGoLife LCG held-out results (899 genomes) | leakage check and runs; curator's truth for non-Zygo genomes | `results/2026-09-28_lcg_holdout/` (running) |
| CAAX-dependent calls: review the unverified label | a labelled set of >=100 non-mating STE3 loci | `results/2026-09-28_validation_f3_f4/NOTE.md`; handoff review-later item |
| Group distant subloci by conserved flanks (mip/beta-fg) | design and test; S. commune Aα–Aβ ~450–550 kb apart | `results/2026-09-28_subloci_literature/NOTE.md` |
| Lichtheimiaceae MAT gene order and flanks | exploration/research; they lack Mucorales flanks | `results/2026-09-27_mucoro_curation_guard/NOTE.md` |
| Ascomycota idiomorph synonym map (MAT vs MATtub vs MATyl vs MATsc vs MTL labels) | curator-supplied mapping | Fable review Q5, `results/2026-09-27_locus_merge/NOTE.md` |
| F9: tests that assert roster state, not intent | rewrite | `results/2026-09-28_fable_review/README.md` |
| Receptor curation for Boletales, Hymenochaetales, rusts | locus deposits or precursor data | `results/2026-09-27_pheromone_positional/NOTE.md` |
| Sporidiobolales A1/A2 receptor classifier | later (separate task) | `results/2026-09-27_receptor_explore/NOTE.md` |
| Receptor copy number at a locus is not visible in reports | report design | `results/2026-09-27_russulaceae_receptor/NOTE.md` |
| Record assembly accession: 92 of 115 records are locus deposits with no assembly | decide how the self-check handles them | `results/2026-09-28_next_fixes/NOTE.md` |
| Shorter homeodomain-only redHD queries to cut HD cost | test | `results/2026-09-27_hd_prescreen/NOTE.md` |
| C. auris flank outside the idiomorph | curation | `results/2026-09-26_gap_zygosity_validation/` |
| Homothallism candidates (D. hansenii CBS767, C. kikuchii) | reads or literature | `docs/publication-notable-findings/` |
| Discovery-only lineages (Mortierellomycota, Kickxellomycota, Zoopagomycota, Glomeromycota, Blastocladiomycota, Chytridiomycota, Endogonales) | discovery projects, not curation | `results/2026-09-26_flank_synteny_ED/NOTE.md` |
| Pucciniales genome GCA_025617555.3 timed out | longer per-genome limit (GENOME_TIMEOUT) | `results/2026-09-26_basidiomycota_full/ANALYSIS.md` |
