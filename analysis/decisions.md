# Curator decisions (J. Stajich), in date order

Source: the session ruling log (kept in the assistant's project memory),
`docs/HANDOFF-2026-09-26.md` and commit messages. Evidence paths are given in
each entry. Earlier rulings replaced by later ones are listed first.

## Superseded rulings

| Earlier ruling | Replaced by |
|---|---|
| Flank-carried: core hit inside flank span ±3 kb (09-26) | Strongest core hit, per-family window (09-26), then bitscore >= 39 instead of E<=1e-5 (09-27) |
| Rhizopus Plus skew is likely sampling bias (09-26) | Labelling artefact from the btbA vote (09-26) |
| Lichtheimiaceae/Syncephalastraceae record only at UFBoot >= 95 (09-26) | Classifier margin >25 bits plus Mucorales gene order (09-27) |
| Fallback confidence floor target 40% (09-27) | No identity floor (09-27) |
| CAAX unverified label for 4 Agaricales families (09-27) | All CAAX-dependent calls unverified; review at >=100 non-mating STE3 (09-28) |
| Group A: never merge Aalpha with Abeta, provisional (09-28) | Subloci ruling: one A and one B call, subloci as evidence (09-28) |
| btbA Plus-only (09-20); review pending (09-27) | btbA non-voting and present in both idiomorphs (09-27) |
| Fixed MAT-gene gate threshold 100 bits (09-28) | Threshold set by each classifier build from 189 paralog negatives (09-29) |
| Secondary-undetermined guard (09-27) | Dropped: the MAT-gene gate made it redundant (09-29) |
| Two-idiomorph genomes flagged as homothallic candidates (proposal, 09-29) | Neutral `two_idiomorphs` field with possible causes; homothallism never asserted (09-29) |
| Circinella/Thamnostylum Plus-labelled strains called Minus: detection error suspected (09-29) | Real sexM; labels disputed, origin unknown (09-29) |
| MAFFT `--auto` in classifier builds (09-26) | MAFFT L-INS-i set explicitly, single-threaded, deterministic (09-29) |

## Log

### 2026-09-26 (first round) They unblock the Basidiomycota launch.

1. Per-family polish cap default = 6 (genes-first rank).
2. Flank-carried calls (no modelled core gene): core hit inside flank span +-3 kb -> keep, cap low, partial_locus, idiomorph_unmodelled flag; outside -> withhold.
3. Phyla with no curated family: default routing is "not searched". Search only with an explicit override. No more default exhaustive routing.
4. Scoring plan approved: score genome genotypes against three truth kinds (curated deposits/holdout; independent calls such as Zygo abspres, C. albicans reads, Saccharomyces cassettes; negative controls). Report recall with no_reference cases removed, idiomorph accuracy on found loci, and per-confidence-tier precision. Culture-collection (+)/(-) strain labels are a candidate truth source (not yet checked).
5. BFD suppress list is `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt` (CSV `ASMID,REASON`). There is no `suppress.txt` at the BFD root. Nine PRJEB104476 "DAH" assemblies (single Sanger amplicons, 356-655 bp, NCBI "Complete Genome", method Chromas) were added on 2026-09-26. The curator said GCA_986280975.1 should be removed completely as a broken NCBI upload.
6. MATPredict must honour a local suppress list too.
7. Start the flank-ortholog synteny check for Mortierellomycota/Kickxellomycota (results/2026-09-26_flank_synteny_ED/).

Why: these were the three blockers in docs/HANDOFF-2026-09-26.md plus data hygiene found in the early-diverging scan.


### 2026-09-26 scoring details
- Unit = the whole genome's mating type. A homothallic genome with two MAT loci is scored as its own genotype class.
- Zygo absence/presence calls count as independent truth. Caveat: they came from a BLAST-based method, so shared errors could inflate agreement.
- Culture-collection strain labels count as truth. They are usually in the strain name: "+", "-", "Plus", "Minus", rarely "P" or "M".
- Chytrids are a negative control for the current curated DB only. Chytrid mating-type discovery is a long-term goal (idiomorph or allelic diversity).
- The six GCF_TEST stub files in the BFD input_clean_genomes folders were deleted on 2026-09-26. The samples.csv rows are NOT edited: samples.csv is regenerated from NCBI Datasets, so suppress.txt is the only exclusion mechanism.
- Rhizopus arrhizus Plus skew: WITHDRAWN as sampling bias. It is largely a labelling error: btbA (marked Plus-only from one record) outvotes sexM by bitscore in idiomorph_candidates; 33 sexM-only Rhizopus calls labelled Plus (found 2026-09-26).
- Commit and push PR #9 after the code changes pass their tests.

### 2026-09-26 third round
- Flank-carried rule: option (c). Audit the Ascomycota calls it drops or withholds before any sweep uses it. Ascomycota may need different parameters than Serinales.
- C. albicans single-idiomorph calls: report zygosity unknown. A better treatment is needed later.
- Report `assembly_gap_at_locus` when an N run sits at the flank-anchored MTL position.
- C. auris reference-built isolates: read-based confirmation is deferred until paper writing.
- Mark reference-built assemblies and exclude them from scoring sets.

### 2026-09-26 fourth round
- Basidiomycota HD allele labels: use "A1", "A2", ... counting when the source does not number alleles. A future HD A-gene phylogeny may resolve them (not scoped).
- MIP1/beta-fg stay as OPTIONAL flanks (+36% runtime, corroboration only). They hold in Agaricomycetes only; the curator's belief that anchors work in Pucciniomycotina was not supported (0/7).
- homothallic_candidate must not fire in Basidiomycota: Agaricomycotina are not homothallic in the same way, and HD1+HD2 / bE+bW are partner genes at one locus.
- Flank-carried fix: strongest core hit, E <= 1e-5, per-family window (3 kb MTL, 20 kb SLA2/APN2/COX13 families and Mucoromycota).
- Mortierellomycota/Kickxellomycota calls: label unverified (no synteny support; not clear what they are yet).
- A C. auris flank outside the idiomorph is queued for curation.
- basidio-anchors was rebased onto polish-scope-cuts (c6ebcda) and pushed as its own branch.

### 2026-09-26 fifth round
- SXI1/SXI2 slot is a BONUS slot (accepted after first approving a required slot).
- PH4021C (GCA_025531435.1) and CBS 132 are AD hybrids; their a+alpha calls are expected, not homothallism.
- Phycomyces NRRL_1554: a real regression (Plus/high on 09-21, lost by 09-23). Being traced (results/2026-09-26_phycomyces_trace/).
- The four Dothideomycetes records (Cochliobolus C4/C5, Zymoseptoria IPO323/IPO94269) are signed off. Branch curation-mucor-dothideo is pushed.
- Diplodia sapinea: skipped (indirect evidence; Botryosphaeriales already 8/8).
- Curate Leptosphaeria, Parastagonospora, Pseudocercospora and Fulvia next.
- Lichtheimiaceae/Syncephalastraceae: propose tier 2 only if the sexM/sexP tree places one of their HMG genes at UFBoot >= 95.

### 2026-09-26 evening
- btbA: idiomorph_informative false (it outvoted sexM by bitscore).
- Model BOTH sexM and sexP (curator: they are hard to tell apart); decide on the two models.
- Fix the relaxed-pass polished_genes defect (6bc985d); explore the medium cap on data.
- Polish cap: test a mixed rank (top 5 by genes + top identity cluster) by replay; adopt only if it loses nothing.
- Zygo regression runs on both inputs: scaffolds+proteome and contigs.fsa (mapped via .agp), reported as two scores.
- Build sexM/sexP HMMs from confident sets and test discrimination on MODELLED proteins (avoids the earlier 6-frame failure; see matpredict-hmm-typing-rejected). The sexM clade is not supported (UFBoot 57), so sexM training cannot come from tree membership.
- Pucciniales genomes run ~43 min median each (some hit the 1 h per-genome limit); the 0.8 h job estimate was wrong. Resubmitted as 5 jobs of 20.
- Tree strategy (curator, 2026-09-26): iterate with FastTree while candidate sampling changes; run full ML (IQ-TREE and/or RAxML-NG) only once sampling is settled. The MAFFT sensitivity IQ-TREE job 29114153 was cancelled at 8.5 h. Cluster modules: fasttree/2.1.11, raxml-ng/2.0.2, hmmer/3.4, mafft/7.505.
- HMM idiomorph classifier (curator, 2026-09-26): add it as a general step that runs after gene modelling. It scores each modelled core protein against per-idiomorph HMMs; the call goes to the higher score; close calls under min_margin are undetermined. The bitscore vote remains the fallback. Files go in db/<Phylum>/classifiers/<family>/ with manifest.yaml, built only by scripts/build_idiomorph_hmms.py (never edited by hand), using pyhmmer. Training option (a): curated DB refs + supported sexP-clade members only; Zygo 23 stays OUT as an independent test. Start with Mucoromycota; extend to other groups if it works. Test result: 38/38 leave-one-genus-out, worst margin 36.6 bits full-length (results/2026-09-26_sexMP_hmm/).

### 2026-09-26 night, Basidiomycota (full run analysed: results/2026-09-26_basidiomycota_full/ANALYSIS.md):
- No full Basidiomycota re-run. A cap-off test on ~50 fallback-order genomes decides whether the fallback orders get re-run (results/2026-09-26_basidio_capoff_test/).
- Curate Sporidiobolales and Pucciniales now; check Wallemia honestly (branch curation-puccinio, from basidio-anchors).
- Pheromone-receptor (PR) curation is REQUIRED but queued. There is no receptor reference outside Agaricales, so Boletales/Polyporales/Russulales are HD-only. Receptor loci are also called as Balpha/Bbeta (106), aLocus (185) and Tremellales MAT (385); the Basidiomycota:PR family has one record (Coprinopsis B43, Agaricales scope). First step when started: count Basidiomycota genomes with a BFD funannotate proteome; hmmsearch Pfam STE3 (PF02076); build a receptor tree with Coprinopsis/Schizophyllum mating receptors. If they form a clean clade, build a receptor classifier HMM; if not, curate deposits. The hard part is non-mating STE3 paralogs, not sensitivity.
- Per-genome timeout becomes GENOME_TIMEOUT (default 3600 s) in run_clade_panel.slurm. Size bins with a 750 Mb cutoff (curator changed it from 500): <750 Mb on short 2 h jobs; >=750 Mb on epyc with GENOME_TIMEOUT=14400. Implement in PR #9 AFTER the HMM-classifier fork finishes.

### 2026-09-27 Pucciniomycotina
- 11 records signed off: Sporidiobolales redPR (7, tier 1, Coelho 2010/2011) and Pucciniales rustHD (4, tier 2, Cuomo 2017). Held-out: Sporidiobolales 4 -> 15/15, Pucciniales 0 -> 7/8. Branch curation-puccinio is pushed.
- Sporidiobolales stay receptor-only; Sporidiobolales HD is queued for curation. Test (c): extend rustHD scope to Sporidiobolales (results/2026-09-27_puccinio_followup/).
- Wallemia: tier-2 record labelled PUTATIVE from the W. mellicola CBS 633.66 annotation. Synteny and mating-type mix are reported, not used as gates.
- Rust receptors (STE3.2/STE3.3) join the receptor queue. The tiny pheromone precursors need a dedicated plan; start from docs/notes/2026-09-20_short-pheromone-orf-detection.md.
- Fallback confidence floor: target 40%, set from data (counts at 30/35/40% plus a synteny sample; results/2026-09-27_fallback_confidence_floor/). The rule goes into PR #9 after the HMM classifier fork finishes.
- Fallback confidence: option (a), no identity floor (identity does not separate good from bad fallback calls; results/2026-09-27_fallback_confidence_floor/). Microbotryales added to the Pucciniomycotina curation queue.
- Rhodotorula MAT (P/R + HD) from the curator's own preprint (Liu...Stajich, bioRxiv 10.1101/2025.09.11.675505 v2): PAUSED. The curator will provide the manuscript text. bioRxiv blocks automated download. Local data: shared/projects/Rhodotorula/{MAT,MAT_search} and ctsai085/projects/Rhodotorula_comparative_genomics/ (HD1_HD2_tree, STE3_tree, synteny_map). The agent was stopped before any commit.
- Wallemia (2026-09-27): add a second tier-2 putative record from a genome carrying the other version (e.g. W. mellicola EXF-1262). The SXI1 HD1/HD2 class is deferred; it is too tenuous without a tree. Sign-off of the first Wallemia record is not yet stated.

### 2026-09-27 afternoon and later
- Wallemia v1 and v2 records signed off; both pushed on curation-puccinio. Publication highlight logged.
- Tier rule: variant B-prime with the closeness guard is being implemented. An allele-absent gene is ignored only when it is <50% identity AND >=10 points below the called allele's best modelled core gene. Unguarded variant B is unsafe: 108 strong-gene promotions, including S. cerevisiae cassettes and Rhizopus btbA. Replay: results/2026-09-27_tier_rule_replay/.
- btbA's present_in_idiomorphs [Plus] is probably wrong (66-99% identity in 29 Rhizopus Minus loci). Review is PENDING; not changed.
- 11 Rhodotorula records signed off (6 redPR + 5 new redHD, tier 2, preprint source). P/R agrees with the study for 61/62 held-out strains; HD 0 -> 217. Pushed. The unpublished manuscript sits at resource/MBE_202608 (git-excluded; never commit its content).
- redHD cost (5.9x runtime): exploring an HMM prescreen of HD candidate clusters before modelling (results/2026-09-27_hd_prescreen/).
- 2026-09-27: btbA is present in BOTH idiomorphs (restriction removed; still optional and non-voting). homothallic_candidate calls are excluded from the allele-absent confidence rule. The tier rule landed as 1e69622 with a whole-call guard; the real-code replay gave 86 risers.
- HD prescreen: not worth building for speed. The HD cost is mostly the genome-wide tblastn (+21 s of +28 s), not modelling. The untested lever is shorter homeodomain-only redHD queries or a stricter e-value (queued idea).
- Receptor queue order (2026-09-27): (b) pheromone-precursor and positional receptor rule, read-only (results/2026-09-27_pheromone_positional/); then (c) curate receptor/pheromone deposits against that rule; (a) the Sporidiobolales A1/A2 receptor classifier is separate and later. Step 1 found that Agaricomycete mating receptors do NOT separate by tree (3-7 STE3 copies per genome), while Sporidiobolales A1/A2 do (results/2026-09-27_receptor_explore/).
- Receptor curation (2026-09-27): Heterobasidion, Trametes and Grifola B/PR records signed off (PR scope widened to Polyporales and Russulales; Polyporales 0 -> 9/20). The 4.1x Polyporales cost is accepted for now. The CAAX positional filter is a FUTURE test (curator: "CAAX is pretty strong signal"). A tier-2 Russulaceae record is being built (results/2026-09-27_russulaceae_receptor/). Known limitation: reports keep one entry per gene name, so receptor copy number at a locus is not visible.
- Curator (author) allows naming the putative Rhodotorula hybrid strains CCT 0783 (GCA_016808315.1) and RIT389 (GCA_002250355.1) in repo notes and findings (2026-09-27). Other unpublished manuscript content (text, tables, per-strain assignments) stays out of the repo.
- Russula nobilis PR record signed off, but not yet callable: tandem receptors share one gene name, and the 25-aa precursor gets no tblastn hit. The fix is option (ii), a strict-CAAX precursor finder (+-10 kb of receptor hits, per-family roster setting, on for Basidiomycota:PR). It is being implemented; the curator asked for a broad Agaricales test (~100 genomes across families) to confirm it works generally (results/2026-09-27_caax_precursor/).
- DISCOVERY-ONLY lineages (curator, 2026-09-27): Mortierellomycota, Kickxellomycota, Zoopagomycota, Glomeromycotina/Glomeromycota, Blastocladiomycota, Chytridiomycota, and Endogonales (no deposit). They are not curation targets; any MAT work there is a discovery project.
- Mucoromycotina button-up (2026-09-27): trimmed final ML tree (IQ-TREE + RAxML-NG, ~200 tips) running; Umbelopsis regression found (13/14 called at 634dda4 vs 8/14 at 4457c3a, the lost calls look like the Minus ones) being bisected; the 12 uncalled R. arrhizus are being diagnosed. The curator has another set of unscanned Mucoromycotina genomes for testing when ready (path not yet given).
- 2026-09-27: flank-carried rule floor becomes bitscore >= 39 (all families; real 7/8, known noise 0/52; results/2026-09-27_flank_bitscore_floor/). When no core protein is modelled, the classifier scores the tblastn HSP fragment (min_margin kept). Both are being implemented on branch flank-bitscore off PR #9, then fast-forwarded to PR #9. Also check the 12 newly passing Ascomycota calls for synteny.
- The Umbelopsis loss (5 calls) was caused by the genome-size-dependent E floor, not the window. The R. arrhizus GL-series (12 uncalled) have sexP present at 98.4% but split onto its own contig: a split-locus rule is pending a curator decision.
- 2026-09-27: build tier-2 Umbelopsis records for both idiomorphs, to evaluate recoverability (branch curation-umbelopsis). The split-locus rule is APPROVED: a single modelled core gene >=95% near a contig end, with family flanks on other contigs, becomes partial_locus, low confidence, flagged split. Start it after flank-bitscore merges into PR #9 (same pipeline code).
- 2026-09-27 evening: the Lichtheimiaceae/Syncephalastraceae rule is REVISED. A tier-2 record needs a classifier margin >25 bits plus the Mucorales gene order; the tree is reported alongside, and UFBoot>=95 is dropped because the 69-column HMG domain gives even the references only UFBoot 62-74. Merge overlapping same-locus calls across families (Schizophyllum PR + Balpha/Bbeta). Umbelopsis: rebase plus a guard against 13 spurious second calls. The split-locus rule landed (44567fc), recovering 11 GL R. arrhizus + R. microsporus 56028; GL5 GCA_011764265.1 has NO sexP hit (the earlier "all 12" claim was wrong). Planned: a Fable-model review of all PR #9 rule changes once the Umbelopsis/Lichtheimiaceae and merge work lands. Still open: Agaricales families with low CAAX enrichment.
- 2026-09-27: CAAX-only PR calls in low-enrichment Agaricales families (Agrocybaceae, Mycenaceae, Physalacriaceae, Galerinaceae) are labelled verification: unverified (option c; calls kept, confidence unchanged). Added to the locus-merge branch.
- 2026-09-28 Fable review (results/2026-09-28_fable_review/README.md): 6 major findings. F1 (relaxed-pass gate before withhold rules; verified at pipeline.py:2690 vs 2758) explains the R. pusillus FCH_5_7 and A. glauca losses. Zygo 23 is saturated and no longer discriminates, so new validation sets are needed. Awaiting the curator's fix priorities.
- 2026-09-28: fix order approved (F1, F2, F5, F7, then F6, then validation F3/F4 with CAAX-only calls labelled unverified meanwhile, then Q5/F9; the curation-umbelopsis sign-off waits until after F1/F2/F5). Group A merge: option (b) is PROVISIONAL. The curator is unsure whether the biology supports separating Aalpha/Abeta (and Balpha/Bbeta) subloci and thinks they may be a different type of duplication. It is behind a roster switch pending the literature review (results/2026-09-28_subloci_literature/). Group B is unchanged.
- 2026-09-28 SUBLOCI RULING (supersedes the provisional group-A option b): following the literature review (results/2026-09-28_subloci_literature/NOTE.md), subloci (Aα/Aβ, Bα/Bβ, C. cinerea groups) are paralogous specificity units inside ONE locus. Report one A call and one B call per haplotype, with the subloci as structured evidence inside. Group A merges HD+Aalpha+Abeta where they overlap. QUEUED: group subloci by conserved flanks (mip/beta-fg), not by distance; S. commune Aα–Aβ sit ~450-550 kb apart, so the overlap-based merge never joins them. Curator: "see how this performs".
- 2026-09-28 F4 ruling: label ALL CAAX-dependent calls verification: unverified (calls kept, confidence unchanged). REVIEW THIS LATER, once a labelled set of >=100 non-mating STE3 loci exists. Measured: 6/9 mating, 0/25 non-mating (upper limit 13.7%), ~13 expected chance admissions of 118 (results/2026-09-28_validation_f3_f4/). Implement after the review fixes land; also write the review-later item into the handoff queue.
- F3 finding: classifier typing is reliable (0/540 wrong) but the margin does NOT separate MAT genes from HMG paralogs; full-protein absolute score >=100 bits separates well (96/108 vs 9/189), fragments do not. The MAT-vs-paralog test is pending a curator decision.
- 2026-09-28 Q2 ruling (a): MAT-vs-paralog gate for classifier families (Mucoromycota). A modelled core protein needs >=100 bits absolute; below that, or for fragment-typed calls, flank support is required. min_margin 25 stays for typing only. Built with F6, the F4 all-CAAX unverified label and a record assembly-accession field on branch next-fixes. Review fixes landed at 3aec88b (F1, F2, F5 general part, F7, subloci evidence, relaxed pass now uses the classifier).
- 2026-09-28: MAT-gene gate landed (PR #9 076afe4; 900 tests). Lichtheimiaceae: option (a) for now; they are called only at >=100 bits because they lack Mucorales flanks by biology. Lichtheimiaceae MAT gene order and flanks are an EXPLORATION/RESEARCH FOLLOW-UP (add to the handoff queue). The S. racemosum NRRL 2496 record now calls its own locus. curation-umbelopsis is being rebased onto 076afe4, measuring first whether the guard is still needed under the gate.

### 2026-09-28 (held-out sets)
Evidence: `results/2026-09-28_mucor_jena_holdout/`, `results/2026-09-28_lcg_holdout/`.
- Mucor_Jena (65 strains) and ZyGoLife LCG (897 genomes) are held-out test
  sets: never used for training, curation or classifier builds. Results are
  keyed by strain; the curator supplies taxonomy and mating types afterwards.
- Zygo 23 is a subset of LCG and stays the known-answer subset.

### 2026-09-29
Evidence: `results/2026-09-29_*/NOTE.md`, `results/2026-09-28_umbelopsis_rebased/NOTE.md`.
- curation-umbelopsis: guard dropped (it withheld nothing under the MAT-gene
  gate); records 41833 (Umbelopsis Plus), 44442 (Umbelopsis Minus) and 13706
  (S. racemosum NRRL 2496 Plus) signed off; Circinella minor traced (lost to the
  classifier rebuild; consistent with the Lichtheimiaceae ruling (a)).
- The MAT-gene gate threshold comes from each classifier build (manifest).
- Strain-name labels count as known answers; convention Plus/Minus ("+",
  "(+)", "plus" and unambiguous "P" map to Plus; likewise Minus).
- A. blakesleeana is a known gap (option c). Core-only admission of a
  best-hit sexP/sexM cluster (option b) and the rnhA-next-to-core flank
  hypothesis go to the Lichtheimiaceae exploration.
- Two-idiomorph genomes are not confirmed homothallics (could be hybrid or
  fusion, duplication, mixed culture or heterokaryon, or assembly artefact):
  test first; neutral report field built.
- Disputed labels (origin unknown), excluded from scoring: Ellisomyces RSA_581-,
  Gilbertella CBS_442.64-, Pirella RSA_622-, Circinella angarensis RSA_198_Plus,
  C. umbellata RSA_505_Plus, Thamnostylum repens RSA_459_Plus, Backusella
  lamprospora NRRL_6044_Plus. Misidentified, excluded from species-level
  scoring: B. ctenidia NRRL 6239, R. arrhizus NRRL 1470, T. repens NRRL 6240
  (and the M. indicus lineage).
- Build R4, a P1 non-MAT HMG paralog class (trained on one sequence).
- Hold the strong-core floor rescue.
- Record the Syzygites two-idiomorph finding (notable finding 023).
- Classifier builds: deterministic, with `--gate-only` / `--paralogs-only`
  fast paths (option c); no full rebuild of a shipped classifier without
  curator approval.
- Aligner: no aligner more accurate; keep MAFFT with L-INS-i set explicitly;
  no trimming; deterministic full rebuild replayed on Mucoromycota, Zygo, LCG
  and Jena before approval. The aligner finding is recorded in
  `2026-09-29_classifier-builds.md`.
- Future classifiers (Sporidiobolales A1/A2 first; then Ascomycota
  MAT1-1-1/MAT1-2-1 and Serinales MTLa/alpha, gated by training diversity;
  Basidiomycota HD later) use the chosen aligner from the start. Pfam models are
  never rebuilt.


## 2026-09-29/30 (J. Stajich)
- Ship the deterministic explicit-L-INS-i Mucoromycota classifier (approved
  after a 0-change replay on 978 genomes and Zygo 23); shipped 9c39e2c.
- curation-umbelopsis: one deterministic full rebuild with its records; merge
  only if every loss falls in the Lichtheimiaceae gap or in a misidentified
  genome. Rule met; merged a2fe1b4. R. microsporus NRRL A-17693 recorded as a
  misidentified C. minor.
- Polish cap: test protecting strong-fragment clusters by replay; one case was
  judged insufficient, so a stress test was run (caps 3/2 vs cap-off). V3
  adopted (rank strong clusters first inside the cap; never add work).
- Adopt a pre-sign-off regression check: every new or changed record, classifier
  rebuild, paralog class, scope or rule change attaches a diff of every changed
  call, gene model, label, confidence and withheld reason.

## 2026-10-01 (J. Stajich)
- V3 signed off on its regression summary (1 call gained, 0 lost; Zygo 23/23);
  merged 52b3ff9.
- Add three record source genomes to the regression panel (now 166).
- Held-out tables: name both `curator_table.tsv`; key taxonomy to the species
  binomial (lineage filled from NCBI); "T"/"T_of_X" = type strain (of synonym X).
- Apply the curator's species/CBS table to Jena and to LCG by collection number;
  keep file names elsewhere. Ellisomyces NRRL 2465 Plus-T = CBS 243.57 (type),
  putatively Plus.
- No curator mating truth for Jena or LCG beyond Plus/Minus in file names.
- Rerun both held-out sets on current code.
- Runtime: check it, but accuracy comes first; only call-neutral optimisations
  (verified with the regression check) are acceptable.
- Lichtheimiaceae: discovery only; the published Lichtheimia SexM is logged as
  notable finding 024. Circinella labels: tree-based treatment accepted
  (`results/2026-10-01_circinella_label_tree/label_treatment.tsv`).
- Circinella group: drop the rnhA-only gate rule (94d9d37). Sign off
  101103_nrrl1351_MAT_Plus and 64656_rsa-1403_MAT_Minus with caveats
  (04edc2d); accept the Phascolomyces RSA 2281 and Absidia sp. NRRL 3163 Plus
  calls at medium confidence; do not record Absidia NRRL 3163 as
  misidentified (one MAT gene tree is too weak); record A-17693 = C. minor in
  misidentified_strains.tsv and annotation-errors C7. Merged into PR #9
  (5b07b7c). Tag ploidy tests for the two-locus Circinella strains.
- Jena unnamed strains: look up the CBS catalogue (form "CBS 169.57");
  the curator confirms names before use.
- Misidentified strains: keep identities in a versioned MATPredict override
  file (db/curation/taxon_overrides.tsv) with basis and status; never edit
  BFD samples.csv.
- Ascomycota locus names: assistant drafts a synonym map; the curator edits
  it before any code uses it.
- Record self-check: search NCBI for a same-strain assembly; with none, report
  "self-check not possible", not a failure.
- Add GENOME_TIMEOUT with size bins (>500 Mb on epyc, 4 h; <500 Mb on short in
  ~1-1.5 h jobs).
- LCG names: check with marker genes, Mucorales first; flags go to the
  curator, no automatic renames.
- Basidiomycota: no full re-run. Merge curation-puccinio and basidio-anchors
  into PR #9 first (regression check and sign-off each), then a cap-off test
  on ~50 uncalled fallback-order genomes.
- B13b regression: wrong-lineage redPR/wallMAT calls in fallback genomes ->
  scope-only families (`fallback_searchable: false` on redPR, redHD, rustHD,
  wallMAT); PR scope adds Boletales. Net result signed off; curation-puccinio,
  basidio-anchors and the fix merged into PR #9 (5885a7e).
- Done 2026-10-01: Jena names applied from NCBI BioSample (curator approved);
  db/taxon_overrides.tsv (4 strains; T. repens NRRL 6240 newly excluded from
  LCG species scoring); GENOME_TIMEOUT + size-binned submitter; record source
  assemblies for 61 records (23 sequence link, 38 verified strain match).
- M. sympodialis ATCC 42132 bLocus record: use the RefSeq assembly GCF_000349305.1 (curator ruling 2026-10-01).
