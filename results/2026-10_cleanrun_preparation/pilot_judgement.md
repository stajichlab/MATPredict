# Clean-run pilot vs v0.6.0: judgement of the changed calls (2026-10-22)

Inputs: 565 pilot reports (`results/2026-10_cleanrun_pilot/*/runs`, frozen worktree `run-8b1a80d`, tag `cleanrun-2026-10`, pilot flags) against `2026-10-03_{ascomycota,basidiomycota}_v060` (`loci.tsv`, `genomes.tsv`, and the per-genome v0.6.0 reports extracted from `reports_all.tar.zst`). I re-ran `pilot_cmp.py` and got the same counts as the caller. All diagnostics are on the HPCC in `results/2026-10_cleanrun_pilot_diag/`: `pilot_changes.tsv`, `dump_all.txt` (old and new gene evidence for all 214 changed genomes), `split.txt`, `spans.txt`, `lcount.txt`, `confup.txt`, `scripts/`, and the `runs_<config>/` re-runs. Nothing in the repo or the frozen worktree was edited.

## How the causes were established (re-runs on `highclock`, 32 genomes covering every group)

| config | code | db | flags | what it shows |
|---|---|---|---|---|
| `pilot` | frozen | frozen | pilot flags | **32/32 identical to the production pilot**. Determinism holds on this set. |
| `v060` | v0.6.0 (`git archive v0.6.0`) | v0.6.0 | defaults | **32/32 identical to the v0.6.0 campaign**. The baseline reproduces, so the diffs are not run noise. |
| `off` | frozen | frozen | none (all defaults) | Shows the effect of code and db changes with every policy off. |
| `nopc` | frozen | frozen | pilot minus `--partial-core-policy` | Isolates the partial-core rule. |
| `v060code_newdb` | v0.6.0 | frozen | defaults | Shows what the db change alone does. |
| `newcode_v060db` | frozen | v0.6.0 | defaults | Shows what the code change alone does. |

(`nohomo` was cancelled. It was redundant once `off`/`nopc` were in, and the account CPU limit was holding the queue.) The cap and single-gene policies are `off` in the pilot, which is their code default, so they cannot cause any difference. Attribution per group:

- **Flag-driven:** the partial-core rule causes the Botrytis/Sclerotinia flips: `off` and `nopc` give the old MAT1-1 call, and `nopc` already gives the homothallic calls. The homothallic policy causes the 30 `undetermined` loci: `off` gives the old call, `nopc` gives `homothallic_candidate`. The scoped roster causes the Pleosporaceae recoveries: CuVerr and Cocsa1 are absent in `off` and `nopc` and present in the pilot.
- **Database-driven (v0.6.0 code + new db reproduces the change; new code + old db does not):** all Aureobasidium flips and confidence drops, both Xylaria flips, both Parastagonospora flips, the Aspergillus class change, all five confidence DOWN, the Eutypa and Plectosphaerella losses, and all non-scoped NEW calls (Komagataella, Yarrowia, Geotrichum, Ceratocystis, Cordyceps, Pestalotiopsis, Hemileia, Ogataea MATog).
- **Code-driven:** the Phyllozyma HD loss is reproduced by new code + old db, so it comes from the 5 kb intron cap (#92). The Lipomycetales/Dipodascales losses need both the new code and the new db: the subphylum fallback guard is `route()` in code plus `fallback_subphylum` in `order.yml`.
- A side observation: v0.6.0 code + the new db searches GAP1/BGL1 **unscoped**. It "recovers" CuVerr and Cocsa1 that way. So the new db must never be run with old code.

## Per-category table

| category | n | cause | verdict |
|---|---|---|---|
| identical | 259 | | |
| uncalled both | 92 | | |
| **NEW call** | 71 | **40 scoped roster** (Bipolaris 26, Curvularia 7, Exserohilum 7), all Pleosporaceae. **8 genetic-code-26 fix**: Nakazawaea 7 and Peterozyma 1 (Alaninales) had *no v0.6.0 report*; the code-26 exonerate defect is fixed. 8 Komagataella `MATkp`. 3 Yarrowia `MATyl`, 2 Geotrichum `MATdc`. 4 Hemileia `rustHD` (bE/bW 63-100% identity, 98-100% coverage). 4 Ceratocystis (new C. fimbriata records). 1 Cordyceps, 1 Pestalotiopsis (new references). | **expected / improvement.** The non-scoped NEW calls all rest on full-length, high-identity genes: Yarrowia matb1/matb2 68-100%, Hemileia bE/bW, Ceratocystis MAT1-1-1 93% at 100% coverage. Two weak ones: Pestalotiopsis GCA_030412595 (`high` from a 16%-coverage MAT1-2-1 plus SLA2) and Cordyceps (medium, gene-only). Komagataella fragmented assemblies (2.1-2.3 kb contigs) report only the MATa cassette, which is half of a switching genome (expected from the assembly). |
| **idiomorph changed** | 49 | **30 homothallic policy** (Pyrenophora 17, Fusarium 5, Sclerotinia 5, Ceratocystis 2, Monilinia 1). **11 partial-core flips** (Botrytis 10, Sclerotinia spermophila 1). 2 Parastagonospora (new `13684` records). 4 Aureobasidium (new `1873960` Mycosphaerellales record lengthens the MAT1-1-1 model). 2 Xylaria (new references move a 30-50 aa HMG fragment vote). | homothallic: **expected**. Botrytis/Sclerotinia: **improvement**. Parastagonospora: **improvement**. Aureobasidium and Xylaria: **needs attention** (neither old nor new label is supported). |
| span changed >1 kb | 40 | 15 GAP1 added (scoped roster: Pyrenophora teres 8, Alternaria 7; MAT1-1-1/MAT1-2-1 models also lengthened). 18 outer cluster span grows or shrinks from weak hits that are not reported, while `core_span` stays at the old extent (Hysterium, Venturia, Botryosphaeria, Eremomyces, Remersonia, Tilletiaria, Ogataea x2, Saccharomycodes, Wickerhamomyces, Alloascoidea, etc.). 4 a distant weak unpolished hit is added and stretches the span 24-35 kb (Scorias and Colletotrichum get MAT1-1-1 hits in MAT1-2 calls; Pleurotheciella MAT1-1-3; Serendipita pheromone). 1 chimeric 15 kb MAT1-1-1 model fixed (Phaeoisaria, intron cap). 2 comparison artefacts: Aureobasidium GCA_004917145, where the real change is a confidence drop on the second locus, and Schizophyllum, where one locus is relabelled `Aalpha` -> `HD` after a 7.7 kb Z model was dropped. | expected (GAP1, intron cap); cosmetic (outer span); **minor attention** (the 4 distant-hit additions) |
| locus count changed | 18 | 5 Ogataea complex -> `MATog` (expected; in scope taxids 460523/1005962/2053387). 2 Komagataella -> `MATkp`. 7 spurious cross-family duplicates removed by the subphylum guard / scope-only families (`MATyl` in Petasospora and O. henricii; `PM` in Pichia x2; `MAT` in O. piceae, Wickerhamomyces, Blastobotrys). 1 Kurtzmanomyces loses a `high` bLocus whose bE model was a 25.7 kb chimera at 27% identity and 29.6 bits (intron cap). 1 Auriscalpium gains a PR locus. 1 Bertia gains a weak third MAT1-2 partial locus. 1 Coniochaeta endophytica loses a relaxed gene-only locus (`polish_capped: true`: outranked under the 6-cluster polish cap). | expected / improvement except Coniochaeta (cap-ranking noise) and Bertia (weak, unverified) |
| confidence up | 11 | 3 Alternaria (GAP1 polished). 2 Yarrowia + Sporopachydermia (MATyl/MATsc references). Metarhizium and Chalaropsis (a gene now polished). **Pyrenula, Xylaria grammica, Tolypocladium**: medium/low -> high from a distant unpolished fragment (an opposite-idiomorph core hit in Pyrenula and Xylaria grammica). | mostly expected; the last 3 are the issue-#129 tie-break pattern **outside** the Pleosporaceae scope |
| LOST call | 10 | 7 subphylum fallback guard (Lipomyces x3, Myxozyma, Tortispora, Saccharomycopsis, Sugiyamaella): `Ascomycota:MAT`/`MATtub` (Pezizomycotina) are no longer searched in Saccharomycotina-type genomes. Saccharomycopsis' old MAT1-2-1 model was also a 37 kb chimera. 1 Phyllozyma HD (old HD2 model 79 kb, 28.5 bits; intron cap). 2 Eutypa and Plectosphaerella sp. AK1849 (db-driven; old calls were weak relaxed `partial_locus`). | guard: **expected** (plan §6.4). The old calls were cross-family labels, and Lipomyces carried MAT1-1 and MAT1-2 on the *same* locus from two families. Phyllozyma: **improvement**. Eutypa/Plecto: **unverified** |
| genes_found changed only | 9 | 4 Ceratocystis (full-length MAT1-2-7, MAT1-1-1, MAT1-1-2 now modelled; new records). 3 Parastagonospora (an unpolished 58-aa MAT1-2-1 cross-hit added to MAT1-1 calls). Pyricularia, Coniochaeta pulveracea (MAT1-1-3 / MAT1-2-1 HMG cross-hits reshuffled). | expected; Parastagonospora's added opposite core hit is cosmetic |
| confidence DOWN | 5 | Aureobasidium x2: the MAT1-1-1 model grows (53.9 -> 81.6 bits), so the vote margin falls 30.8 -> 3.5 (< `min_idiomorph_margin` 5). **Pyricularia T47, Thyridium**: the opposite hit (MAT1-1-3 or MAT1-2-1) now loses an overlap resolution, so only one candidate is left and the margin **falls back to the identity-point margin of that resolution** (4.6 / 4.9 < 5). The call is unchanged. CTG MAG: MTLA1 drops from polished to unpolished (3 -> 2 polished genes). | Aureobasidium: reflects the real ambiguity. Pyricularia/Thyridium: **needs attention** (a systematic artefact of the margin fallback, not evidence against the call). MAG: minor |
| class changed | 1 | Aspergillus fumigatus GCA_060371155: v0.6.0 `homothallic_candidate` came from the old different-domains path with a 1,035 bp MAT1-2-1 model. With the updated `af293` record, that model is 363 bp (37.6% coverage), so neither the old path (needs >= 50% coverage) nor the both-genes rule (needs >= 190 aa) fires, and the locus becomes `mat_locus` MAT1-1 high. MAT1-1-1 is 100% identical at full length. | **improvement**: A. fumigatus is heterothallic, and this is the false positive the homothallic controls (0/40 A. fumigatus) are meant to exclude. The policy is not applied because the HMG gene is a fragment, which is correct. |

## (a) Idiomorph changes in detail

- **The 30 homothallic candidates are all policy effects.** Every one carries `homothallic_note` with `decision: homothallic_candidate`. In 30/30, `vote_idiomorph` equals the v0.6.0 idiomorph. `previous_locus_class` equals the old class in 28. In the 2 Ptr (Ptr5213, Ptr239) it reads `partial_locus` while v0.6.0 said `idiomorph_gene_only`: the locus grew to a strict-pass call, so the note's class is the new run's pre-policy class.
  - Biology agrees in every group. P. tritici-repentis: three BioSamples state "homothallic", and the measurement note cites Lepoint 2010. F. graminearum and S. sclerotiorum are textbook homothallic. Ceratocystis is an unswitched unisexual-switching locus with both MAT1-1-1 and MAT1-2-1 (Wilken 2014). M. vaccinii-corymbosi was already a fused-locus verdict in the flip review.
  - Inconsistency, not a regression: 5 other Ceratocystis loci with full MAT1-1-1 and MAT1-2-1 models stay `MAT1-2` `mat_locus` / `idiomorph_gene_only` (GCA_003449175, GCA_038407085, GCA_038382965, GCA_038407035, GCA_038381545). Two are relaxed-pass. Why the rest fail was not established, because the report does not record which gate failed.
- **Botrytis 10 MAT1-1 high -> MAT1-2 medium: right direction.** All 10 have the same signature: MAT1-1-1 is the 5'-truncated remnant (91-94% identity, **62% coverage**, 410-421 bits) and MAT1-2-1 is complete (99.7-100% identity, **100% coverage**). The profile margin is 95-96 bits for HMG. Three are in the BioSample truth set as MAT1-2 (Vv11 GCA_015148905, Vv9 GCA_015148935, Sl3 GCA_015149565). The other 7 have no truth but carry exactly the MAT1-2 idiomorph structure of the curated SAS405 record (Abate 2018, dMAT1-1-1 fragment). Every Botrytis genome with a full-length MAT1-1-1 (98.6% coverage) stays MAT1-1 high, which is the right control. The `medium` cap after a flip is the provisional rule. Sclerotinia spermophila GCA_038090695 has the same pattern (MAT1-1-1 66% coverage vs MAT1-2-1 99.7%, profile margin 229 bits) and is judged right, though there is no truth for it.
- **Parastagonospora 2 MAT1-1 medium -> MAT1-2 high: improvement.** The new call rests on a complete MAT1-2-1 (100% identity and coverage, 377 bits, from the new `13684_sn436ga98` record). The old MAT1-1 call was a 45-bit hit; in the new run MAT1-1-1 is a 31-bit unpolished fragment. Its sister genomes with a full MAT1-1-1 stay MAT1-1 high.
- **Aureobasidium 4 MAT1-2 -> MAT1-1: needs attention, neither label is right.** Every Aureobasidium locus carries *both* a MAT1-1-1 and a MAT1-2-1 fragment, each covering 30-40% of its reference at 32-53% identity. The vote margins are 0.4-17 bits. The flips come from the new Mycosphaerellales record `1873960_unknown-1`: the MAT1-1-1 model lengthens and its score rises 54-79 -> 82-105 bits. v0.6.0 code with the new db reproduces every flip. The independent evidence points to homothallism: according to PubMed, Gostinčar et al. 2019 found "a homothallic mating locus ... in all of the sequenced genomes" of 50 A. pullulans strains ([doi:10.1111/1462-2920.14693](https://doi.org/10.1111/1462-2920.14693)). The homothallic rule does not fire because both models are 100-200 aa alignment-edge fragments (the measurement note already reports Aureobasidium at 0 of 119). Both MAT1-1 and MAT1-2 are coin tosses here. The gene-model extension work, which translates to the next stop, is the identified fix.
- **Xylaria 2 (FL1030 MAT1-1 -> MAT1-2; Z50 MAT1-2 -> MAT1-1): needs attention.** In both, the vote is MAT1-1-3 against MAT1-2-1, two overlapping 30-50 aa HMG-box fragments at 9-24% coverage, with margins of 1.2 / 0.4 (Z50) and 3.8 bits. New references shifted the bitscores by a few bits (db-driven). No truth exists; these are not in the Robinson & Natvig table. Neither call is supported; both stay `medium`.

## (b) LOST calls

| genome | old | cause | was v0.6.0 right? |
|---|---|---|---|
| Lipomyces x3, Myxozyma, Tortispora, Sugiyamaella, Saccharomycopsis | Pezizomycotina `MAT`/`MATtub` labels via `phylum_fallback` | subphylum guard (curator 2026-10-17); needs new code + new db | The label was not right. These are Pezizomycotina-family labels on Saccharomycotina-type yeasts, built on 7-23%-coverage hits, and Lipomyces got MAT1-1 (`MAT`) and MAT1-2 (`MATtub`) on one locus. An APN2-SLA2 neighbourhood suggests a real locus exists, so this is a recall loss until those lineages have their own families. Sugiyamaella's `high` MAT1-2 rested on one HMG gene at 44% identity with an unpolished SLA2. |
| Phyllozyma linderae (HD) | HD2 "model" 79 kb, 28.5 bits | intron cap (#92), code-driven | no (chimera) |
| Eutypa lata | MAT1-2-1 at 27% identity, 42 bits, relaxed `partial_locus` | db-driven; the new best cluster has only 1 modelled gene | unverified, weak |
| Plectosphaerella / Xylariales sp. AK1849 | MAT1-2-1 135 bp model, relaxed | db-driven; the new best cluster (MAT1-1-1 + APN2 + COX13) is flank-carried with a core hit under 39 bits | unverified. The new evidence even points to MAT1-1, so the old MAT1-2 label is doubtful. |

## (c) Confidence DOWN and class change

See the table above. Main finding: **the single-candidate margin fallback is fragile.** When the only opposite-idiomorph hit loses an overlap resolution, `idiomorph_margin` switches from a bitscore gap to that resolution's identity-point gap. That gap is between MAT1-2-1 and a MAT1-1-3 cross-hit on the same HMG box. It is then compared with `min_idiomorph_margin` = 5 and drops `high` to `medium` (Pyricularia T47, Thyridium). That is not evidence against the call. The same mechanism also explains the medium tier of the Parastagonospora and Xylaria single-candidate calls. The Aspergillus class question is answered above: the HMG side is now a 121-aa fragment, so the policy correctly does not apply.

## (d) NEW-call spot checks (15)

Bipolaris GCA_008452705 and GCA_035078985 (MAT1-2, BGL1+GAP1), Curvularia CuVerr and GCA_002161795, Exserohilum GCA_019453395 and GCA_024221855, and Cocsa1 (B. sorokiniana ND90Pr): all have a polished core gene of one idiomorph plus the scoped flanks. The opposite idiomorph is at most a 14-21%-coverage unpolished cross-hit. This is consistent with the 90/101 recovery audit, which had a correct idiomorph and location in every case. I did not check the BioSample mating type of the Cochliobolus genomes individually.

The other spot checks:
- Yarrowia lipolytica GCA_044705915 (matb1/matb2 100%) and Y. galli: right.
- Hemileia x2: bE/bW 100% identity.
- Ceratocystis platani: MAT1-1-1 93% and MAT1-1-2 99%, full length.
- Geotrichum x2: MATB 53-56% plus sla2.
- Komagataella GCA_023013795: MATa1/MATa2/sla2/dic1.
- Nakazawaea GCA_024271865: MATALPHA1/2 + sla2. Code-26 genomes, never run in v0.6.0.
- Pestalotiopsis GCA_030412595: `high` on a 16%-coverage MAT1-2-1 (51.6 bits) plus SLA2. This one is weak and is the tie-break pattern again.

None contradicts independent evidence. Only the Botrytis/Ptr truth and records were checked; there is no BioSample truth for the others.

## (e) Spans and locus counts

Covered in the table. Two points for the site:
1. The reported `start`/`end` is the raw cluster span. In 18 genomes it grew or shrank by up to 30 kb from weak hits that are not reported, while `core_span` (new in this run) kept the old extent. Site maps should draw `core_span` or `supported_span`.
2. In 7 genomes (Scorias, Colletotrichum, Pleurotheciella, Serendipita, Xylaria grammica, Tolypocladium, Pyrenula) a distant (20-45 kb) unpolished fragment now enters `genes_found`. In Pyrenula, Xylaria grammica and Tolypocladium it lifts the tier.

## (f) Supported vs tie-break, Pleosporaceae `high` calls (pilot)

Extra gene = any gene_evidence entry except the called core gene and GAP1/BGL1. Taking the best entry per gene:
- **supported** = at least one non-opposite extra gene that is polished or reaches >= 50% identity at >= 30% coverage;
- **tie-break** = every extra gene is unpolished with < 30% coverage, or is the opposite core gene;
- anything else = **grey**.

| | high | medium |
|---|---|---|
| supported | **18** (Alternaria 10, P. teres 8) | 15 |
| tie-break | **6** (Curvularia 3, Exserohilum 3; all NEW calls) | 18 (+15 with no extra gene at all) |
| grey | **1** (Curvularia BF9436 GCA_965118015) | 2 (Bipolaris) |

The thresholds separate the groups cleanly. Every supported call has polished APN2 and COX13 at 51-56% identity and 76-100% coverage. The largest non-opposite extra gene in a tie-break call is 27% coverage (unpolished MAT1-1-3 or APN2 at 5-6%). Nothing falls between 30% and 75% coverage, so the result is not sensitive to the exact cut.

The **only grey cases** are a MAT1-1-3 hit that *lost an overlap resolution to the called MAT1-2-1* (unpolished, 24-31% identity, 39-49% coverage): an HMG-box cross-hit on the core gene itself.

Proposal: treat any extra gene with `lost_to` equal to the called core gene as an opposite-idiomorph core hit, which makes it tie-break. With that rule the split is **18 supported / 7 tie-break / 0 grey**. All 7 tie-break `high` calls are new scoped-roster recoveries, as issue #129 predicted.

Note that the same pattern also produces `high` calls outside the scope: Pyrenula, Xylaria grammica, Pestalotiopsis GCA_030412595. It may be worth computing the split campaign-wide, not only in the Pleosporaceae.

## (g) Mucoromycotina (25)

24 identical to `2026-10-03_mucoromycotina_mat/calls.tsv` in idiomorph, confidence, class, contig, span and genes, including the one uncalled genome (Jena CBS564_66). The remaining one (Gonbut1 GCA_025201375) shifts its start by 3 bp with the same call. As the plan says, no policy acts on this phylum.

## Needs attention

1. **Aureobasidium (8 changed genomes: 4 idiomorph, 2 confidence, 2 span)**: both-gene fragments, margins 0.4-17 bits, homothallic per Gostinčar 2019. MAT1-1/MAT1-2 should not be published as calls. Report them `undetermined` or flag them until the gene-model extension lets the homothallic rule see them.
2. **Margin fallback** (Pyricularia GCA_059330725, Thyridium GCF_004353045): a confidence drop caused by an identity margin against a same-HMG cross-hit. It is systematic and will recur campaign-wide; the count is not measured.
3. **Weak distant / opposite-idiomorph fragments inflating `genes_found` and tier** (Pyrenula GCA_051942085, Xylaria grammica GCA_022495055, Tolypocladium GCA_024340795, Pestalotiopsis GCA_030412595, Scorias, Colletotrichum, Pleurotheciella, Parastagonospora x3): this is the #129 mechanism outside the scope.
4. **Xylaria FL1030 / Z50**: HMG-fragment votes within 0.4-6 bits. They should be `undetermined`, not a label.
5. **Unverified losses**: Eutypa lata GCA_022661275 and Plectosphaerella GCA_022432505 (db-driven, weak old calls); Coniochaeta endophytica GCA_007388115 (polish-cap ranking).
6. **Ceratocystis**: 5 both-gene loci are not classed `homothallic_candidate` while 2 others are.
7. **Provenance**: every pilot report says `matpredict_version: 0.6.1`, but the code is at 0.6.3 per the plan. The version string comes from the installed package metadata, not the frozen source. The run manifest must carry the commit and tag; do not rely on this field.

## Recommendation

**The full run can proceed.**
- Every changed category is attributed to a specific approved change: a policy flag, a db record or the guard, or a code fix (intron cap, gencode 26).
- The baseline and the pilot both reproduce exactly on 32 re-runs, so the diffs are not noise.
- The call rate rises (402 -> 463 called genomes). The only losses are the guard losses the plan expects, plus three weak or chimeric calls.
- No change was found that moves a call against independent evidence. Every Botrytis truth case is fixed, and the Ptr, Fusarium and Sclerotinia homothallic candidates agree with the literature.

Items 1-4 do not change what the run computes. They are reporting and interpretation decisions for the post-run diff and the site, and they can be applied to the reports afterwards:
- an Aureobasidium/HMG-fragment `undetermined` rule;
- the margin fallback;
- the supported/tie-break split extended campaign-wide, with the `lost_to` amendment;
- `core_span` for maps.

If the curator wants item 2 or item 3 fixed in the tier itself, that is a code change and must happen before the run.

## Not verified

- No reads or BioSample truth beyond the 3 Botrytis and 3 Ptr genomes. The 7 other Botrytis flips, S. spermophila, the Cochliobolus recoveries and all NEW yeast calls are judged on gene evidence and records only.
- The Aureobasidium conclusion rests on one paper's abstract (read via PubMed), not the full text.
- Why the 5 Ceratocystis loci fail the homothallic gates.
- What exactly sets the outer span in the 18 span-only cases: weak unreported hits in the cluster, shown in one case (Hysterium, cluster of 88 hits), not traced in every genome.
- The campaign-wide frequency of the margin-fallback and distant-fragment effects.
- `nohomo` was not run; homothallic attribution rests on `off` vs `nopc` vs pilot.
