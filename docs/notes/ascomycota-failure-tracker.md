# Ascomycota failure tracker

One row per failure type, with what is established, what is not, and the next step. Update the status when a row
moves; add the evidence file. Started 2026-10-17. "Reported by an agent" means a subagent's finding that has not been
independently re-run here.

| id | failure | evidence | status | next |
|---|---|---|---|---|
| F1 | Parastagonospora labelled MAT1-1 in v0.6.0 (87 genomes) | root cause confirmed in the reports: the HMG gene hit the F. graminearum record's `MAT1-1-3`, `resolve_idiomorph_overlaps` chose it over MAT1-2-1 on 41.07% against 40.0% identity, so MAT1-2-1 was superseded. Records run: 91 MAT1-1 and 91 MAT1-2, all agree with the profiles; independent truth (Syme 2013, Bennett 2003) agrees | fixed for this genus by records; software fix approved 2026-10-17 | draft PR on branch `optional-gene-policy` (policy `core_first`, default unchanged) and a paired legacy/core_first run on the 8,938 genomes the fix can change plus a 500-genome control; PR #86 has the replay |
| F1b | the 4 Parastagonospora genomes with no call | tBLASTn finds the gene in 3 of 4 (reported by an agent: two MAT1-2-1 at 99%, one MAT1-1-1 at 86%); the fourth is an assembly gap | approved 2026-10-17: trace | traced: three gates in series (evidence floor, modelled-gene bar, fraction floor); relaxing the floor alone recovers nothing; proposal is a 'strong single core gene' rule (docs/notes/2026-10-17_parastagonospora-misses.md on `parastagonospora-misses`); the earlier check note is docs/notes/2026-10-17_parastagonospora-check.md on this branch |
| F2 | optional gene wins the vote over a core hit, other genera | the same pattern in v0.6.0: Pyrenophora 59, Corynespora 10, Pyricularia 197, Fusarium 433, Calonectria 78, Verticillium 29 (reported by an agent, not examined) | open | measure with the replay (variant C in PR #86) and check a sample against profile calls |
| F3 | 17 Botrytis genomes stated MAT1-2, called MAT1-1 high | reported by an agent: the 'MAT1-2-1 at 77%' hit is the true B. cinerea MAT1-2-1 (98 to 100% to KF944386.1); the MAT1-1-1 'hit' is a 61% remnant that matches the curated MAT1-1 record at 91 to 95% and wins the bitscore vote; the database has no B. cinerea MAT1-2 record; 32 of 32 PRJNA624742 genomes agree with their BioSample label by gene content; not a mixed isolate (one genome checked at read level) | software call wrong (reference gap), BioSample label right | record curated on branch `botrytis-curation` (commit 4e40985) but it does NOT flip the 17: all still MAT1-1 (the MAT1-2 score rose from about 222 to 326, the 92.7%-identity MAT1-1-1 remnant still scores 416); the stated-MAT1-1 genomes are unchanged; a vote change in `src/` is needed (coverage weighting or an alpha/HMG profile tie-break), not yet approved (notes: docs/notes/2026-10-17_botrytis-mislabels.md on this branch, docs/notes/2026-10-17_botrytis-curation.md on `botrytis-curation`) |
| F4 | Dipodascales: 67% of good-assembly genomes uncalled, fallback routing | reported by an agent: Yarrowia MATB not curated (28 of 30 uncalled Yarrowia carry it, CLIB122 at 100%); Geotrichum MATA/MATB deposits (HF558448/9) match at 95 to 100% beside SLA2; 34 of 113 called genomes are called more than 30 kb from SLA2; Wickerhamiella and Starmerella have no reference | decisions 2026-10-17: extend `MATyl` to the Yarrowia genus (4951), curate Yarrowia MATB, add Geotrichum, address the wrong-contig calls; remaining questions go to the curator | an agent is curating on branch `dipodascales-curation` (records, `order.yml`, before/after on the wrong-contig genomes); notes on branch `dipodascales-refs` |
| F5 | Microascales (Ceratocystis fimbriata) uncalled; 11 genomes recoverable from the deposits at about 92% identity | rescue test (2026-10-05) | decisions made 2026-10-17: add a `unisexual_switching` system value to the schema; verify on a same-strain assembly before the records are accepted; hold the 11 genomes out as the test | an agent is on branch `ceratocystis-curation`: schema change with tests, then the assembly search, then the records if an assembly qualifies |
| F6 | Xylariales 71% uncalled | known absence (Robinson and Natvig 2019: no canonical MAT in 35 Xylariales); synteny work found an empty or rearranged flank block in NC1011 | expected biology, now shown: 0 of 184 uncalled Xylariales genomes carry a MAT-like domain next to the flanks (analysis/2026-10-18_uncalled-triage.md); 11 of 12 called Xylariales controls also show none, so some Xylariales calls may rest on a non-MAT HMG protein (unchecked) | inspect the called Xylariales genomes |
| F7 | held-out misses | Taphrina deformans `PM` missed at every radius; S. cerevisiae HMRa and K. lactis MATa called as the active locus; Serinales `MTL` 1 of 6 and Mucorales 2 of 7 at order radius; N. crassa MAT1-1-1 only a short fragment with its own records excluded | open | each needs a reference or a scoring decision |
| F8 | PR hygiene | #71 (tracks, GenBank input) and #76 (profiles) were stacked and merged into branches that had already gone to main; #71's content is still not on main; #76 is re-opened as #86 | #71 content re-opened as #88; #76 re-opened as #86; branches `botrytis-mislabels` and `parastagonospora-check` on the remote are empty: the two earlier agents committed their notes in the main checkout, so those commits landed on `uncalled-triage` (7c166fc, 453c0d4, 2363a9f) | merge #86 and #88; decide what to do with the empty branches |

## Work started 2026-10-18 (curator confirmed the decisions; every branch is a draft or notes, nothing merged)

| work | branch | decisions it follows |
|---|---|---|
| single strong core gene rule | `single-gene-rule` (draft PR) | `idiomorph_gene_only` at low confidence, identity gate 85, off by default, two opposite loci both reported and flagged through `two_idiomorphs_report`; split-locus limits unchanged |
| Komagataella MATa (CBS 7435) and MATalpha (GS115) records | `komagataella-curation` | same-strain assembly anchored; the eight Komagataella test genomes stay out of the records |
| P. nodorum flank records | `pnodorum-flanks` | split-locus limits unchanged |
| polish-cap fix | `polish-cap-fix` (draft PR) | its own PR and campaign; does not hold the Dipodascales guard |
| review of the 45 unverified flips; margin calibration | `flip-review` | default stays off until done; both-strong loci stay `undetermined`; the four genera checked here |
| general confidence cap (partial winner) plus weak-identity criterion | `confidence-cap` (draft PR) | separate follow-up after #95, with a campaign |
| remnant status paper | `remnant-status` | revisit after #95 (merged) |
| Ceratocystis: decisions, re-validate, recall test, draft PR | `ceratocystis-curation` | switched record stays `needs_review`; neutral labels; the 11 held-out genomes no longer test recall |
| optional-gene policy and paired campaign | `optional-gene-policy` (PR #89, draft) | waiting on the campaign; Sl8 interaction reported there |
| Dipodascales: guard, scope-only families, records, Pichiales deep dive | `dipodascales-curation` | held until the single-gene rule is built and evaluated; no idiomorph_class for MATdc until a Dipodascales profile exists; Orbiliales losses accepted |
| deferred | | trace of the other uncalled Pleosporales (after the single-gene rule is evaluated); Dipodascales-specific profile; Lipomycetales/Blastobotrys/Sugiyamaella families (52 and 20 genomes lose calls to the guard) |

## If the session ends (state as of 2026-10-18, late)

Cluster jobs keep running without the session; agents may not. Work in flight and where it lands:

| item | cluster state | output location |
|---|---|---|
| optional-gene paired campaign (PR #89, draft) | array 30031368, 557 of 941 tasks complete, 58 running | `/bigdata/stajichlab/jstajich/projects/optgene_runs` (core_first and legacy) |
| re-run of the 63 linked triage genomes on current main | array 30034651, queued | `/bigdata/stajichlab/jstajich/projects/main-now/out/<genome>/detection_report.yaml`; list in `main-now/linked63.tsv` |
| single-gene rule | `sgr-eval` running | branch `single-gene-rule`, worktree `MATPredict-sgr` |
| polish-cap fix | `cap-sample`, `cap-panel` running | branch `polish-cap-fix`, worktree `MATPredict-cap` |
| general confidence cap | branch has code | branch `confidence-cap`, worktree `MATPredict-conf` |
| Komagataella and P. nodorum records | `kom-tests` running | branches `komagataella-curation`, `pnodorum-flanks` |
| remnant status step 1 and step 2 | `rem-s1` running | branches `remnant-status-step1`, `remnant-status-step2` |
| curate-db validate wrapper (2 GB) | done | PR #99 |
| Ceratocystis | done | draft PR #98 (merge recommended by its agent) |
| Dipodascales guard, records, Pichiales deep dive | done | branch `dipodascales-curation` (no PR; held until the single-gene rule is evaluated) |
| flip review | done | branch `flip-review` (no PR) |

Every branch pushed so far holds its work up to its last commit. Uncommitted work in a worktree is not on the remote.

## Revisit list (decisions deferred on purpose, 2026-10-18)

Each item says what was decided now, what triggers the revisit, and what data to bring.

| item | decided now | revisit when | bring |
|---|---|---|---|
| partial-core rule profile margin | stays 25 bits (borrowed from Mucoromycota, never calibrated) | the companion-gene veto's results are in (branch `companion-gene-veto`) | the flip-review margins: only two loci fall at 25-70 bits (Monilinia 40.5 wrong, Stromatinia 66 right); S. nivalis is wrong at 83 while Stromatinia (66) and Sporothrix (89) are right; the reviewer proposed 50 with about 40% confidence from four wrong calls in one clade; check what the veto leaves at 25-70 bits and test more clades (docs/notes/2026-10-17_flip-review.md on `flip-review`) |
| partial-core rule default | OFF | the veto is in and the unverified flips are re-checked | the flip review: 41 of 46 flips right; safe now only as a clade allow-list (Botrytis, Botryotinia); Ophiostomatales and Stromatinia probably right once the veto exists; not for Sclerotinia or Monilinia |
| both-strong veto (set-aside >= 200 bits and winner >= 150) | not implemented (thresholds fitted on 49 loci) | more loci with both idiomorphs are available | whether the companion veto already covers them |
| fused or homothallic loci (Sclerotinia, Monilinia, Polytolypa, Calogaya) | stay `undetermined` | the veto is in | whether `full_length_models` should stop excluding MAT1-1-3/MAT1-2-1 cross-matches so `homothallic_candidate` can fire (47 of 49 loci blocked by that exclusion) |
| general confidence cap (`high` to `medium` when the winner is partial) plus weak identity | separate follow-up | `confidence-cap` evaluation is back | how many `high` calls drop and whether the dropped ones are wrong more often than right |
| remnant status | step 1 (metadata) now; step 2 (remnant proteins leave the references) only after a measured rerun | `remnant-status-step2` measurement is back | the Afum `homothallic_candidate` label (does it clear), clades that lose their only reference |
| single strong core gene rule | `idiomorph_gene_only`, low confidence, gate 85, off by default | the control run on called genomes is back | whether any called genome gains a locus; the 85 vs 90 trade-off |
| Dipodascales: `idiomorph_class` for MATdc | left unset | a Dipodascales-specific profile is built from the recovered Geotrichum genomes | HMG-box and alpha-box scores of MATA and MATB against that profile |
| Dipodascales guard | kept; merges after the single-gene rule is evaluated | the single-gene rule is evaluated | 52 Lipomycetales and 20 Dipodascales genomes that lose calls; each clade needs its own family |
| uncalled Pleosporales trace (49 flank-carried, 32 zero-modelled) | deferred | the single-gene rule is evaluated | the genomes the rule does not recover |
| optional-gene policy default | off (legacy) | the paired campaign (array 30031368) finishes | label changes by genus, agreement with truth, the Sl8 interaction with `core_first` |
| Orbiliomycetes MATsc fragments (2 genomes) | guard's removal accepted without inspection | Orbiliomycetes get their own references | GCA_047671915.1 and GCA_047716685.1 |
| B. rouxiae GCF_048999995.1 MAT1-1-1 | explained: a 33.9-bit fragment now named MAT1-2-7 (PR #98) | a second genome shows the same pattern | the locus alignment |
| S. spermophila verdict (flip review: unclear, low confidence, margin 229) | not checked | the companion-gene veto results are in | self-fertility literature and, if needed, reads; a MAT1-2-only Sclerotinia would be unexpected because S. trifoliorum switching loses the MAT1-2 genes (doi:10.1038/srep27083) |
| read-level checks of Monilinia (RL-1) and the Sclerotinia fused structure | skipped: with the companion veto those loci become `undetermined` whatever the reads show | someone wants the biology | miniprot ranges so far, not full gene models |
| A. fumigatus edge-kept contig | folded into the remnant step 2 measurement (branch `remnant-status-step2`) | that measurement is back | the genes' coverage, the reference each model was built from, the effect of step 2 on its label and `homothallic_candidate` |
