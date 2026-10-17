# Clades that need their own curation (2026-10-18)

A working list for the curator, assembled from notes and analyses already in the repository. **Nothing here was re-run.**
Counts are the ones those documents give, with the source named in each row; an item marked *(agent)* is a subagent's
finding that was not re-checked. The campaign they come from is v0.6.0 (2026-10-03) unless a row says otherwise, so
"uncalled" in particular may be out of date on current main; the planned clean full-BFD run will refresh it.

Priority below is **my judgement** (genomes affected, how clear the evidence is, how much curation it takes), not a
measurement. Status is as of this date.

## How to read the kinds of gap

| kind | meaning | typical action |
|---|---|---|
| **A. No family** | no family covers the clade, so its genomes are called by a distant family or not at all | curate records and declare a lineage-scoped family |
| **B. Thin references** | a family exists but the clade has one or no records, so recall depends on distance | add records from well-assembled genomes |
| **C. Profile or class gap** | calls are made, but the idiomorph class or the profile cannot be assigned | build a lineage-specific profile |
| **D. Expected absence** | the clade is believed not to carry a canonical MAT locus | record the absence; check existing calls for non-MAT hits |
| **E. Check before curating** | the evidence of a locus is weak or unverified | verify first |

## Ascomycota: Saccharomycotina

| clade | kind | evidence | what exists | what is missing | source |
|---|---|---|---|---|---|
| **Lipomycetales** (*Lipomyces*) | A | 59 genomes in the retest set; 52 lost a `MAT`/`MATtub` call to the subphylum guard (PR #106). The 96 lost calls hold a core gene plus SLA2 and APN2. | none | records and a family scoped to Lipomycetales. Whether those loci are real MAT genes is **unverified** (SLA2-anchored hits, median 33% identity, Pfam not run). | `2026-10-17_dipodascales-curation.md` |
| **Blastobotrys, Sugiyamaella, Wickerhamiella, Starmerella and other Dipodascales** | A | 20 Dipodascales genomes lose a call within 30 kb of SLA2; the retest set holds 36 *Blastobotrys*, 24 *Sugiyamaella*, 58 *Wickerhamiella*, 35 *Starmerella* genomes. *Wickerhamiella* and *Starmerella* have no reference. | `MATyl` (*Yarrowia*) and `MATdc` (*Geotrichum*) | records and families per genus group; a Dipodascales-specific alpha-box and HMG-box profile so `MATdc` can get an `idiomorph_class`. | same note; tracker F4 |
| **Pichiales other than *Komagataella*** | A | 575 Pichiales genomes were routed by `phylum_fallback` in v0.6.0 | `MATkp` (*Komagataella* only; branch `komagataella-curation`, not yet merged) | a Pichiales/Pichiaceae profile. The MATa2 profile built from seven Saccharomycotina proteins separates only some lineages (margins +22.6 *Geotrichum* to −1.8 *Clavispora* bits), so it is supporting evidence only. | `2026-10-18_komagataella-mata-analysis.md` on `komagataella-curation` |
| **Ogataea** | B, C | 112 genomes; `MATog` records (deposit AJ617305.1 and DL-1) exist on `komagataella-curation`, **scope deliberately held back**. MATa2 is unresolved (divergent or split ORF, 14.0 and 18.1 bits against the profile). | two records | leave-own-record-out (running) to say whether the 45–75% identity hits deserve `high` confidence; records from more species. | same note |
| **Phaffomycetales, Ascoideales, Saccharomycodales, Trigonopsidales, Sporopachydermiales, Alloascoideales** | A | 187 Phaffomycetales genomes and most of the rest of the 1,472 `phylum_fallback` genomes in v0.6.0 | none | checking what the fallback calls (if any) are before curating; for the group as a whole, the fallback is searching every Ascomycota family on them. | dipodascales note, section on `phylum_fallback` |
| **Serinales `MTL`** | B | held out at order radius, 1 of 6 found with one reference left | `MTL` | more references across Serinales genera. | `analysis/2026-10-17_holdout-all-phyla.md` |
| **`MATsc` silent cassettes** | C | *S. cerevisiae* HMRa and *K. lactis* MATa score as wrong because the runner scores the record's own idiomorph but the detector calls the active locus | `MATsc` | a scoring decision, not a record: score each cassette by its own role. | same analysis |

## Ascomycota: Pezizomycotina

| clade | kind | evidence | what exists | what is missing | source |
|---|---|---|---|---|---|
| **Xylariales** | D, E | 0 of 184 uncalled genomes carry a MAT-like domain next to the flanks, as Robinson and Natvig 2019 report. But **11 of 12 called Xylariales controls also show none**, so existing calls may rest on a non-MAT HMG protein. | none | inspect the called genomes; decide whether to withhold or mark them. | `analysis/2026-10-18_uncalled-triage.md` on `uncalled-triage`; tracker F6 |
| **Microascales** (beyond *C. fimbriata*) | B | 47 uncalled; 19 have a linked MAT-like gene, 17 strong. Clearest recall target in the triage. | *C. fimbriata* unswitched record (PR #98); a switched candidate is `needs_review` | records from further Microascales genera. | triage; tracker F5 |
| **Pleosporales other than *Parastagonospora*** | B, E | 150 uncalled; 11 linked, 9 strong. 49 flank-carried and 32 zero-modelled uncalled in the earlier trace. | *Parastagonospora* fixed by records; *P. nodorum* flank records are on `pnodorum-flanks`, **not merged** | the trace of the other uncalled Pleosporales you deferred until the single-gene rule was evaluated (the rule stays off for now). | triage; tracker F1b |
| **Hypocreales, Chaetothyriales** | E | the flank search, not the MAT search, is the limit: Hypocreales 32 and Chaetothyriales 6 uncalled genomes have no flank found | | better flank references for the failing genera. | triage |
| **Mycosphaerellales, Magnaporthales, Pezizales** | E | Mycosphaerellales 16 linked (1 strong); Magnaporthales 3; Pezizales 1 | | the 37 strong candidates across orders are not individually inspected. | triage |
| **Pezizales (`MATtub`)** | C | the Pezizomycotina MAT profiles do not fit `MATtub` genes, so the triage controls fail 0 of 12 | `MATtub` records | a Pezizales-specific profile. | triage |
| **Sclerotiniaceae fused or homothallic loci** (*Sclerotinia*, *Monilinia*, *Polytolypa*, *Calogaya*) | E | the flip review's four wrong flips are fused loci; `full_length_models` blocks `homothallic_candidate` in 47 of 49 loci | *Botrytis* record (SAS405) | reference structures for the fused arrangement; the cross-match exclusion question. | tracker, revisit list |
| **Orbiliomycetes** | A | 67 Orbiliales genomes in the `phylum_fallback` set; two `MATsc` fragment calls were removed by the guard without inspection | none | references; GCA_047671915.1 and GCA_047716685.1 if someone wants to look. | tracker, revisit list |

## Ascomycota: Taphrinomycotina

| clade | kind | evidence | what is missing |
|---|---|---|---|
| ***Taphrina*** | B | *T. deformans* `PM` is missed at every holdout radius; the only references left are *Pneumocystis* | a *Taphrina* or Taphrinales reference with a locus map. (`analysis/2026-10-17_holdout-all-phyla.md`) |

## Basidiomycota

| clade | kind | evidence | what exists | what is missing | source |
|---|---|---|---|---|---|
| **Sporidiobolales (*Rhodotorula* and relatives)** | E | `redHD` (5 records) and `redPR` (13) are scope-only; `redPR` is 5 of 5 at genus radius with no reference at wider radii; **`redHD` has no held-out result**. Which orders the *Rhodotorula* genomes fall in is being checked. | records | held-out results; genomes outside the scope that are never searched (agent running). | `db/Basidiomycota/order.yml` |
| **Pucciniales (`rustHD`)** | E | 4 records, held-out 4 of 4 at genus and family radius, no reference at order radius | records | a test across the 126 Pucciniales genomes (agent running). | holdout analysis |
| **Trichosporonales and Cantharellales** | A | 126 and 97 genomes with no curated record, per the 2026-09-20 note; `db/Basidiomycota` still has no directory for either. | none | records. | `2026-09-20_basidiomycota-order-testing.md` |
| **`aLocus` at family radius** | B | found at genus radius but never when the family is withheld (0 of 8) | records | more families. | holdout analysis |
| **Lineages that the receptor references mislabel** (*Rhizoctonia*, *Auricularia*, *Fomitiporia*, *Dacryopinax*, *Exobasidium*, *Microbotryum*, *Leucosporidium*) | A | STE3 receptor homology gave false `redPR` and `wallMAT` calls there, which is why `redHD`/`redPR` are scope-only | | these lineages' own loci. | `order.yml` comment on `redHD` |

## Mucoromycota

Not a family gap: the family covers Mucorales and Umbelopsidales. The open items are reference gaps inside it
(`results/2026-09-29_two_idiomorphs/`, `analysis/2026-10-17_mucoromycota-min-margin-floor.md` on the floor-study branch):

| item | evidence |
|---|---|
| ***Zygorhynchus heterogamus* sexP** | the published sexM–sexP pair, 5.3 kb apart, is called Minus only; sexP is a detection miss |
| **Probable unrecovered sexM** | *R. nodosus*, *A. zychae* and *A. spinosa* var. *biappendiculata* sit at the top of the paralog control; the *Syzygites megalocarpus* Minus locus is withheld at the fraction floor |
| **11 fragment-typed weak Minus calls** | *M. genevensis* and *M. philippovi* have no model to score; you decided no second paralog class |

## Suggested order

My judgement, to be changed freely:

1. **Lipomycetales and the Dipodascales genera.** They are the largest known loss (72 genomes lose calls), the cause is
   understood (PR #106 guard), and the first step is a check, not curation: do the SLA2-anchored hits carry an alpha-box
   or HMG-box?
2. **Pichiales (with *Ogataea* once the leave-own-record-out result is in).** 575 genomes; `MATkp` is the start.
3. **Microascales.** The clearest recall target of the triage, and a family that has one real record to build on.
4. **Xylariales: inspect the called genomes.** It is a check of existing calls, and cheap.
5. **Basidiomycota: Trichosporonales and Cantharellales**, after the scope-only family check reports.
6. **Pezizales `MATtub` profile and the Sclerotiniaceae fused loci**, which need method work more than records.

## What this list does not cover

- **Genomes outside the sampled orders.** The triage covered ten lineage-routed orders (754 uncalled genomes); the others were not triaged.
- **Staleness.** "Uncalled" is the v0.6.0 status. Several fixes since (polish-cap, remnant step 2, the `Botrytis`/partial-core rule, the curated records) change which genomes are uncalled. Re-derive this list from the clean run.
- **Verification of the agent-reported items.** Most gene counts for Botrytis, the Pleosporales trace and the Dipodascales investigation were produced by subagents and were checked only where a note says so.
