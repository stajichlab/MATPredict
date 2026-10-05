# MATPredict campaigns: overview, numbers and where to look

Status: summary of existing runs (2026-10-05); no new compute. Interactive version: the campaign dashboard
(`results/2026-10-05_campaign_overview/dashboard.html`, built by `make_dashboard.py`; published as a private artifact).
Figures: `results/2026-10-05_campaign_overview/fig1_*.svg` to `fig3_*.svg`; numbers: `campaign_summary.tsv`.

## Campaigns on v0.6.0 (BFD genomes)
| Campaign | Genomes | With a call | Loci | High / medium / low | Directory in `results/` |
|---|---|---|---|---|---|
| Ascomycota | 19,415 | 17,314 (89.2%) | 21,011 | 13,740 / 7,112 / 159 | `2026-10-03_ascomycota_v060` |
| Basidiomycota | 3,270 | 2,982 (91.2%) | 5,900 | 3,813 / 2,077 / 10 | `2026-10-03_basidiomycota_v060` |
| Mucoromycotina (BFD rows) | 289 | 256 (88.6%) | 265 | 204 / 60 / 1 | `2026-10-03_mucoromycotina_mat`, re-runs `2026-10-04_mucoro_bfd_*` |
| Alaninales re-run (20 code-26 genomes) | 20 | 15 (exonerate string) | n/a | n/a | `2026-10-04_alaninales_v061*` |
Ascomycota also had 2,066 genomes without a call, 20 failures (the code-26 bug, since fixed) and 15 with no report. The Mucoromycotina table also holds
LCG (656) and Jena (66) held-out rows; the BFD rows are used here. Older code versions: `2026-09-26_basidiomycota_full`, `2026-09-26_early_diverging`
(Mucoromycota 227 of 293, Mortierellomycota 6 of 100, Kickxellomycota 7 of 190), `2026-09-21_tremellales334`, `2026-09-21_zygo23`. Existing write-up of the
v0.6.0 runs: `analysis/2026-10-04_basidiomycota-ascomycota-v060-campaign.md`.

## What "called" and "recovered" mean
- A genome is called if it has at least one MAT call. That is detection, not truth. Most clades have no curated reference, so an uncalled genome is a reference
  gap until shown otherwise. Orders without their own record are searched through the phylum fallback.
- Locus class (per locus): full locus (`mat_locus`) 14,657 Ascomycota, 2,721 Basidiomycota, 218 Mucoromycotina; idiomorph gene only 4,466 / 2,955 / 1;
  partial locus 1,531 / 224 / 46; homothallic candidate 357 Ascomycota. Half of the Basidiomycota loci are the idiomorph gene without a recognised locus.
- Basidiomycota calls include PR-only calls, 1,360 of them from the strict-CAAX scan alone (`verification: unverified`). Without them the rate falls sharply in
  the phylum-fallback orders (Cantharellales 91% to 17%, Trichosporonales 47% to 8%, Cystofilobasidiales 50% to 17%, Filobasidiales 46% to 20%).

## Call rates by group
- Ascomycota classes with 20 or more genomes (14 classes, 99% of genomes): Schizosaccharomycetes 98.5%, Eurotiomycetes 96.5%, Sordariomycetes 93.0%,
  Dothideomycetes 89.4%, Saccharomycetes 85.1%, Pezizomycetes 79.7%; fallback-routed: Dipodascomycetes 28.7%, Orbiliomycetes 9.0%.
- Basidiomycota orders with 40 or more genomes (16 orders): lineage-routed orders 83% to 100% (Ustilaginales and Wallemiales 100%, Agaricales 98.8%,
  Pucciniales 82.9%).
- Group counts for choosing figure sets (genomes per level with 20 or more): Ascomycota 14 classes, 41 orders, 70 families (of 284); Basidiomycota 9 classes,
  18 orders, 34 families (of 177). The top 20 families cover 76% (Ascomycota) and 67% (Basidiomycota) of genomes, so family level is for drill-downs, not for
  one panel per group.

## Locus size and content
- Ascomycota full generic-MAT loci, median length: Sordariomycetes 17.1 kb, Eurotiomycetes 13.5, Leotiomycetes 13.1, Lecanoromycetes 12.7, Dothideomycetes 8.9.
  Share carrying both APN2 and SLA2: 82 to 91% in all of these except Dothideomycetes (6%).
- Mucoromycotina (BFD): median 14.5 kb overall; Umbelopsis 41.5 kb, Rhizopus 14.5, Mucor 14.8, Backusella 9.2, Cunninghamella 10.4, Apophysomyces 7.8,
  Syncephalastrum 5.3. btbA appears only in Rhizopus (77% of its loci). Per-genus gene presence is a table in the dashboard.

## Xylariales synteny (from PRs #32 and #34)
- 257 genomes, 164 species; 41 genomes are *Eutypa lata* strains. By species: outgroup-like order (SAC) 69, COX13 beside SLA2 (SCA) 73.
- Two derived junctions relative to outgroups (CIA30 beside SLA2; COX13 beside SLA2 in SCA); COX13 and APN2 always stay together. See
  `analysis/2026-10-05_xylariales-synteny.md`.

## Exemplars and outliers (more in the dashboard)
Exemplars: Wallemiales 0 to 51 of 51 from one record; Rhodotorula P/R 61 of 62 held-out strains; Zygo 23 of 23. Outliers: Orbiliomycetes 9% and
Dipodascomycetes 29% (fallback); PR-only calls in fallback orders; Dothideomycetes flank genes split; early-diverging fungi with few calls; Xylaria NC1011 with
intact flanks and no MAT gene; Mycotypha sexM 153 kb from sexP.

## Where the full reports are
- Per campaign, `reports_all.tar.zst` holds one `detection_report.yaml` per genome (and `wall_seconds`); `genomes.tsv`, `loci.tsv` and the by-class or
  by-order tables are the fast entry. Per-genome GFF3 and FASTA are not committed; single runs write `detected_loci.gff3` and `detected_loci.fasta`.
- Drawing: `matpredict curate-db draw-locus` and `draw-synteny` for curated records (`docs/HANDOFF-visualization-followups.md`). Campaign figures are scripts in
  `results/`. The aggregate run report and per-locus report in the packaging spec (`docs/superpowers/specs/2026-10-04-packaging-service-and-reports-design.md`,
  workstreams C and D) are planned, not built.

## Should the Basidiomycota campaign be re-run?
Not yet, on the evidence in the repo.
- Nothing that changes Basidiomycota calls has merged since the v0.6.0 run: the code changes since 2026-10-03 are the code-26 fix (Ascomycota), frameshift-aware
  classification and the Mycotypha record and rebuild (Mucoromycota), the *F. oxysporum* records (Ascomycota), an *Umbelopsis* override and the offline
  taxonomy. The 33-genome Basidiomycota panel in the *F. oxysporum* regression shows 0 changed loci.
- The open Basidiomycota questions are decisions and curation, not code: whether any order may lose the `unverified` label (the draft PR #33 reports 6 of 9
  independent mating receptors flagged against 1 of 46 non-mating copies, with 179 of 1,359 unverified calls at chance level, and changes no label), and new
  records for the fallback orders (Trichosporonales, Cantharellales, Cystofilobasidiales, Filobasidiales, Sebacinales). `analysis/decisions.md` already records
  "no full Basidiomycota re-run".
- A re-run is cheap (the whole phylum took three jobs of about 1 hour and one of 32 minutes on `exfab`), so cost is not the barrier. Re-run the affected orders after new records or a
  label ruling, gated by the regression panel, and the whole phylum only if detection code changes.

## Limits
Counts are descriptive. Several runs used earlier code (see dates). Genome counts reflect uneven strain sampling (for example 41 *Eutypa lata* genomes).
The dashboard was checked at desktop and phone width; it shows committed results only.
