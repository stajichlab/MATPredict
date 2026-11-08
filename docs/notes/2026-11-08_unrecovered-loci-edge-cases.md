# MAT loci not recovered, and what could improve them (v0.8.1 run census)

Date: 2026-11-08. Data: the v0.8.1 full run (`results/2026-11-05_v081_*`), 23,579 reports; scripts and tables in
`results/2026-11-07_edge_cases/` on the HPCC (`edge.py`, `edge2.py`, `edge.out`, `edge2.out`, `nocall_by_genus.tsv`, `fdb_nocall.py`).
This is an analysis of what the pipeline reports, not a ground-truth benchmark: a genome with no call may truly lack a recoverable MAT locus
(asexual lineages, lost or deleted loci, an assembly gap). The only independent checks used here are the FungiDB annotation of 313 genomes and
biology I name explicitly as "to verify".

## 1. How many are not recovered

| Phylum | Genomes searched | No call | Share |
|---|---|---|---|
| Ascomycota | 19,285 (58 not searched: no taxid or phylum) | 1,667 | 8.6% |
| Basidiomycota | 3,263 | 268 | 8.2% |
| Mucoromycotina | 973 | 105 | 10.8% |
| All | 23,521 | 2,040 | 8.7% |

Calls that exist but are weak: of 20,985 Ascomycota loci, 4,054 are `idiomorph_gene_only` and 1,425 `partial_locus`; the idiomorph is
`undetermined` in 473 Sordariomycetes and 226 Saccharomycetes loci. In Basidiomycota 4,143 Agaricomycetes loci are `undetermined` (expected
for a multiallelic A/B system) and 3,050 loci are `idiomorph_gene_only`.

## 2. The unrecovered genomes almost all have a hit that was withheld

Of the 2,040 no-call genomes, only 27 have no hit at all. The rest carry at least one withheld locus (`suppressed_loci`). The shape of the best
withheld locus in each genome:

| Best withheld locus | Genomes | Share | What it means |
|---|---|---|---|
| 2+ genes found, exactly 1 polished, `modelled_gene_bar` | 1,258 | 62% | one real gene model plus other hits no tool could model |
| 2+ genes found, 0 polished, `modelled_gene_bar` | 534 | 26% | both core genes found by search, neither modelled |
| `flank_carried_core_outside_flank_span` | 203 | 10% | flanks modelled, the core gene lies outside the flank span |
| no hit at all | 27 | 1% | genuinely nothing found |
| other (one gene, mat_gene_gate, fraction floor) | about 18 | under 1% | |

The reporting bar is `MIN_POLISHED_GENES = 2` (`detect/pipeline.py`, curator ruling 2026-09-22: loci with exactly one polished gene gave no
high-confidence call). So 88% of the unrecovered genomes (the first two rows) are cases where the search found the genes but the polishing
tools (exonerate, miniprot) built fewer than two gene models. The limiting step is gene **modelling** of divergent MAT genes, not search.

## 3. The failure is concentrated in a small number of groups

43 genera with at least 8 genomes have at least half of their genomes unrecovered (full list in `edge2.out`):

- **Dipodascales yeasts:** *Wickerhamiella* 58/58, *Starmerella* 35/35, *Sugiyamaella* 24/24, *Magnusiomyces* 15/15, *Entelexis* 12/12,
  *Zygoascus* 8/8, *Blastobotrys* 30/36, *Geotrichum* 21/91. **Lipomycetes:** *Lipomyces* 37/37, *Myxozyma* 11/14.
- **Orbiliomycetes (nematode-trapping fungi):** *Orbilia* 28/28, *Arthrobotrys* 16/17, *Dactylellina* 9/10.
- **Xylariales:** *Eutypa* 40/41, *Hypoxylon* 31/33, *Daldinia* 18/18, *Annulohypoxylon* 16/16, *Xylaria* 24/56. In *Hypoxylon* 12 of 31 are the
  `flank_carried_core_outside_flank_span` shape.
- **Chaetothyriales:** *Phialophora* 21/21, *Fonsecaea* 7/8, *Capronia* 7/9, *Exophiala* 19/114.
- **Dothideomycetes:** *Neofusicoccum* 37/42, *Ophiosphaerella* 11/11, *Sphaerulina* 20/28. **Leotiomycetes:** *Podosphaera* 11/13 (powdery mildew).
- **Saccharomycetes:** *Hanseniaspora* 53/61, *Nakaseomyces* 31/110, *Saccharomyces* 161/1,508 (the largest single group by count), *Torulaspora*
  18/79, *Saccharomycopsis* 28/72. **Serinales:** *Yamadazyma* 29/67 (25 are the `flank_carried` shape; compare #260).
- **Taphrinomycetes:** *Taphrina* 12/13. **Pyricularia** 27/520.
- **Basidiomycota:** *Microbotryum* 29/38, *Cutaneotrichosporon* 23/37, *Trichosporon* 19/41, *Apiotrichum* 16/28, *Naganishia* 14/28,
  *Filobasidium* 6/8, *Cystobasidium* 11/16.
- **Mucoromycotina:** 105/973 (86 are the one-polished shape; 7 stopped by `mat_gene_gate`).

## 4. FungiDB cross-check (313 genomes, independent annotation)

31 searched genomes had no call (8 more were not searched because GmG gave no phylum). Their best withheld loci fall into the same groups:

1. **Both core genes found, none modelled:** the Chaetothyriales (*Cladophialophora*, *Exophiala*, *Fonsecaea*, *Hortaea*) and *Bipolaris*
   (maydis, sorokiniana) and *Curvularia clavata* carry `MAT1-1-1`, `MAT1-2-1` and often `MAT1-1-2/3/5` in withheld loci with 0 polished genes.
   *Bipolaris* species are well-studied heterothallic fungi with curated references in the database, so a call is expected: this is the first
   case to take apart (to verify: whether the genes are present but modelled below the bar, or the heterothallic pair is a single idiomorph
   plus a paralog hit).
2. **Flanks modelled, core not:** *Podosphaera aphanis* (APN2, SLA2, COX13 all polished; MAT1-2-1 and MAT1-1-2 hit but unmodelled) and
   *Lodderomyces elongisporus* (PAP1, OBP1, PIK1 polished, MTLA2 found, `flank_carried_core_outside_flank_span`).
3. **Family pieces split across loci:** *S. pombe*, *S. octosporus*, *S. japonicus* (`matMc` and `matPi` in separate withheld loci, one polished
   each), *Pneumocystis oryctolagi* (`matMc` and `matPi` apart), *Taphrina deformans*, and *Puccinia striiformis* (`bE` and `bW` apart).
4. ***Lichtheimia* (Mucoromycota):** `sexM`, `sexP`, `tptA` found, 3 genes polished, still withheld (7 modelled_gene_bar, 1 mat_gene_gate,
   16 below_fraction_floor in L. corymbifera): a case where the gate, not the modelling, decides.
5. **Not searched:** *Coprinopsis cinerea* (2), *Malassezia globosa*, *Naganishia albida*, *Candida haemulonis*, *Rhizophagus irregularis* (3):
   GmG made no phylum call. In the BFD run a taxid routes these, so this is a FungiDB-run effect, not a detector gap.

## 5. What could improve recovery, ranked

Effort is the work; risk is the chance of a false MAT call. Each action keeps the repository's rule: opt-in, measured, genome-by-genome review.

| # | Action | Targets | Effort | Risk | Notes |
|---|---|---|---|---|---|
| 1 | Take apart *Bipolaris maydis* C4 and *B. sorokiniana* (why a curated, heterothallic genus gives no call) | Dothideomycetes, 137 no-call | small | low | decides whether the cause is polishing, the roster, or the both-core-genes cross-hit |
| 2 | Polishing for divergent MAT genes: a miniprot-first or protein-profile model path for genes with search identity under about 50% | the 534 zero-polished and many of the 1,258 | medium | low | the limiting step is modelling; improving it changes no rule, but must pass the identity test |
| 3 | `--strong-single-gene-policy` gate check (agreed: modelled identity >= 85 and coverage >= 90) on the one-polished genomes | the 1,258 | small (check running, #261) | medium | the rule exists and is off; this measures how many it would recover |
| 4 | Curate references for the unreferenced groups: Dipodascales and Lipomycetes yeasts, Orbiliomycetes, Xylariales, Chaetothyriales, Botryosphaeriales, Podosphaera | about 1,000 genomes | large (curation) | low | each record raises search sensitivity where none exists; follows the #189/#202/#260 pattern |
| 5 | Review the one-polished-gene ruling for a flank-supported tier: one polished MAT gene plus a polished flank gene reports `low`, not nothing | part of the 1,258 | small | medium | needs a curator decision; the 2026-09-22 ruling was made when no tier below the bar existed |
| 6 | `flank_carried_core_outside_flank_span`: widen the span or allow the core gene beyond it when the flanks are polished | the 203 (Hypoxylon, Yamadazyma, Pyricularia) | medium | medium | overlaps #260 |
| 7 | Family-structure fixes: join `matMc`/`matPi` pieces of the PM family and `bE`/`bW` of rustHD across a split | S. pombe group, Pneumocystis, Taphrina, rust | medium | medium | relates to the silent-pair work (#252 revert) |
| 8 | A per-genome "why no call" field naming the best withheld locus and the reason | all | small | none | report-only; turns this census into a column |

Suggested order: 1, 3, 8 (cheap, informative), then 2 and 4 (largest effect), then the curator-gated 5 to 7.

## 6. What this census cannot say

- Whether a no-call genome truly lacks a MAT locus. The Orbiliomycetes, Dipodascales and Lipomycetes groups at 100% no-call may be missing
  references (a detector gap) or may use a structure the families do not describe (a curation gap); only curated records or the literature settle it.
- The share that a relaxed bar would recover and at what precision. That needs an opt-in run with measured precision, as for every other rule.
- Anything about genomes with a call but the wrong idiomorph; that is a separate question (#138, #145, #255, #275).

## 7. For the curator

1. Rule on action 5 (a tier for one polished gene plus a polished flank, instead of withholding).
2. Choose which groups to curate first (action 4); the Dipodascales yeasts and Xylariales have the most genomes.
3. Say whether the per-genome "why no call" field (action 8) belongs in v0.9.0.
