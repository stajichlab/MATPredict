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
tools (exonerate, miniprot) built fewer than two gene models. Whether the limiting step is gene **modelling** or the **sensitivity of the search** is not settled by these counts: the unmodelled hits are mostly weak fragments (section 5a measured this for the best cases: 23-47% identity, 12-53% coverage, 28-41 bits), so a divergent gene that the search barely detects looks the same as a gene the polishing tools cannot build.

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

### 4a. Two causes tested on those 31 genomes (frozen v0.8.1 code, same database)

| Test | Result |
|---|---|
| Drop the annotated proteome (`--proteins`), genome-only | 2 of 31 recover (*Cladophialophora carrionii* KSF, medium; *Hortaea werneckii* EXF-2000, high); 29 stay uncalled |
| Also remove the polish cap (`--max-polished-clusters-per-family 0`) | identical: the same 2 recover, 29 stay uncalled |

So the annotated-proteome path costs a small number of calls (2 of 31), and the polishing budget is **not** why these genomes are missed (a
hypothesis from *B. maydis* C4's 35 withheld clusters, tested and rejected). *Bipolaris* and *Curvularia* are called in 82 of 82 BFD genomes in
v0.7.0, v0.8.0 and v0.8.1, so the *B. maydis* C4 failure is specific to that assembly: its MAT1-1-1 + SLA2 and MAT1-2-1 + SLA2 clusters (on MU825198
and MU825193, the second with one polished gene) are the candidate locus, and neither reaches two modelled genes. The cause is not determined here.
The 29 stay in the "modelling" group of section 2.

## 5. What could improve recovery, ranked

Effort is the work; risk is the chance of a false MAT call. Each action keeps the repository's rule: opt-in, measured, genome-by-genome review.

| # | Action | Targets | Effort | Risk | Notes |
|---|---|---|---|---|---|
| 1 | Take apart *Bipolaris maydis* C4 and *B. sorokiniana* (why a curated, heterothallic genus gives no call in these two FungiDB assemblies when 82 of 82 BFD *Bipolaris* and *Curvularia* genomes are called) | Dothideomycetes, 137 no-call | small | low | the proteome path and the polish cap are ruled out (section 4a); next: the genomic window of the two SLA2-flanked clusters, hand-modelled against a curated record |
| 2 | Polishing for divergent MAT genes: a miniprot-first or protein-profile model path for genes with search identity under about 50%. **Check first (section 5a): in the best cases the unmodelled partner is only a weak search fragment, so search sensitivity may matter more than polishing** | the 534 zero-polished and many of the 1,258 | medium | low | the limiting step is modelling; improving it changes no rule, but must pass the identity test |
| 3 | `--strong-single-gene-policy` (gate: modelled identity >= 85 and coverage >= 90), **measured for #261: not the lever** | the 1,258 | done | medium | on 3,100 candidate genomes the rule recovers 5 loci (each `low`, `idiomorph_gene_only`, no flank genes, none with a truth label) and changes no existing call; 117 lone core hits at >= 85% in 96 genomes never reach the model gate because they are 90-430 bp fragments; one of the 5 recovered models is defective (a frameshift, no stop), so any use needs a defective-model gate (`results/2026-11-07_issue261/FINDINGS.md`) |
| 4 | Curate references for the unreferenced groups: Dipodascales and Lipomycetes yeasts, Orbiliomycetes, Xylariales, Chaetothyriales, Botryosphaeriales, Podosphaera | about 1,000 genomes | large (curation) | low | each record raises search sensitivity where none exists; follows the #189/#202/#260 pattern |
| 5 | Report a locus with exactly one polished gene at `low` instead of withholding it. **Measured, 2026-11-08: not worth a general rule** (section 5a) | the 1,264 | done | high | 92% of the polished core genes are under 60% identity; a narrow gate recovers 12 genomes |
| 6 | `flank_carried_core_outside_flank_span`: widen the span or allow the core gene beyond it when the flanks are polished | the 203 (Hypoxylon, Yamadazyma, Pyricularia) | medium | medium | overlaps #260 |
| 7 | Family-structure fixes: join `matMc`/`matPi` pieces of the PM family and `bE`/`bW` of rustHD across a split | S. pombe group, Pneumocystis, Taphrina, rust | medium | medium | relates to the silent-pair work (#252 revert) |
| 8 | A per-genome "why no call" field naming the best withheld locus and the reason | all | small | none | report-only; turns this census into a column |

Suggested order: 1, 8 (cheap, informative), then 2 and 4 (largest effect), then the curator-gated 5 to 7. Action 3 is measured and small (5 loci), so the rule does not address the one-polished-gene genomes: those have other, unpolished hits and no lone modelled gene to promote.

### 5a. The one-polished-gene measurement (action 5)

The bar counts any polished gene, so a locus with a polished core gene and a polished flank gene already passes; what is withheld is a locus with
**exactly one** polished gene. With `--withheld-detail-policy on` (report-only) the 1,935 no-call Ascomycota and Basidiomycota genomes of the v0.8.1 run
were rerun with the same code (`results/2026-11-08_withheld_detail/`; calls unchanged by construction). Of 162,417 withheld `modelled_gene_bar` loci,
160,051 (98.5%) have no polished gene and 2,366 have exactly one.

| | Loci | Genomes (best withheld locus) |
|---|---|---|
| the one polished gene is a flank gene (SLA2, APN2, COX13, PAN6, FAO1 ...) | 1,202 | 557 |
| the one polished gene is a core gene | 1,164 | 707 |

For the 707 genomes with a polished core gene: identity under 60% in 653 (92%), 60-85% in 23, 85% or more in 31. In 545 the locus also holds another
core hit (a possible opposite idiomorph) and in 495 it holds no flank hit at all. With identity at least 60% and coverage at least 60%, 35 genomes
remain, and 23 of those rest on `RPL39`, a conserved ribosomal gene that is not mating-type specific. The other 12 are *Taphrina* `matMc`
(100% identity, `matPi` unpolished; 5), *Sphaerulina musiva* `MAT1-1-1` (3), *Flavoplaca navasiana* `MAT1-1-1` (1), *Asteromyces cruciatus* and
*Sporothrix epigloea* `MAT1-2-1` (2) and *Melampsora aecidioides* `bW` (1).

So a general `low` tier would add about 700 loci that are mostly weak models with no flank support (a false-call risk the 2026-09-22 ruling was made to
avoid), and a gate that keeps only well-modelled specific genes recovers 12 genomes, about 0.6% of the unrecovered ones. The "unpolished partner" in these 12 cases is not a gene the polishing tools failed to build: it is a weak search fragment (identity 23-47%,
coverage 12-53%, 28-41 bits; for example *Taphrina deformans* `matPi` at 33%, 18% coverage and 29 bits beside `matMc` at 100%, *Melampsora aecidioides* `bE`
at 37.5%, 16% and 32 bits beside `bW` at 75%). So where one gene is modelled and its partner is not, the partner is barely detected at all, and the
lever is the **sensitivity of the search** for the divergent partner (a profile or HMM search, more references, or a lower threshold in a window already
known to hold a modelled MAT gene), not the polishing tools. That belongs with action 4 (references) and a search-sensitivity experiment, not action 2.

## 6. What this census cannot say

- Whether a no-call genome truly lacks a MAT locus. The Orbiliomycetes, Dipodascales and Lipomycetes groups at 100% no-call may be missing
  references (a detector gap) or may use a structure the families do not describe (a curation gap); only curated records or the literature settle it.
- The share that a relaxed bar would recover and at what precision. That needs an opt-in run with measured precision, as for every other rule.
- Anything about genomes with a call but the wrong idiomorph; that is a separate question (#138, #145, #255, #275).

## 7. For the curator

1. Rule on action 5 (a tier for one polished gene plus a polished flank, instead of withholding).
2. Choose which groups to curate first (action 4); the Dipodascales yeasts and Xylariales have the most genomes.
3. Say whether the per-genome "why no call" field (action 8) belongs in v0.9.0.
