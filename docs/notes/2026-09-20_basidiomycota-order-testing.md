# Basidiomycota order-by-order detection testing (2026-09-20)

Branch `basidiomycota-order-testing`, based on `detect-targeting` (34117e9, the
PR #2 merge). Baseline before any change: **463 tests passing**.

All detection work before this note was validated on Mucoromycota only. This
note is the first Basidiomycota measurement. Read the audit section first: it
bounds what any recall number can possibly mean.

**Sample size warning, stated once and applying throughout.** Basidiomycota has
**10 accepted curated records across 5 species in 3 orders** (Agaricales 6,
Tremellales 2, Ustilaginales 2). Every per-locus number below is a count out of
1 or 2. Nothing here supports a claim about Basidiomycota as a whole.

---

## 1. Audit: what can be detected in principle

### 1.1 Roster genes with no reference protein

Measured with `detect.reference_fasta.searchable_genes_by_family` against the
live `db/`, narrowed to the 9 Basidiomycota families (40 reference proteins).

| locus | curated records | roster genes | searchable | **unsearchable** |
|---|---|---|---|---|
| Aalpha | 1 | 2 | 2 | — |
| **Abeta** | **0** | 2 | **0** | **`Y`, `Z`** |
| Balpha | 1 | 4 | 3 | `pheromone_receptor` |
| Bbeta | 1 | 9 | 8 | `pheromone_receptor` |
| HD | 2 | 2 | 2 | — |
| MAT | 2 | 9 | 9 | — |
| PR | 1 | 4 | 4 | — |
| aLocus | 1 | 3 | 3 | — |
| bLocus | 1 | 2 | 2 | — |
| **total** | **10** | **37** | **33** | **4** |

**These 4 are NOT the `algA`/`glrA` defect.** That defect was an unsearchable
gene inflating the scoring denominator of a family that *did* produce hits.
Both mitigations are already in place here:

* `scoring.score_cluster` skips any family with zero hits (`if not found:
  continue`), so **Abeta**, having no reference protein at all, can never be
  scored and never inflates anything. It appears in the report only as a
  `not_detected` entry: *"no reference-protein hits found for this family in
  this genome"*.
* For Balpha/Bbeta, `searchable_genes_by_family` puts `pheromone_receptor` in
  `genes_not_searchable`, which `score_cluster` drops from the
  `fraction_found` denominator. Confirmed in the live report: the Balpha entry
  at `JAAGWA010000010.1:1806650-1826460` reports
  `genes_found=['bar3','bap3-3']`, `genes_missing=[]`,
  `genes_not_searchable=['bap3-1','pheromone_receptor']`.

**But the fix was applied to `fraction_found` only, not to
`tiering.assign_tier`** — see §4.1. That is a real bug and it is the same
defect class, one function downstream.

### 1.2 Abeta is undetectable, and worse, it is *mis*detectable

Abeta has zero curated records and zero reference proteins, so no Abeta call
can ever be made. That alone is a gap.

The sharper problem: Abeta declares gene names `Y` and `Z` — **the same two
names Aalpha declares**, for the same species (*S. commune*, taxid 5334).
`search._attribute` keys a hit on the curated record id, so a real Abeta `Y`
gene in some genome will be found via the *Aalpha* record's protein and
reported as **Aalpha**. The absence is therefore not silent, it is
mislabelling. Curating one Abeta record is the fix.

### 1.3 No Basidiomycota family declares a single flanking gene

| phylum | families | families with >=1 `flanking_conserved` gene |
|---|---|---|
| Ascomycota | 9 | 4 |
| Mucoromycota | 1 | 1 (`tptA`, `rnhA`) |
| **Basidiomycota** | **9** | **0** |

Only `aLocus.rba1` is non-core, and it is `flanking_variable`. Every other one
of the 37 Basidiomycota roster genes is `core_MAT`. Two machinery consequences,
both read straight from the code and both confirmed in the live report:

**(a) `tiering.assign_tier` loses its top discriminator.**
`has_flanking_conserved(family)` is `False` for all 9 families, so the
`elif has_flanking_conserved(family)` branch is dead and the tier is `high`
as soon as every core gene is found and polished. In Mucoromycota `high`
means "core genes *plus* a conserved flank". In Basidiomycota it means only
"every core gene of this family was found". The bar is strictly lower and the
same word is used for it.

**(b) `idiomorph.classify_locus` is inert for 8 of the 9 families.**
`classify_locus` returns `idiomorph_gene_only` only when `restricted and not
unrestricted`, where `restricted` is the set of found genes carrying
`present_in_idiomorphs`. **Only the Tremellales `MAT` family declares
`present_in_idiomorphs` at all.** For the other 8, `restricted` is always
empty, so the function falls through to `return LOCUS_CLASS_MAT`
unconditionally — a lone single-gene hit is classed `mat_locus` exactly like a
complete locus. The docstring's "`mat_locus` -- a core gene with at least one
flanking gene" is not what the code does in this phylum; it is the default
return.

Measured on *C. cinerea* (65 detected entries): `locus_class` is `mat_locus`
for 60 and `idiomorph_gene_only` for 5 — and all 5 are the Tremellales `MAT`
family misfiring on a Coprinopsis genome. `detection_pass` is `strict` for all
65. **Neither field discriminates.** `confidence` is the only field that does,
and per (a) its own top tier is weakened.

### 1.4 How much of the idiomorph work applies here

`idiomorph.assign_idiomorph` returns `undetermined` for every
`vocabulary_type: pattern` family by design. **8 of 9 Basidiomycota families
are `pattern`**; only Tremellales `MAT` is `enum` (`alpha`/`a`). So:

* `undetermined` is the correct and expected answer for HD, PR, Aalpha, Abeta,
  Balpha, Bbeta, aLocus and bLocus. It must not be scored as a failure.
* The whole sexM/sexP-style resolution stack — `resolve_idiomorph_overlaps`,
  `min_idiomorph_margin`, the `_rank` proteome-over-identity tiebreak — is
  reachable in exactly **one** of nine families. The Mucoromycota tuning
  (23/23 idiomorph correct) transfers to essentially none of this phylum.
* `classify_locus`'s `homothallic_candidate` branch is likewise reachable only
  in `MAT`, and only with a proteome (`_PROTEOME_METHODS` gate).

### 1.5 Curated ground truth, with source assemblies

| record | order | locus | taxid | source accession | contig |
|---|---|---|---|---|---|
| `5346_a43-b43-okayama-7_HD_A43` | Agaricales | HD | 5346 | **GCA_016772295.1** | JAAGWA010000001.1 |
| `5346_a43-b43-okayama-7_PR_B43` | Agaricales | PR | 5346 | **GCA_016772295.1** | JAAGWA010000010.1 |
| `5334_h4-8_Aalpha_4` | Agaricales | Aalpha | 5334 | **GCF_000143185.2** | NW_026089539.1 |
| `5334_h4-8_Balpha_3` | Agaricales | Balpha | 5334 | **GCF_000143185.2** | NW_026089548.1 |
| `5334_h4-8_Bbeta_2` | Agaricales | Bbeta | 5334 | **GCF_000143185.2** | NW_026089548.1 |
| `192523_h97_HD_A1` | Agaricales | HD | 192523 | NW_006267344.1 | NW_006267344.1 |
| `40410_jec21_MAT_alpha` | Tremellales | MAT | 40410 | AF542531.2 | AF542531.2 |
| `40410_jec20_MAT_a` | Tremellales | MAT | 40410 | AF542530.2 | AF542530.2 |
| `5270_521_aLocus_a1` | Ustilaginales | aLocus | 5270 | U37795.1 | U37795.1 |
| `5270_521_bLocus_b1` | Ustilaginales | bLocus | 5270 | NC_026478.1 | NC_026478.1 |

Only **2 of 10** records cite a whole-genome assembly accession. The other 8
cite locus-level GenBank deposits or a single RefSeq scaffold, so they can
never accession-match a BFD assembly.

All 6 relevant assemblies are present in the BFD library:

| species / strain | BFD assembly | **BFD taxid** | **curated taxid** |
|---|---|---|---|
| *C. cinerea* Okayama-7 AmutBmut | `GCA_016772295.1_ASM1677229v1` | **1132390** | 5346 |
| *S. commune* H4-8 | `GCF_000143185.2_Schco3` | **578458** | 5334 |
| *A. bisporus* H97 | `GCF_000300575.1_Agabi_varbisH97_2` | **936046** | 192523 |
| *C. deneoformans* JEC21 | `GCF_000091045.1_ASM9104v1` | **214684** | 40410 |
| *C. deneoformans* JEC20 | `GCA_056621545.1_ASM5662154v1` | 40410 | 40410 |
| *M. maydis* 521 | `GCF_000328475.2_Umaydis521_2.0` | 5270 | 5270 |

---

## 2. The self-consistency harness cannot see the reference genomes

`benchmark.match_ground_truth` gates on **exact taxid equality**
(`if doc["taxonomy"]["taxid"] != taxid: continue`), where `taxid` is parsed out
of the batch runner's `<taxid>_<accession>` directory name and therefore comes
from the BFD manifest. Four of the six reference assemblies carry a
**strain-level** taxid in BFD while the curated record carries the
**species-level** taxid. Measured, calling `match_ground_truth` directly:

| genome_id as BFD names it | curated records matched |
|---|---|
| `1132390_GCA_016772295.1` (*C. cinerea* Okayama-7) | **0** |
| `578458_GCF_000143185.2` (*S. commune* H4-8) | **0** |
| `936046_GCF_000300575.1` (*A. bisporus* H97) | **0** |
| `214684_GCF_000091045.1` (*C. deneoformans* JEC21) | **0** |
| `40410_GCA_056621545.1` (*C. deneoformans* JEC20) | 1 exact, 1 ambiguous |
| `5270_GCF_000328475.2` (*M. maydis* 521) | **2 exact** (aLocus + bLocus) |

Substituting the curated species taxid recovers them immediately, by accession
match alone:

| genome_id with the curated taxid | curated records matched |
|---|---|
| `5346_GCA_016772295.1` | **2 exact** — HD and PR |
| `5334_GCF_000143185.2` | **3 exact** — Aalpha, Balpha, Bbeta |
| `192523_GCF_000300575.1` | 1 ambiguous (record cites the scaffold `NW_006267344.1`, not the assembly) |
| `40410_GCF_000091045.1` | 2 ambiguous (records cite `AF542530/1.2`, not the assembly) |

Sweeping every BFD genome whose manifest taxid equals some curated record's
taxid, and asking `match_ground_truth` about each, gives the whole picture:

| curated record | best reachable status | via |
|---|---|---|
| `5270_521_aLocus_a1` | **exact** | `5270_GCF_000328475.2` (521) |
| `5270_521_bLocus_b1` | **exact** | `5270_GCF_000328475.2` (521) |
| `40410_jec20_MAT_a` | **exact** | `40410_GCA_056621545.1` (JEC20) |
| `40410_jec21_MAT_alpha` | ambiguous | `40410_GCA_056621545.1` (JEC20 — wrong strain) |
| `5346_..._HD_A43` | ambiguous | `5346_GCA_982397435.1` (T48-F — wrong strain) |
| `5346_..._PR_B43` | ambiguous | `5346_GCA_982397435.1` (T48-F — wrong strain) |
| `5334_h4-8_Aalpha_4` | ambiguous | `5334_GCA_019143615.1` (14-01_S62 — wrong strain) |
| `5334_h4-8_Balpha_3` | ambiguous | `5334_GCA_019143615.1` (wrong strain) |
| `5334_h4-8_Bbeta_2` | ambiguous | `5334_GCA_019143615.1` (wrong strain) |
| `192523_h97_HD_A1` | ambiguous | `192523_GCA_006491665.1` (ARP23 — wrong strain) |

**Exact-reachable: 3 of 10.** The 7 ambiguous rows are not near-misses against
the right genome — they are the *wrong strain*, reached only because it happens
to carry the species taxid, while the record's own reference assembly is
invisible.

So the failure is not "no ground truth exists". The accession evidence is
already sitting in the record and is already unambiguous — the taxid gate
throws it away before the accession is ever compared.

`route()` solved exactly this problem for family routing by matching the
NCBI lineage rather than the bare taxid, after an audit found 50 of 61 curated
records falling through to the exhaustive rule. The same reasoning applies
here and has not been carried across.

**Recommended fix (not yet applied): compare the source accession BEFORE the
taxid gate.** An assembly-accession match is already conclusive on its own —
the taxid adds nothing to it — and this is a strictly smaller change than
making the matcher lineage-aware. It would raise exact-matched Basidiomycota
records from 3/10 to 8/10 with no new false pairings. Making the taxid gate
lineage-aware as well would additionally let *A. bisporus* and *C. deneoformans
JEC21* be reached, but those would still come out `ambiguous` for a separate
reason (their records cite a sub-assembly accession), so it is the smaller half
of the win.

---

## 3. Tetrapolar detection: does a run find BOTH unlinked loci?

*Coprinopsis cinerea* Okayama-7 (`GCA_016772295.1_ASM1677229v1`), genome-only
(`tblastn`, no proteome), `--phylum Basidiomycota`, 11 min 13 s wall, 162 MB
peak RSS, 40 reference proteins.

### 3.1 Yes — both loci are found, and they are the only two high-confidence intervals

| locus | curated | detected | conf |
|---|---|---|---|
| **HD** | `JAAGWA010000001.1:1625680-1631546` | `JAAGWA010000001.1:1625683-1631394` | **high** |
| **PR** | `JAAGWA010000010.1:1806154-1826859` | `JAAGWA010000010.1:1806650-1826460` | **high** |

HD is 3 bp off at the start and 152 bp inside at the end; both genes (`HD1`,
`HD2`) found, `genes_missing=[]`. PR is 496 bp inside at the start and 399 bp
inside at the end; **all four** distinct gene names found
(`pheromone_receptor`, `pheromone_B43`, `pheromone_B44`,
`fungal_mating_type_pheromone`), `genes_missing=[]`. Both are on the correct,
different contigs. `reference_records` correctly names the matching curated
record in each case.

**Tetrapolar two-locus detection works, on this one genome, without a
proteome.**

### 3.2 But the report does not make the two-locus structure legible

The run emitted **65 `detected` entries over 45 distinct intervals**.

| field | distribution over the 65 entries |
|---|---|
| `detection_pass` | `strict` 65 |
| `locus_class` | `mat_locus` 60, `idiomorph_gene_only` 5 |
| `confidence` | `low` 54, `medium` 8, **`high` 3** |
| `idiomorph` | `undetermined` 61, `alpha` 3, `a` 1 |

The two true loci are 2 of the 3 `high` entries. Filtering on `confidence ==
"high"` yields exactly the right answer plus one duplicate. Filtering on
`locus_class` or `detection_pass` yields nothing at all.

**8 intervals are claimed by more than one family**, the same genomic region
emitted once per family as an independent top-level result:

| interval | families claiming it |
|---|---|
| `JAAGWA010000010.1:1806650-1826460` (**the true PR locus**) | PR(high), Balpha(medium), aLocus(low), MAT(low), Bbeta(low) |
| `JAAGWA010000004.1:2155782-2187065` | PR(medium), Bbeta(medium), Balpha(low), aLocus(low), MAT(low) |
| `JAAGWA010000010.1:1873681-1878716` | Balpha, aLocus, MAT, PR, Bbeta (all low) |
| `JAAGWA010000001.1:1625683-1631394` (**the true HD locus**) | HD(high), **Aalpha(high)**, bLocus(medium) |
| `JAAGWA010000004.1:3525580-3545056` | MAT(medium), Aalpha(low), Balpha(low) |
| `JAAGWA010000003.1:2195478-2195666` | HD, Aalpha, bLocus (all low) |
| `JAAGWA010000007.1:1465149-1465280` | HD, Aalpha (both low) |
| `JAAGWA010000011.1:369126-392390` | MAT, bLocus (both low) |

The mechanism is `gene_class`, not a bug in `_attribute`. HD1/HD2 is the
declared class of HD's `HD1`/`HD2`, Aalpha's `Z`/`Y`, bLocus's `bW`/`bE` and
MAT's `SXI1`/`SXI2`; pheromone_receptor/pheromone_precursor is shared by PR,
Balpha, Bbeta, aLocus and MAT. These reference proteins are genuine homologs of
each other, so one real locus lights up every family in its class group.
`_attribute` keying on the record id correctly prevents *gene-name* confusion;
it does not and cannot prevent *family* confusion.

**The report already records the collision** — the HD entry carries
`ambiguous_with: ["Basidiomycota:Aalpha", "Basidiomycota:MAT",
"Basidiomycota:bLocus"]`. What is missing is any grouping: the five rows for the
true PR locus are five peer entries in `detected`, with nothing saying they are
one locus seen five ways. A reader, or a downstream consumer, sees 65 loci in a
genome that has 2.

**This is the design gap.** The report needs a per-interval collapse — one
locus object per genomic interval, with a best-scoring family and the losing
families demoted to `also_matched` — and, for a tetrapolar genome, a
genome-level statement that two unlinked loci of complementary class (one
HD-class, one pheromone/receptor-class) were found. Neither exists. The
`ambiguous_with` field is the right raw material for the first.

### 3.2a Correct taxid routing fixes nearly all of §3.2 — measured

The same *C. cinerea* genome, same code, same database, run with
`--taxid 1132390` (BFD's own strain taxid) instead of `--phylum
Basidiomycota`. Lineage routing matches 5346 as an ancestor and narrows the
run to the HD and PR families:

| | `--phylum Basidiomycota` | `--taxid 1132390` |
|---|---|---|
| families searched | 9 | **2** (HD, PR) |
| detected entries | 65 | **14** |
| distinct intervals | 45 | 14 |
| intervals claimed by >1 family | 8 | **0** |
| `high` entries | 3 | **2** |
| spurious Cryptococcus `MAT` calls | 5 | **0** |
| wall time | 11 min 13 s | **46.8 s** |

The two `high` entries are exactly HD `JAAGWA010000001.1:1625683-1631394` and
PR `JAAGWA010000010.1:1806650-1826460` — **the correct tetrapolar pair, on
different contigs, and nothing else**. The report is legible without any new
grouping machinery. Cost falls 14x.

*S. commune* reproduces it. `--taxid 578458` routes by lineage to the four
*S. commune* sublocus families:

| | `--phylum Basidiomycota` | `--taxid 578458` |
|---|---|---|
| detected entries | 96 | **26** |
| intervals claimed by >1 family | 10 | **1** |
| wall time | 21 min 21 s | **3 min 14 s** |

All three curated loci are still recovered in full — Aalpha 2/2 `high` at
`1821010-1827453`, Balpha 3/3 and Bbeta 8/8 `medium` — and the single
remaining multi-family interval is Balpha + Bbeta on the true B locus, which
is a real sublocus co-occurrence rather than noise. The spurious HD-family
claims on the Aalpha intervals are gone, because HD is no longer searched.
`Abeta` is routed and, as predicted, contributes nothing.

**So the noise in §3.2 is not a property of the detector. It is the cost of
`--phylum`.** The per-interval collapse recommended there is still worth
building, because it is what makes a *phylum-level* run readable, but it is no
longer the first thing to fix.

**The first thing to fix is `taxonomic_scope`.** Every Basidiomycota family
scopes to a single species taxid — 5334, 5346, 5270, 40410, 192523 — whereas
Mucoromycota's scopes are broad clade taxids. Lineage routing therefore works
beautifully for the five curated species and their strains, and not at all for
anything else: of the 3,174 Basidiomycota genomes in BFD, all but a handful
fall through to `phylum_fallback` and get the 9-family, 65-entry treatment.
Widening the scopes to real clades (HD/PR to an Agaricomycete or Agaricales
taxid, `MAT` to Tremellales, `aLocus`/`bLocus` to Ustilaginales, the four
*S. commune* sublocus families to Schizophyllaceae or Agaricales) would give
a Boletales or Polyporales genome the same narrow, fast, clean run this
*C. cinerea* genome just got. That is curation data and needs the curator's
ruling on each clade boundary, not a guess here.

### 3.3 The Tremellales `MAT` family produces confident false idiomorph calls outside Tremellales

Four of the five `idiomorph_gene_only` entries in this *Coprinopsis* genome are
Cryptococcus `MAT` calls with a **definite** idiomorph:

* `JAAGWA010000004.1:3525580-3545056` — `alpha`, `MFalpha1/2/3`, medium
* `JAAGWA010000012.1:107651-110349` — `alpha`, `MFalpha1/2/3`, medium
* `JAAGWA010000011.1:369126-392390` — `alpha`, `SXI1`, low
* `JAAGWA010000010.1:910864-911052` — **`a`**, `SXI2`, low

Three `alpha` and one `a` in the same haploid genome is self-contradictory. The
`MAT` family's `taxonomic_scope` is `[40410]`, one *Cryptococcus* species, and
`--phylum Basidiomycota` deliberately bypasses scope entirely. The idiomorph
call itself is doing what it is told; the fault is that a whole-phylum run
searches a one-species family against 3,174 genomes across 30+ orders. There is
a fifth, larger artefact of the same kind: a single 261 kb `MAT` cluster at
`JAAGWA010000001.1:1368880-1629971` holding **all nine** `MAT` genes — both
`SXI1` and `SXI2`, all three `MFalpha` and all three `MFa` — which lands on
`undetermined` only because both idiomorphs are indicated at once.

The narrow `taxonomic_scope` values are the root cause and they are worth
re-examining: every Basidiomycota family scopes to a single species taxid
(5334, 5346, 5270, 40410, 192523), whereas Mucoromycota's scopes are broad
clade taxids. Any Basidiomycota genome that is not one of those five species
falls to `phylum_fallback` and searches all nine families. For the 3,174-genome
library that is the normal case, not the exception.

---

## 3A. Per-order, per-locus self-consistency (genome-only)

Scored by coordinate overlap: did the run emit a `detected` entry for the
curated record's own family that overlaps the curated span, on the curated
contig? Genome-only, `--phylum Basidiomycota`, no proteome.

### Agaricales — 5 of 5 curated loci hit, every gene recovered

| genome | locus | curated span | detected span | conf | genes |
|---|---|---|---|---|---|
| *C. cinerea* Okayama-7 | HD | `JAAGWA010000001.1:1625680-1631546` | `1625683-1631394` | **high** | **2/2** |
| *C. cinerea* Okayama-7 | PR | `JAAGWA010000010.1:1806154-1826859` | `1806650-1826460` | **high** | **4/4** |
| *S. commune* H4-8 | Aalpha | `NW_026089539.1:1821007-1827456` | `1821010-1827453` | **high** | **2/2** |
| *S. commune* H4-8 | Balpha | `NW_026089548.1:241497-248430` | `190191-273652` | medium | **3/3** |
| *S. commune* H4-8 | Bbeta | `NW_026089548.1:216536-233781` | `190191-273652` | medium | **8/8** |

Every gene each record marks present was found in every case, including all
eight *S. commune* Bbeta pheromone genes. Aalpha lands 3 bp inside the curated
span on each side.

**The two 100%-complete calls are the two capped at `medium`.** Balpha found
3 of 3 and Bbeta 8 of 8, with `genes_missing=[]` in both — and both are
`medium`, purely because `assign_tier` still counts the unsearchable
`pheromone_receptor` alias as an unmet core requirement (§4.1). This is the
measured cost of that bug, and it is the worst possible case: it demotes
exactly the calls that are perfect.

**The Balpha/Bbeta sublocus structure is lost.** The two subloci are curated
7,716 bp apart (Bbeta ends 233781, Balpha starts 241497). Neither family
declares `max_cluster_gap_bp`, so both use the 25 kb default, which chains
them into one 83.5 kb cluster `190191-273652` — 26 kb of overhang below the
true start and 25 kb above the true end. Both families then report that same
merged interval. The genes are still correctly partitioned by family, so the
call is right and only the boundary is wrong, but the report shows one B
region where the curation describes two subloci.

This is the same class of finding as the Mucoromycota 50 kb ruling:
`max_cluster_gap_bp` is per-locus curation data and **no Basidiomycota locus
declares one**. 25 kb is demonstrably too wide for the *S. commune* B subloci.

### Tremellales — the genes are found at 100% identity, the locus call is not

*C. deneoformans* JEC21 (`GCF_000091045.1_ASM9104v1`, a MATalpha strain),
1 min 38 s, 68 detected entries, **zero at `high`**.

The real MATalpha locus **is** found. In the `MAT` family entry's own
`gene_evidence`:

| gene | location | identity | method | status |
|---|---|---|---|---|
| MFalpha1 | `NC_006686.1:1534710-1534823` | **100.0** | exonerate_refine | polished_agree |
| MFalpha2 | `NC_006686.1:1534710-1534823` | **100.0** | exonerate_refine | polished_agree |
| MFalpha3 | `NC_006686.1:1535437-1535550` | **100.0** | exonerate_refine | polished_agree |
| STE3 | `NC_006686.1:1577472-1578789` | **100.0** | exonerate_refine | polished_agree |

That is the genuine locus, at perfect identity, correctly polished. **It is
not reported as a Tremellales MAT locus.** Two things go wrong:

**(a) The `MAT` entry is a cross-chromosome chimera.** It carries
`fragmented: true` and three `segments`:

```
NC_006670.1:1974456-1987792
NC_006686.1:729509-1578789      (849 kb)
NC_006686.1:729509-1578789      (duplicate)
```

The top-level `contig`/`start`/`end` report **only segment 0**,
`NC_006670.1:1974456-1987792` — a different chromosome from the one holding
almost all of the evidence. A consumer reading those three fields, as the
GFF3 and the summary do, is handed a location that contains four of the
entry's ten gene-evidence rows' worth of nothing. The 849 kb second segment is
not a locus either.

Consequences that follow mechanically: `fragmented` downgrades the tier, so
the best MAT call in a genome whose MAT genes match at 100% is **`low`**. And
because a stray `SXI2` hit on another chromosome survives alongside `SXI1`,
`expected_genes_for_idiomorph` sees two idiomorphs, declines to narrow, and
the entry is `idiomorph: undetermined` — **for a strain that is definitively
MATalpha**. The 51 `idiomorph_resolutions` on the entry all resolve correctly
(MFalpha beats MFa, 100.0 vs 54.8, `idiomorph_margin` 42.6); the resolution
machinery is not what fails. What fails is that a locus was allowed to span
two chromosomes.

**(b) The correct interval is reported under the wrong family.** The span
`NC_006686.1:1534710-1578789` — exactly the real locus — is emitted as
`aLocus` (the *Ustilago* a locus family) at `medium`, found via `mfa1`/`pra1`.
The pheromone/receptor gene_class collision of §3.2 again.

**The mirror strain reproduces it exactly.** *C. deneoformans* JEC20
(`GCA_056621545.1_ASM5662154v1`, a MAT**a** strain) gives the same shape with
the idiomorphs swapped:

| gene | location | identity | method | status |
|---|---|---|---|---|
| MFa1 | `CM152755.1:1555944-1556069` | **100.0** | exonerate_refine | polished_agree |
| MFa3 | `CM152755.1:1555944-1556069` | **100.0** | exonerate_refine | polished_agree |
| MFa2 | `CM152755.1:1564451-1564576` | **100.0** | exonerate_refine | polished_disagree |
| STE3 | `CM152755.1:1567768-1569091` | **100.0** | exonerate_refine | polished_agree |

and the entry again reports `CM152752.1:1986764-2000100` — the wrong
chromosome — at `low`, `idiomorph: undetermined`, for a strain that is
definitively MATa.

**Tremellales result: 2 of 2 genomes find the true locus genes at 100%
identity, and 2 of 2 report it as a `low`-confidence, `undetermined`,
cross-chromosome chimera.** The detection works; the locus assembly and the
reporting do not.

`fragmented` cross-contig merging is the single highest-value thing to look at
for Tremellales. A cap — same contig only, or a maximum span — would turn
these two entries from `low`/`undetermined` into correct, idiomorph-resolved
calls.

### Ustilaginales — 2 of 2, both high

*M. maydis* 521 (`GCF_000328475.2_Umaydis521_2.0`), genome-only.

| locus | curated | detected | conf | genes |
|---|---|---|---|---|
| bLocus | `NC_026478.1:1690007-1693902` span curated on the same contig | `NC_026478.1:1690007-1693902` | **high** | **2/2** (`bE`, `bW`) |
| aLocus | `U37795.1` — a GenBank deposit, **not an assembly contig** | `NC_026482.1:1023525-1028134` | **high** | **3/3** (`mfa1`, `pra1`, `rba1`) |

The aLocus cannot be scored by coordinate overlap because its curated record
cites a locus-level GenBank accession rather than an assembly contig, so the
automated scorer files it `FAMILY-ONLY`. Read directly, the call is right:
all three curated genes at high confidence on `NC_026482.1`, which is *U.
maydis* chromosome 5, where the a locus lives. Counting it as found is a
judgement, not a coordinate match — stated so it is not mistaken for one.

### Per-order summary (genome-only, no proteome)

| order | genomes | curated loci | correctly reported | tiers |
|---|---|---|---|---|
| **Agaricales** | 3 | 6 | **6 / 6** | 4 high, 2 medium |
| **Ustilaginales** | 1 | 2 | **2 / 2** | 2 high |
| **Tremellales** | 2 | 2 | **0 / 2** | best MAT entry is `medium` and reports the *wrong idiomorph* |

**8 of 10 curated Basidiomycota loci are found and correctly reported without
a proteome.** Every gene of every hit locus was recovered — 6/6 loci in
Agaricales at 100% gene recall, including all eight *S. commune* Bbeta
pheromone genes.

The 2 failures are both Tremellales and both the same diagnosable cause
(cross-contig `fragmented` merging), not a homology or sensitivity failure:
the genes themselves matched at 100% identity in both genomes.

Caveat on the Tremellales rows: for JEC21, a MAT**alpha** strain, the
highest-confidence `MAT` entry reports `MFa1/MFa2/MFa3` — the **a**-idiomorph
pheromones — at `medium`. Anything consuming "best entry per family" gets the
opposite mating type from the truth.

### Cross-family noise scales with genome, not with truth

| genome | detected entries | distinct intervals | intervals claimed by >1 family | high | medium | low |
|---|---|---|---|---|---|---|
| *C. cinerea* | 65 | 45 | 8 | 3 | 8 | 54 |
| *S. commune* | 96 | 77 | 10 | 3 | 14 | 79 |

*S. commune* has 3 curated loci and produced 96 entries. As in *C. cinerea*,
`detection_pass` is `strict` for all 96 and `locus_class` is `mat_locus` for
81 — neither discriminates.

**The unreferenced Abeta family shows up as predicted (§1.2).** *S. commune*
NW_026089539.1 carries two high-confidence HD-class intervals:
`1821010-1827453` (the true Aalpha, also claimed by HD) and
`1245862-1269816` (claimed by HD high, Aalpha medium, bLocus medium). With no
curated Abeta record there is nothing to attribute a real Abeta sublocus to,
so any such locus can only surface under Aalpha or HD. Whether the second
interval is Abeta or an unrelated homeodomain pair **cannot be decided from
this data** — 550 kb from Aalpha is farther than a sublocus should sit, and
resolving it needs a curated Abeta record, not more inference.

The Tremellales `MAT` family again fires spuriously outside its clade: 17
`MAT` entries in *S. commune*, including seven `medium` calls with a definite
`a` idiomorph from `MFa1/2/3` triplets on six different contigs, plus one
`alpha`. Same root cause as §3.3.

---

## 4. Bugs found

### 4.1 `tiering.assign_tier` ignores `genes_not_searchable`

```python
core_genes = {g["name"] for g in expected_genes_for_idiomorph(family, score.genes_found)
              if g["role"] == "core_MAT"}
core_found = core_genes.issubset(set(score.genes_found))
```

`score` is a `FamilyScore` and carries `genes_not_searchable`, which
`assign_tier` never reads. A gene with no reference protein anywhere is
therefore counted as a core requirement the genome failed to meet.

`scoring.score_cluster` was fixed for exactly this and drops those genes from
the `fraction_found` denominator. `assign_tier` was not.

**Measured cost, on true calls.** In *S. commune* H4-8, the two curated
B-locus records are recovered completely and both are demoted:

| locus | genes_found | genes_missing | genes_not_searchable | tier |
|---|---|---|---|---|
| Balpha | `bar3`, `bap3-1`, `bap3-3` (**3/3**) | `[]` | `pheromone_receptor` | **medium** |
| Bbeta | all 8 `bbr2`/`bbp2-*` (**8/8**) | `[]` | `pheromone_receptor` | **medium** |

Nothing searchable was missed in either case, `fraction_found` is 1.0 for
both, and both are capped at `medium` solely by an alias gene that no curated
record anywhere in `db/` carries a protein for. These are the two most
complete locus calls in the run.

(Also seen as a harmless under-call on the *C. cinerea* Balpha false positive
at `JAAGWA010000010.1:1806650-1826460`.)

Interaction to handle in any fix: because no Basidiomycota family has a
`flanking_conserved` gene (§1.3), excluding unsearchable genes from
`core_genes` makes `high` reachable from a *single* found gene in a family
whose other core genes are all unsearchable. The `isolated_single_hit`
demotion currently only applies in the `not core_found` branch and would need
to apply here too.

### 4.2 `benchmark.match_ground_truth` taxid gate

See §2. Not a Basidiomycota-specific bug, but Basidiomycota is where it bites:
4 of 6 reference assemblies are unmatchable, including both tetrapolar
genomes.

---

## 5. Cost, and why no large sweep was run

Measured, one node, 4 cores, genome-only:

| genome | size | wall time | detected entries |
|---|---|---|---|
| *S. commune* `GCF_000143185.2` | 39 MB | **21 min 21 s** | 96 |
| *C. cinerea* `GCA_016772295.1` | 39 MB | **11 min 13 s** | 65 |
| *U. maydis* `GCF_000328475.2` | 20 MB | 2 min 10 s | — |
| *A. bisporus* `GCF_000300575.1` | 31 MB | 1 min 48 s | — |
| *C. deneoformans* JEC21 | 19 MB | 1 min 38 s | 68 |
| *C. deneoformans* JEC20 | 20 MB | 1 min 37 s | — |

**A 13x spread over six genomes of comparable size.** The cost is dominated by
the serial `exonerate` polish and scales with candidate count, not with genome
size — the two 39 MB Agaricomycete genomes take 11 and 21 minutes while a
31 MB one takes under two. The two slow genomes are precisely the two that are
closest to the short-pheromone reference records: *S. commune*'s own Bbeta
record contributes 7 `bbp2-*` precursors plus `bbp2_a`/`bbp2_b`, all short, and
`--phylum Basidiomycota` puts all of them in the query set for every genome.

Mean of these six is about 6.6 min. Agaricales alone (778 BFD genomes) is
therefore of order 85 core-hours if the mix resembles this sample, but the
sample is 6 genomes and the spread is 13x, so that figure should not be
planned against tightly. That is a SLURM batch, not an interactive run, and it
should not be launched before §4.1 and §2 are fixed and before the per-interval
collapse in §3.2 exists — otherwise it produces ~65 entries per genome with no
way to tell the 2 real ones from the 63 duplicates at scale.

---

## 6. Runs in this note

Outputs under `$SCRATCH/basidio/`, summaries only in the repo.
Reference command:

```
matpredict detect --genome <g>.fna [--proteins <g>.faa] \
  --phylum Basidiomycota --out-dir <out> \
  --evidence-diagnostics <out>/evidence_diagnostics.jsonl
```

Proteomes for *C. cinerea* (16,862 proteins) and *S. commune* (16,193) were
built from the NCBI FTP `*_protein.faa.gz` + `*_genomic.gff.gz` pair into
`search.PROTEOME_DEFLINE_FORMAT` (`>{protein_id} {contig}:{start}-{end}:{strand}`,
CDS envelope, strand from GFF column 7) and live at
`$SCRATCH/basidio/proteome/`. Every defline parses; every contig name is
present in the matching `.fna`. Note for the annotation-gap question: the
*C. cinerea* annotation **does** cover both MAT loci — 3 proteins inside the
curated HD window and 14 inside the curated PR window, several of them
195-210 bp, i.e. pheromone-precursor sized. Basidiomycota does not obviously
repeat the Mucoromycota annotation gap, at least here.
