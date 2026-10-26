# Basidiomycota, v0.7.0 definitive run (`cleanrun-2026-10c`): independent judgement

> **Provenance.** Independent second read written by a separate Opus agent on 2026-10-26, read-only, from the committed tables, the
> cluster reports and the control reruns. The scripts and intermediate tables it names (`parse.py`, `a*.py`, `p070_*.tsv`) were in a
> session scratchpad and are not in the repository. Checked afterwards against the tables by the session that ran the campaign:
> the 7 *Austropuccinia* / *Phakopsora* genomes (called in v0.6.0, now "no report"), the 794-of-1,191 core-gene share of defective
> flags, and the conflicting controls for GCA_011057475.1 (cap-off reproduces the first run, the first run's own code does not).
> Not re-checked: the 102 new PR loci, the 6 long-intron models, the 198 high HD loci on one region, the empty `two_idiomorphs`.
> The changelog wording was corrected for the three checked points. The curator list below is open.


Reviewer: second read (Opus), 2026-10-26. Read-only. Inputs: `origin/main` copies of `results/2026-10-26_v070_basidiomycota/`
(genomes.tsv, loci.tsv, by_class.tsv, reports_all.tar.zst), `results/2026-10-03_basidiomycota_v060/` (tables and tarball),
all 3,254 clean-run 10b reports (per-chunk `reports.tar.zst` from the cluster), the cap5k and old10b control reruns, and the
cluster run logs and suppress lists. Local `main` is 37 commits behind `origin/main`; files were read from `origin/main` via `git archive`.
Scripts and intermediate TSVs are in the scratchpad (`parse.py`, `a*.py`, `cmpall.py`, `cmp36.py`, `p070_*.tsv`, `p060_*.tsv`, `p10b_*.tsv`).

## Verdict: pass with a curator list

The run itself is clean: 3,254 reports, every one stamped `cleanrun-2026-10c`, genetic code 1, no genome errors or anomalies in
the 16 rollout summaries, and the 21 missing reports are the suppress list. The changes against the earlier clean run (10b) are
explained by the bLocus intron cap and by tie-break variance, and I reproduce the 5 call and 26 model differences (plus one
label-only difference not on the list). The changes against v0.6.0 are mostly intended (the 5 kb cap removing chimeric models,
new records). None of what I found is a run failure. But several things in the tables are not stated in the CHANGELOG and a person
should decide on them before the numbers are quoted: 7 rust genomes that v0.6.0 called are dropped by an annotation-pipeline
suppress list; 102 new pheromone-receptor (PR) loci come from one record's pheromone proteins; 6 bLocus models (not 1) carry
5-10 kb introns; and a large share of "high" HD loci rest on one homeodomain fragment that two references both hit.

## 1. Rates by class (my counts)

Source: `a25.py` (genomes.tsv v070 and v060; 10b `tables/genomes.tsv` from the cluster).

| class | listed | no report | called v0.7.0 | called 10b | called v0.6.0 |
|---|---|---|---|---|---|
| Agaricomycetes | 1856 | 5 | 1798 | 1797 | 1798 |
| Tremellomycetes | 559 | 1 | 429 | 427 | 429 |
| Microbotryomycetes | 318 | 0 | 273 | 273 | 273 |
| Ustilaginomycetes | 133 | 0 | 133 | 133 | 133 |
| Pucciniomycetes | 131 | 14 | 104 | 104 | 107 |
| Exobasidiomycetes | 92 | 0 | 87 | 87 | 88 |
| Malasseziomycetes | 85 | 1 | 83 | 83 | 83 |
| Wallemiomycetes | 51 | 0 | 51 | 51 | 51 |
| Cystobasidiomycetes | 28 | 0 | 7 | 7 | 7 |
| 11 small classes | 22 | 0 | 12 | 12 | 13 |
| **total** | 3275 | 21 | **2977** | 2974 | 2982 |

Status totals match the CHANGELOG (2,977 / 277 / 21; loci 5,986 = 3,802 high / 2,173 medium / 11 low). Lowest orders (n >= 10):
Kriegeriales 2/10, Microbotryales 9/38, Cystobasidiales 6/18, Filobasidiales 19/41, Trichosporonales 59/127, Cystofilobasidiales 21/42;
all within one genome of v0.6.0. Genera with no call: *Austropuccinia* (0/5, all no report; v0.6.0 called 4), *Pascua* 0/7, *Tausonia* 0/6,
*Erythrobasidium* 0/4, *Phenoliferia* 0/4 (the last four were also uncalled in v0.6.0).

## 2. Lost and new against v0.6.0 (`change_vs_baseline`, checked report by report)

Category counts (genomes.tsv): identical 2,785; uncalled both 271; locus count changed 102; genes_found only 41; span >1 kb 37;
no new report 21; NEW 8; LOST 6; confidence DOWN 2; class changed 2.

| genome | species | v0.6.0 call | now | my reading of the reports |
|---|---|---|---|---|
| GCA_001600235.1, GCA_002897395.1 | *Vanrija humicola* | bLocus | uncalled | v0.6.0 bE/bW 26-50 kb models with 23-26 kb "introns", 26-28% identity: chimeric. Real loss of a false call. |
| GCA_049862275.1 | *Parajaminaea* sp. | bLocus | uncalled | same pattern (per intron-cap README; not re-read by me) |
| GCA_058224995.1 | *Phyllozyma linderae* | HD | uncalled | v0.6.0 HD2 model 79 kb with a 78 kb intron, 31% identity: chimeric. Now HD2+MIP1 withheld at 1 modelled gene. |
| GCA_900009495.1 | *Microbotryum violaceum* | PR | uncalled | v0.6.0 receptor 10.9 kb, 8 kb intron, 23.7% identity: chimeric. Two bLocus clusters withheld with 1 modelled gene each. |
| GCA_033819875.2 | Microbotryomycetes sp. JL221 | PR | uncalled | v0.6.0 receptor 39% identity + CAAX ORF. Now three PR clusters withheld with 1 modelled gene; HD and bLocus also 1 each. Detector miss, not a real loss (inferred; I did not find which model dropped). |
| GCA_003057935.1, 030280995.1, 036250445.1, 054948715.1 | *Hemileia vastatrix* | uncalled | rustHD high | new record (#116); GCA_030280995 bE/bW 99-100% identity, 100% coverage. Sound. |
| GCA_022506105.1, GCA_900008855.1 | *Microbotryum violaceum* | uncalled | PR medium | receptor 37% identity, 19% (or 54%) coverage + a CAAX ORF; weak. |
| GCA_982397355.1 | *Trichosporon asahii* | uncalled | PR medium | receptor 40% identity, 29% coverage + CAAX ORF; weak. |
| GCA_002222225.1 | *C. neoformans* | not in v0.6.0 | MAT high | not in the v0.6.0 table; FAO1/SXI1/STE3/MFalpha 83-95% identity. Sound. |

The 7 Pucciniales genomes that v0.6.0 called and that now have no report (*Austropuccinia psidii* x4, *Phakopsora pachyrhizi* x3; all
`large_000`) are counted as "no new report", not as lost. The large-job log says "suppressed" and the reasons in the BFD list are about
annotation: "Rust genome EVM never completes", "same EVM failure", "funannotate_too_few_models:15", "Temp suppress as it is too large"
(`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt`). MAT detection runs on the assembly and v0.6.0
called rustHD in all 7 (0.56-1.27 Gb assemblies). The other 14 no-report genomes are tiny or fragmented assemblies (30 kb to 2.4 Mb,
or funannotate training failures); none was called in v0.6.0. The clean-run plan's statement that GCA_002222225.1 is on the suppress
list is out of date: it is not on either list now and it has a report.

Locus-level, family-aware matching (`a9.py`, genomes called in both runs): 102 new PR loci, 30 PR / 18 HD / 5 bLocus spans changed >1 kb,
13 Aalpha -> HD relabels, 6 PR and 2 bLocus loci gone, 1 confidence change. Causes I traced:
- **102 new PR loci (94 genomes, 87 in Agaricales).** All 102 were withheld in v0.6.0 (99 `below_fraction_floor`, 3 `modelled_gene_bar`) and
  all 102 now carry a `fungal_mating_type_pheromone` tblastn hit from `5346_a43-b43-okayama-7_PR_B43` (`a10.py`, `a13.py`). That record gained
  its accession-less pheromone proteins in the build-gff fix (#121/#134); its pheromone hits went from 92 to 362 entries. All 102 are medium.
  The CHANGELOG's Basidiomycota line ("6 lost and 8 new, 2 confidence down") does not mention them. The same record also pushed 2 Megacollybia
  PR loci below the fraction floor (inferred mechanism: larger denominator).
- **13 Aalpha -> HD** (10 *Schizophyllum commune*, 3 *Termitomyces*): the v0.6.0 Z models had 5.3-9.4 kb introns (chimeric); the 5 kb cap removed
  them, so the locus loses Y/Z and is labelled HD. Correct direction, but *Schizophyllum* now gets the less specific label.
- **2 bLocus loci gone** (GCA_008065275.1, GCA_054471085.1): v0.6.0 bE models 52 kb / 25.7 kb with 27.6 / 25.4 kb introns. Correct.
- **Confidence:** only one locus moved (GCA_002718315.1 HD medium -> low: HD1 and HD2 are 6% and 14% coverage tblastn fragments; the call rests on
  MIP1 and beta_fg). The second "confidence DOWN" genome (GCA_019142975.1) is an Aalpha -> HD relabel. `tier_note` is empty for all 5,986 loci.
- `change_vs_baseline` puts these relabels under "span changed", "class changed", "confidence DOWN" or "genes_found changed only"; there is no single category for them.

## 3. Differences against the earlier clean run (10b) and their causes

My own comparison of all 3,254 report pairs (`cmpall.py`; call = family, coordinates, idiomorph, confidence, polished_genes, genes_found;
model = gene coordinates, exons, method, reference, status): 3,175 identical, **5 call**, **27 model**, 47 differ only in withheld (suppressed) loci.
The 5 call genomes equal `diff_basidio_call.txt`; 26 of my 27 equal `diff_basidio_model.txt`; the extra one is GCA_928852645.1 (*Mycosarcoma
mackinlayi*), where one bE model has the same coordinates and identity but `method` changes from `exonerate_refine` to `exonerate_unrefined`.
All 98 withheld-locus differences in the 47 are bLocus (cap), and none reaches the tables.

Controls (`cmp36.py`; same signature):

| pattern (new vs 10b / cap5k vs 10b / old-code rerun vs 10b) | genomes | my reading |
|---|---|---|
| differ / same / same | 4 call + 14 model | cap. Includes 10 bLocus model changes and 4 *Cryptococcus* MFa copy switches (see below) |
| differ / differ (= new) / differ (= new) | 1 call + 10 model | variance: the old code on rerun gives the new result |
| differ / same / differ (= new) | 1 (GCA_011057475.1) | variance, not cap: the old code reproduces the new result |
| differ / differ (= new) / same | 1 (GCA_020975605.1) | ambiguous: variance, or a non-cap change in 10b..10c (the pair-loser narrowing #140 is also in that range) |

The CHANGELOG says 19 cap and 12 variance; I count 18 clearly cap, 12 clearly variance and 1 ambiguous, so the split is off by one at most
(inferred: GCA_011057475 was counted as cap because cap5k matched 10b). Both odd ones are *Cryptococcus* MFa/MFalpha copy switches between
identical pheromone copies (90% / 87% identity, same score), which is the textbook equal-score tie.

Four of the "cap" genomes (GCA_023650555.1, 056632815.1, 056632945.1, 057715595.1) change only which MFa copy is reported, in the MAT family,
which keeps the 5 kb cap. One rerun each returned the 10b copy at 5 kb, so the cap code path changed the tie-break order rather than the MAT
alignment itself (inferred; one rerun per condition cannot separate this from chance). The one variance call change, GCF_022376455.1
(*Epithele typhae*), adds a low HD locus resting on MIP1 and beta_fg with a 9%-coverage HD2 fragment; the old code makes the same call on rerun.

Confidence in this account: high for the 4 restored b-locus calls and the 10 bLocus model changes (direct intron evidence plus controls);
moderate for the variance group (one rerun each). What would change my mind: a second cap5k rerun that differs from the first on the MFa
genomes would move them to variance; a rerun of 10c code that still differs from 10c in the variance group would point to a non-deterministic
step outside polishing.

**Long introns let in by the 10 kb bLocus cap.** In v0.7.0, 6 polished bLocus models in 5 genomes have an intron over 5 kb (none in 10b; 0
in any other family): *Pucciniomycotina* sp. GCA_058411105.1 bW 9.6 kb and 6.5 kb introns, *Quambalaria eucalypti* GCA_004016185.1 bW 7.7 kb,
*Jaminaea rosea* GCF_003144245.1 bE 8.5 kb, *Microstroma glucosiphilum* GCF_003144135.1 bE 5.7 kb, *Tilletia fusca* GCA_045933405.1 bE 5.1 kb
(`a16.py`). All are 12-33% coverage, all in loci already called high in 10b; models over 2 kb intron in bLocus went from 11 to 23. The
CHANGELOG mentions one such model. The three restored calls rest on bE 10-13% coverage (114-236 bp) and bW 17-26% coverage with 2.2-3.0 kb
introns, and are reported high.

## 4. Idiomorph and allele sanity

| family | loci | idiomorph values |
|---|---|---|
| HD 2529, PR 1569, bLocus 348, redHD 253, aLocus 185, rustHD 172, Aalpha 168, Bbeta 33, Balpha 29 | 5286 | all `undetermined` (no allele typing for b/HD/PR) |
| MAT (Tremellomycetes type) | 399 | alpha 271, a 68, undetermined 60 |
| redPR | 250 | A2 191, A1 59 |
| wallMAT | 51 | v1 33, v2 18 |

- 164 of the 339 MAT a/alpha labels use the identity fallback margin (#137, known); 9 labels are `label_unsupported` (bit margin <20), 7 of them in
  *Mrakia*/Tremellales with margins of 0.8-2.7 bits.
- Genomes holding both labels: 5 with a + alpha (two *C. neoformans x deneoformans* AD hybrids, plausible; two *C. neoformans* GCA_025531435.1 and
  GCA_025531515.1, where the a locus is SXI2 at 100% identity to JEC20; *Mrakia blollopis* with two clusters 25 kb apart on one contig at 1-3 bit
  margins) and 4 with A1 + A2 redPR. `two_idiomorphs` is empty in all 3,254 reports, so the report-level flag never fires for these.
- `homothallic_candidate`: 0 loci (the policy targets alpha-box/HMG pairs, so none is expected here).
- MAT family outside Tremellomycetes: 14 loci (*Ceratobasidium* 7, *Thecaphora* 4, *Exobasidium*, *Neolentinus*, *Microbotryum scorzonerae* 1 each),
  all on STE3 plus PAN6 or RPL39. GCA_047651735.1 is high on PAN6 at 10.9% coverage and STE3 at 37.6%.

## 5. Defective-gene flag (`a5.py`, `a6.py`)

1,038 of 5,986 loci (847 genomes) carry at least one defective gene; 1,191 of 19,243 gene entries. 1,144 are `exonerate_refine` models
(7.8% of that method), 47 `miniprot_refine` (25% of 188). 546 frameshift only, 645 internal stop only, none both.
Unlike the Ascomycota wording in the CHANGELOG ("mostly low-identity flank-gene models"), here **794 of 1,191 are core genes**: HD1 337,
beta_fg 270, HD2 199, MIP1 114, pheromone_receptor 83, bW 50. The rate is not concentrated at low identity (5.6% at 30-40%, 9.6% at 40-50%,
3.3% above 90%). By order: Wallemiales 43%, Boletales 36%, Russulales 31%, Hymenochaetales 26% of loci, against 0-8% in most yeast orders.
Boletales HD1 237/978 flagged (median identity 45%, coverage ~55%, 1 stop or 2 frameshifts). *Wallemia* STE3: 18/51 genomes carry the same
model (64.7% identity, 97.5% coverage, 4 exons) with one internal stop, identical across isolates; that pattern fits either a real feature of
these isolates or a systematic splice-site mis-model, and I cannot tell which from the reports.

## 6. Other checks

- Locus spans: 8 loci are `span_exceeds_plausible_bound` (all MAT, 202-269 kb). The largest, GCA_006352305.1, ends at a 34.7%-identity, 11%-coverage
  SXI2 tblastn hit 85 kb past the next gene. No exact duplicate loci and no overlapping same-family loci. 15 span groups are called under two or
  three families with identical coordinates (11 MAT+PR, 1 HD+bLocus, others); counted as separate loci in the 5,986.
- Overlapping core models: in 1,063 of 2,529 HD loci the HD1 and HD2 models overlap by at least half the shorter one (one homeodomain region hit by
  both references), and 356 PR loci count a CAAX ORF and a pheromone tblastn hit on the same ORF. **198 high HD loci have one distinct core region
  and no polished flank gene** (e.g. GCA_019142975.1, GCA_003316405.1: HD1 and HD2 at the same 200-650 bp). This is v0.6.0 behaviour, not a change.
- 240 genomes are called only through PR (all medium; 20 with best receptor coverage under 30%). 1,447 of 1,569 PR loci are `unverified`.
- 49 genomes report an assembly gap at a locus; 24 loci span more than one segment.
- The README validation row (Pucciniomycotina 13.2% called) does not match this run (Pucciniomycetes 104/117 with reports).

## 7. The 277 uncalled genomes (`a19.py`)

Classes: Tremellomycetes 129, Agaricomycetes 53, Microbotryomycetes 45, Cystobasidiomycetes 21, Pucciniomycetes 13, others 16. 208 are
phylum-fallback. No uncalled genome has a withheld cluster with 2 or more modelled genes: 239 have a best cluster with exactly 1 (reasons over
all their withheld clusters: 899 `modelled_gene_bar`, 193 `below_fraction_floor`), 38 have none. Uncalled Tremellomycetes, Microbotryomycetes
and Cystobasidiomycetes have better assemblies than the called ones (median N50 0.5, 1.3 and 1.1 Mb), so these are reference or detector
gaps, not assembly problems. Uncalled Agaricomycetes and Pucciniomycetes are fragmented (median N50 9.0 kb and 4.9 kb). Under the rules as
set (strong-single-gene off) none of them should have been called; whether a 1-gene bLocus or HD cluster in a well-assembled yeast should be
is the open #97 question, not a run fault.

## 8. Reports read (YAML, gene evidence and locus fields)

Changed: GCA_001444555.1 *Naganishia albida* (whole file), GCA_004016185.1, GCF_022376455.1, GCA_002718315.1, GCA_019142975.1, GCA_051942515.1
*Mycena purpureofusca*, GCA_022506105.1, GCA_030280995.1; withheld and not-detected sections of GCA_033819875.2, GCA_058224995.1, GCA_900009495.1,
GCA_001600235.1; gene evidence across 10b / 10c / cap5k / old10b of GCA_011057475.1, GCA_020975605.1, GCA_928852645.1 and the 18 cap genomes.
Unchanged: GCA_016772295.1 *Coprinopsis cinerea* (HD1/HD2, receptor and B43/B44 pheromones at 100%; sound), GCA_004917895.1 *Wallemia*,
GCA_003316145.1 *Neoboletus*, GCA_003316405.1 *Thelephora*, GCA_025531515.1 *C. neoformans*, GCA_047651735.1 *Ceratobasidium*, GCA_006352305.1, GCA_000950635.1 *Mrakia blollopis*.
Reports carry no confidence-reason field; I judged confidence from polished_genes, coverage and identity.

## Curator list

1. **Rust genomes on the BFD suppress list.** 7 genomes called in v0.6.0 (4 *Austropuccinia psidii*, 3 *Phakopsora pachyrhizi*) are skipped for
   annotation-pipeline reasons, and so are a fifth *A. psidii* and *Hemileia* GCA_004125335.1 ("too large"), both uncalled in v0.6.0; no
   *A. psidii* genome has a report. Decide whether MATPredict should keep them out (and say
   so in the CHANGELOG as 7 lost calls) or move them to a MATPredict-only allow list and run them.
2. **102 new PR loci from the B43 pheromone proteins.** Accept as intended (and document), or require that the pheromone hit and the CAAX ORF are
   not the same ORF before they count as two genes (356 loci do this).
3. **bLocus 10 kb cap.** 6 models with 5.1-9.6 kb introns in 5 genomes, not 1. Accept with a known-issue line that names all 6, or withhold models
   over 5 kb from GFF/locus exports while keeping the calls.
4. **High confidence on one fragment.** 198 HD loci are high with one homeodomain region and no flank; the 3 restored b-locus calls are high on
   10-13% bE coverage. Decide whether "high" should require two distinct core regions.
5. **Defective flag text.** In Basidiomycota it is mostly core genes (HD1, HD2, receptor) in Boletales/Agaricales, and a fixed STE3 stop in *Wallemia*.
   Adjust the CHANGELOG wording for this phylum; decide whether to check the *Wallemia* STE3 stop against the sequence.
6. **Mixed idiomorphs.** 5 genomes with a + alpha and 4 with A1 + A2, while `two_idiomorphs` is always empty. Decide whether the flag should cover
   Basidiomycota MAT and redPR, and look at the two *C. neoformans* and the *Mrakia blollopis* cases.
7. **MAT family outside Tremellomycetes** (14 loci on STE3 + PAN6/RPL39). Keep, demote, or label as unverified.
8. **Run-to-run attribution.** Correct the 19/12 split to 18 cap / 12 variance / 1 ambiguous (or say why GCA_011057475.1 is cap), and add
   GCA_928852645.1 (method label only) to the diff list or say why the comparator ignores `method`.
9. **README validation row** for Basidiomycota is from an older campaign; update or date it.

## Not verified

- Whether any new or lost call is biologically right: no truth set for Basidiomycota was used; judgements of "chimeric" rest on intron length and
  identity only.
- The *Wallemia* STE3 stop and the Boletales HD1 defects against the genome sequence (I did not open any FASTA).
- Loci on very small contigs: reports give a distance to one contig edge, not contig length, so I could not count loci on tiny contigs.
- The cause of the JL221 loss (which model stopped counting) and of the tie-break change in the 4 *Cryptococcus* MFa genomes (inferred, one rerun each).
- The *Parajaminaea* loss (taken from `results/2026-10-25_intron_cap/`, not re-read).
- Mucoromycotina and Ascomycota were out of scope; the 5 Ascomycota run-to-run genomes were not checked.
- I did not rerun anything; all statements about code effects come from existing control runs and the git log between the tags.
