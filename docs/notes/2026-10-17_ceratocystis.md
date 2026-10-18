# Ceratocystis fimbriata curation (2026-10-17)

Branch `ceratocystis-curation`. Curator decisions of 2026-10-17: (1) a new `system` value for unisexual species with
mating-type switching, (2) add the C. fimbriata records only after verifying the locus on a same-strain assembly (B10),
(3) 11 genomes are held out and were never used for a record, profile or tuning.

## 1. Literature (opened via PubMed/NCBI on 2026-10-17; abstracts and GenBank flat files only, not full texts)

* Wilken PM, Steenkamp ET, Wingfield MJ, de Beer ZW, Wingfield BD (2014) DNA loss at the Ceratocystis fimbriata mating locus
  results in self-sterility. **PLoS ONE** 9(3):e92180, PMID 24651494, doi 10.1371/journal.pone.0092180. (The brief called this
  a Fungal Genet Biol paper; PubMed and the GenBank reference lines both give PLoS ONE.) Abstract: three mating genes
  (MAT1-1-1, MAT1-2-1, MAT1-1-2) in an atypical organisation of the MAT locus of the type species; self-sterile offspring
  carry a 3,581 bp deletion that includes the entire MAT1-2-1 gene, flanked by two 260 bp direct repeats; the loss is permanent
  (uni-directional switching).
* Witthuhn RC, Harrington TC, Wingfield BD, Steimel JP, Wingfield MJ (2000) Deletion of the MAT-2 mating-type gene during
  uni-directional mating-type switching in Ceratocystis. Curr Genet 38:48-52, PMID 10953881, doi 10.1007/s002940000131.
  Abstract: MAT-2 HMG fragments are present in self-fertile isolates and self-fertile progeny of C. virescens, C. pinicola and
  C. fimbriata and absent in their self-sterile (MAT-1) progeny; C. eucalypti is strictly heterothallic.
* Harrington TC, McNew DL (1997) Self-fertility and uni-directional mating-type switching in Ceratocystis coerulescens.
  Curr Genet 32:52-59, PMID 9309170, doi 10.1007/s002940050247. NOTE: this paper is on **C. coerulescens**, not C. fimbriata;
  it shows only MAT-2 strains self and half of the selfed progeny are MAT-1, and hypothesises a MAT-2 gene deletion.
* Wilken et al. 2013, IMA Fungus 4:357-358 (PMID 24563841, doi 10.5598/imafungus.2013.04.02.14): draft nuclear genome of
  C. fimbriata, 29,410,862 bp. I did not establish which NCBI assembly it corresponds to.
* Not found/not opened: a Wingfield paper beyond the above; the Wilken et al. paper on the C. eucalypti locus was not looked up
  (the brief's "C. fimbriata/C. eucalypti mating-type locus" attribution is not something I could confirm).

## 2. What the two deposits contain (read from the GenBank records)

Both: Wilken et al. (PLoS ONE 2014), submitted 10-MAY-2013 (FABI, Pretoria), `/isolate="CBS114723"`,
`/type_material="culture from epitype of Ceratocystis fimbriata"`, Sanger + 454; they differ only in `/clone`.

| | KF033902.1 | KF033903.1 |
|---|---|---|
| clone | unswitched | switched |
| length | 46,416 bp | 42,835 bp (46,416 - 3,581) |
| MAT1-1-1 | yes (complement 27087-29023, AHV84685.2) | yes (27087-29023, AHV84693.2) |
| MAT1-1-2 | yes (complement 34661-36014, AHV84687.1) | yes (complement 31080-32710, AHV84694.2) |
| MAT1-2-1 | yes (30431-32100, AHV84686.1) | absent |
| MAT1-2-7 | yes (complement 29231-29883, WOK41352.1, /gene="MAT1-2-7") | absent (lies in the deleted interval) |
| flank (by product) | COX13, APN2, SLA2 (+APC5) | COX13, APN2, SLA2 at identical coordinates |

MAT1-2-7 is annotated only in the 2023 revision of KF033902; it is not named in the 2014 abstract.
Flank names rest on the depositors' product text (cytochrome c oxidase subunit VIa = COX13, DNA lyase = APN2, cytoskeleton
assembly control protein = SLA2), not on my own homology search.

## 3. Same-strain assembly (B10)

NCBI Assembly lists three assemblies of "Ceratocystis fimbriata CBS 114723" (BioSample SAMN01985854): GCA_000389695.1 (CFim_1.0),
.2 (Cfim_2.0) and .3 (Cfim3.0). Only .3 is in the BFD library (`GCA_000389695.3_Cfim3.0.fa.gz`); none is on the held-out list.
`scripts/find_record_assemblies.py verify` (job 30031365, Cfim3.0): both records pass,
unswitched 100.0% of the core segment, APWK03000076.1:59345-88276, 99.9% identity; switched 100.0%, :59345-88553, 99.88%.
Whole-deposit megablast: KF033902 maps colinearly over 52,199-98,750 at >=99.8% (46,552 vs 46,416 bp). KF033903 maps in two blocks with
the ~3.3 kb interval APWK03000076.1:81,476-84,795 uncovered: **Cfim3.0 holds the unswitched locus.** I did not run the `search`
subcommand (it sweeps every record without an assembly); the same-strain match rests on the NCBI Assembly strain field and the 99.9% blastn.
The other BFD genomes were not blasted against the deposits for this purpose (the 11 are held out; others are not CBS114723).

## 4. Decisions and assumptions

* `system: unisexual_switching` added to the schema, duckdb schema comment and tests (own commit). The only code enumerating
  system values is the schema; `build_duckdb.py` stores the string, and the holdout scripts copy it. No report text changed.
* Records: `5158_cbs114723_MAT_unswitched` ACCEPTED (db/Ascomycota/Microascales, idiomorphs MAT1-1 and MAT1-2, both gene sets
  present, `assembly_accession_basis: strain_match_verified`). `5158_cbs114723_MAT_switched` kept in db/candidates as
  `needs_review` (the schema has no `draft` status; needs_review is what existing drafts use). Reason: it passes the B10 rule
  mechanically, but the assembly carries the unswitched allele, so a self-call against it would see both idiomorphs, not MAT1-1.
* `MAT1-2-7` added to the Ascomycota MAT vocabulary in order.yml (core_MAT, MAT1-2, optional). APC5 omitted (not in vocabulary).
* validate: all 7 and 5 proteins match their records' own CDS translations at 100%/100%. `taxonomy_current` was first set by hand
  (the pre-PR-93 taxonkit call returned nothing). After merging main (PR #93: taxonkit failures are explicit errors) I re-ran
  `curate-db validate` on both records: on the login node it reported `taxonkit exited -9` (killed, no flag written); in a Slurm job
  (`cerato-validate`, 8 GB) the check passes and sets `taxonomy_current: true` for both. The value now comes from the check.
* Record ids use `5158_cbs114723_MAT_unswitched|switched` (no existing convention for a state; question for the curator).

## 4b. Curator rulings of 2026-10-17 (applied here)

1. The switched record stays in db/candidates as `needs_review` (not accepted).
2. The switching state stays in the record id (`..._MAT_unswitched`, `..._MAT_switched`).
3. Report labels stay neutral (MAT1-1 / MAT1-2); switching is mentioned only in record notes, no report text.
4. `idiomorphs: [MAT1-1, MAT1-2]` for the unswitched record.
5. MAT1-2-7 stays in the shared MAT vocabulary, `optional: true`.
6. The 11 genomes are scored on gene recovery and confidence only.

## 5. The 11 held-out genomes: NOT a recall test

**All 11 were already called on current main before the new record existed. The 2026-10-05 "uncalled" status is stale, so these genomes
no longer test whether the new records recover anything.** They are scored only on gene recovery and confidence. Table: `detect --taxid
<own taxid>`, MATPREDICT_HTML=0, main db (origin/main after PR #93) vs this branch's db, same merged code. Taxa: 5 C. fimbriata,
1 Huntiella omanensis, 4 Berkeleyomyces rouxiae, 1 C. atlantica.

| held-out genome (taxon) | before (main db) | after (new db) |
|---|---|---|
| GCA_000833645.1_ASM83364v1 (Huntiella omanensis) | no call | partial, MAT1-2, medium; found SLA2, MAT1-2-7; missing MAT1-2-1, APN2, COX13 |
| GCA_009914735.1_ASM991473v1 (C. fimbriata) | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, SLA2, MAT1-1-2; missing APN2, COX13 | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, SLA2, MAT1-1-2, MAT1-2-7; missing APN2, COX13 |
| GCA_017160715.1_Cfim0.1 (C. fimbriata) | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2, COX13 |
| GCA_026151695.1_PM34_polca_canu (Berkeleyomyces rouxiae) | mat, MAT1-1, medium; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-3 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-3 |
| GCA_032173455.1_Cf_hybrid (C. fimbriata) | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-4 | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-2-7, MAT1-1-4 |
| GCA_033035045.1_Cfim_AS236_NCSU (C. fimbriata) | mat, MAT1-2, medium; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-4 | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-2-7, MAT1-1-4 |
| GCA_036245355.1_LT1_01 (Berkeleyomyces rouxiae) | mat, MAT1-1, medium; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-3 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-3 |
| GCA_038382305.1_ASM3838230v1 (C. atlantica) | partial, MAT1-1, medium; found MAT1-1-1, SLA2; missing APN2, COX13 | partial, MAT1-1, medium; found MAT1-1-1, SLA2; missing APN2, COX13 |
| GCA_038406785.1_ASM3840678v1 (C. fimbriata) | mat, MAT1-2, medium; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2 | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-2-7 |
| GCA_048999915.1_ASM4899991v1 (Berkeleyomyces rouxiae) | mat, MAT1-1, medium; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-3 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-2, MAT1-1-3 |
| GCF_048999995.1_ASM4899999v1 (Berkeleyomyces rouxiae) | mat, MAT1-2, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13, MAT1-1-3 | mat, MAT1-2, high; found MAT1-2-1, APN2, SLA2, COX13, MAT1-1-3, MAT1-2-7; missing MAT1-1-1 |

Reading: labels never change where a call existed; confidence rises medium to high in 5 genomes; MAT1-2-7 is found in 6; Huntiella
omanensis goes from no call to a weak partial MAT1-2 call (SLA2 + a ~42%-identity MAT1-2-7 hit, MAT1-2-1 missing; the label
rests on that gene name). Several genomes carry MAT1-1-1 and MAT1-2-1 under one idiomorph label. GCF_048999995.1 loses MAT1-1-1: see 5b.

### 5b. Why GCF_048999995.1 (Berkeleyomyces rouxiae) loses MAT1-1-1

Five runs of this genome (`--taxid 2035830 --evidence-diagnostics`, same merged code, jobs cerato-inv): main db; branch db; branch db
with `--max-polished-clusters-per-family 0`; main db plus only the MAT1-2-7 vocabulary entry (= branch with the new record removed);
branch db with the old order.yml (new record, no vocabulary entry).

| run | MAT1-1-1 in genes_found | MAT1-2-7 |
|---|---|---|
| main | yes | no |
| branch | **no** | yes |
| branch, uncapped | no | yes |
| main + vocabulary entry only | yes | no |
| branch, record without the vocabulary entry | yes | no |

* **Not the polish cap.** Uncapped gives the same call; the `polish_scope` event shows 7 polish pairs attempted on the branch (6 on main),
  nothing dropped.
* **Cause: the two changes together.** The vocabulary entry alone changes nothing; the record alone (its MAT1-2-7 gene dropped as out of
  vocabulary) changes nothing; with both, the new record's MAT1-2-7 protein hits NW_027396493.1:1157712-1158366 at 62.1% identity, 100%
  coverage (bitscore 218), and the evidence log gains an `idiomorph_resolution` row where MAT1-2-7 (60.2%) beats MAT1-1-1 (32.1%) at
  100% overlap.
* **The lost "MAT1-1-1" was a weak fragment, identical in both runs:** 1157973-1158140, 32.1% identity, 18% coverage, bitscore 33.9
  (floor 33.0), from the ipo323 MAT1-1-1 reference, lying inside the region that is a MAT1-2-7 homolog, on a locus the call labels MAT1-2
  (HMG gene MAT1-2-1 present; MAT1-1-1 had already lost an identity resolution to MAT1-2-1 on main, 33.3 vs 37.4). On the evidence this
  is the gene being named correctly, not a lost real MAT1-1-1. Not checked: a manual alignment of the region or a look at the annotation
  of that genome.
* The polish-cap fix would not cover this, since the cap is not involved.
* Recommendation (not a merge): the change is acceptable as it is; the curator may prefer to confirm by eye first.

## 6. Recall test

**No recall test is currently possible**: every genome that could have tested the new record (the 11 held-out genomes, and the 9 other
strong-linked uncalled Microascales in results/2026-10-18_uncalled_triage on origin/uncalled-triage: 8 Scedosporium apiospermum and 1
Thielaviopsis ethacetica) is already called on current main, with high confidence, and the call is identical after the new record
(table below, the 9 genomes). This will be re-opened if the re-run of the 63 triage-linked genomes on current main finds uncalled
Microascales genomes.

| held-out genome (taxon) | before (main db) | after (new db) |
|---|---|---|
| GCA_001599055.1_JCM_6961_assembly_v001 (Thielaviopsis ethacetica) | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2, COX13 |
| GCA_002158515.1_SaHDO1_v1 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 |
| GCA_040285505.1_ASM4028550v1 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2; missing COX13 | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2; missing COX13 |
| GCA_052249755.1_ASM5224975v1 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 |
| GCA_058095615.1_UMB_ScApio477_1.0 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 |
| GCA_058095685.1_UMB_ScApio478_1.0 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 |
| GCA_058131905.1_UMB_ScApio476_1.0 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 |
| GCA_982397395.1_T61-F (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, APN2, SLA2, COX13 |
| GCF_000732125.1_ScApio1.0 (Scedosporium apiospermum) | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 | mat, MAT1-1, high; found MAT1-1-1, MAT1-2-1, APN2, SLA2, COX13 |


## 7. Tests

HPCC pixi env, `pytest tests/db tests/detect --continue-on-collection-errors`, after merging origin/main: this branch 31 failed, 1078 passed,
4 skipped, 2 errors; clean origin/main 31 failed, 1072 passed, 4 skipped, 2 errors. Same failures (classifier/paralog tests); the 2 errors
are tests/db/test_draw.py and test_synteny.py (pygenomeviz not installed). The +6 passes are the new schema/duckdb tests.
