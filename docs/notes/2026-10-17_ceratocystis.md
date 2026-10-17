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
* validate: all 7 and 5 proteins match their records' own CDS translations at 100%/100%. `taxonomy_current` came back false
  because taxonkit given a taxid with no trailing newline returns nothing in this environment (pre-existing; `taxonkit` on 5158
  gives Microascales/Ceratocystidaceae/Ceratocystis), so I set it true by hand.
* Record ids use `5158_cbs114723_MAT_unswitched|switched` (no existing convention for a state; question for the curator).

## 5. Held-out genomes, before (main db) and after (this db)

`detect --taxid T`, MATPREDICT_HTML=0, jobs `cerato-detect` (30031408). Taxids from NCBI Assembly. The 11 genomes are
not all C. fimbriata: 5 are C. fimbriata, 1 Huntiella omanensis, 4 Berkeleyomyces rouxiae (taxid 2035830), 1 C. atlantica.
On current main they are **already called** (the 2026-10-05 "uncalled" status is stale), by other families' records.

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

Reading: the new record is a reference_record in every "after" call. Idiomorph labels do not change where a call existed before.
Confidence rises medium to high in 5 genomes, MAT1-2-7 is found in 6. Huntiella omanensis goes from no call to a partial MAT1-2
locus built on SLA2 and a 42%-identity MAT1-2-7 hit, MAT1-2-1 missing; that label rests on the gene name of the 42% hit and is weak.
GCF_048999995.1 (B. rouxiae) loses MAT1-1-1 in the after run (found before). Most "found" gene lists include MAT1-1-1 and MAT1-2-1 together
yet the idiomorph is a single label; the report does not say "unisexual". Not a recall measurement: before/after differ in whether
a call exists at all for only one genome, because main already calls the rest.

## 6. Tests

HPCC pixi env, `pytest tests/db tests/detect --continue-on-collection-errors`: this branch 31 failed, 1035 passed, 2 errors;
clean origin/main 31 failed, 1029 passed, 2 errors (same failures: tests/detect/test_idiomorph_classifier.py, test_paralog_class.py;
collection errors tests/db/test_draw.py, test_synteny.py need pygenomeviz). The +6 passes are the new schema/duckdb tests.
