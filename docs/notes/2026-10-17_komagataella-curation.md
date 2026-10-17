# Komagataella curation, 2026-10-17

Branch `komagataella-curation` (from `origin/main` 5b9a9e0). Nothing is merged and no PR is open. Curator decision of 2026-10-17: curate a MATa record from *Komagataella phaffii* CBS 7435 and a MATalpha record from GS115, each anchored on a same-strain assembly. This note records the evidence, the family decision, the before/after detection runs, and what could not be decided or verified.

## What changed

| File | Change |
|---|---|
| `db/Ascomycota/order.yml` | New family `MATkp` appended at the end of `loci` (scope genus *Komagataella* 460517, `fallback_searchable: false`, idiomorphs `MATa`/`MATalpha`, genes MATa1, MATa2, MATalpha1, MATalpha2, optional sla2 and dic1). Appended at the end, not beside `MATyl`, to keep the merge with `dipodascales-curation` (which edits `MATyl`/`MTL`) trivial. |
| `db/Ascomycota/Pichiales/460519_cbs-7435_MATkp_MATa/` | New record (FR839631.1:1666-7892). |
| `db/Ascomycota/Pichiales/460519_gs115_MATkp_MATalpha/` | New record (NC_012966.1:1646885-1652702). |
| `tests/detect/test_fallback_searchable.py` | The test pinned the exact set of scope-only families; `MATkp` added to it (renamed `..._scope_only_families`). `dipodascales-curation` edits the same test (adds MATyl/MATdc); expect a one-line merge conflict. |

`propose`, `validate`, `accept`, `build-gff` were run on the HPCC (on a compute node: on the login node `taxonkit` was killed, "taxonkit exited -9", and validate then leaves `taxonomy_current: null` with a `taxonomy_error`, as PR #93 intends; on a compute node `taxonomy_current: true` for both). Both records: every gene `pass`, 100% identity, 100% coverage.

## Evidence

### The locus (all sequences re-fetched from NCBI, 2026-10-17)

K. phaffii has two complete MAT cassettes on chromosome 4, one near the chromosome end and one 135 kb in, in inverted orientation within a 138 kb region; inversion switches mating type. Opened and read: Hanson, Byrne & Wolfe 2014 (PMID 25349420, doi 10.1073/pnas.1416014111, abstract), Heistinger et al. 2018 (PMID 29061733, doi 10.1128/MCB.00398-17, full text from PMC5748462), Feng et al. 2020 (PMID 31756769, abstract), Coughlan et al. 2016 (PMID 27497317, abstract).

| Strain / accession | Region | Layout (strand) |
|---|---|---|
| CBS 7435, FR839631.1 (GCA_000223565.1) | 1666-7892 | DIC1-1 (+) 1666-2583; MATa2 (-) 2740-3276; MATa1 (+) 4058-4735; SLA2-2 (-) 4755-7892 |
| CBS 7435, FR839631.1 | 135891-141707 | SLA2-1 (+); MATalpha1 (-) 139129-139701; MATalpha2 (-) exons 140525-140535 + 139957-140467; DIC1-2 (-) 140790-141707 |
| CBS 7435, LT962479.2 (GCA_900235035.2) | 1770-7997 and 136012-141828 | same, +105 bp offset |
| GS115, NC_012966.1 (GCF_000027005.1) | 1646885-1652702 | DIC1 (+); MATalpha2 (+) exons 1648057-1648067 + 1648125-1648635; MATalpha1 (+) 1648891-1649464; SLA2 (-) 1649568-1652702 |

The earlier agent's layout (MATa near 2,732-5,018, other configuration 136-142 kb, GS115 1,646,885-1,652,702) is confirmed by annotation. Corrections: `FR839631` is not a small deposit but the complete chromosome 4 of CBS 7435 (1,820,458 bp, the Küberl 2011 assembly, PMID 21575661); the GS115 cassette is on the reverse strand relative to CBS 7435 (it lies 126 kb from the end of NC_012966.1, i.e. at locus 2). The flanks DIC1/SLA2 are an inverted repeat shared by both cassettes (blastn: 99.97% over 3.1 kb of SLA2, 100% over 1 kb of DIC1), so only the MATa1/MATa2 or MATalpha1/MATalpha2 interval distinguishes the cassettes.

In CBS 7435 the MATa cassette is at the chromosome end (locus 1) and MATalpha at locus 2; GS115's assembly has only the locus-2 MATalpha cassette. blastn of the CBS 7435 MATa cassette against NC_012966.1 hits only the SLA2/DIC1 repeats at that cassette, no MATa1/MATa2 anywhere. Not verified: whether GS115 lacks the other cassette or only the assembly does (the NC_012966.1 record says finishing was in progress).

### B10 anchor (strain match, then blastn)

- CBS 7435: strain label CBS 7435 on both assemblies (BioSample SAMEA2272267 for GCA_000223565.1, which holds FR839631.1; GCA_900235035.2 organism "K. phaffii CBS 7435"). blastn of FR839631.1:1665-7892 (6,228 bp) against LT962479.2: 100.000% over 6,228/6,228 bp at LT962479.2:1770-7997. Recorded `assembly_accession: GCA_900235035.2`, `strain_match_verified`, `assembly_location` LT962479.2:1770-7997, identity 100.0. The segment itself is the deposit FR839631.1.
- GS115: the segment is NC_012966.1 (RefSeq of GCF_000027005.1), `assembly_accession_basis: sequence_link`. The GS115 cassette is 99.98% identical over 5,818 bp to the CBS 7435 MATalpha cassette (one extra base, below).

### Gene models and introns

- MATa2 (178 aa), MATa1 (225 aa), MATalpha1 (190 aa), DIC1, SLA2: single exon, identical in the two independent CBS 7435 annotations (FR839631.1, 2011/2016; LT962479.2, 2022; the 2022 proteins CAH2450367.1/CAH2450368.1 are identical to the 2011 ones). I looked for no introns beyond that and did not re-examine RNA-seq.
- MATalpha2 (173 aa) has one 57 bp intron after an 11 bp first exon (GS115 NC_012966.1: exon 1 1648057-1648067 ATGATTTCTCA, intron 1648068-1648124 GTAAGT...CATTAG, exon 2 1648125-1648635). Read through, the ATG is blocked by a stop after four codons (MISQ*, more stops at codons 12 and 18); the spliced ORF is open, has a homeodomain (PF00046, 27.0 bits, GA 22.9) and hits MATalpha2 of S. cerevisiae (65.6 bits), Z. rouxii (63.4) and Candida (59.3). Same model in both CBS 7435 annotations. The curated record validates it: translation of the GS115 exons equals SCV12330.1 (100%/100%).
- The GS115 RefSeq models here are truncated and are NOT used: PAS_chr4_0878/XP_002494327.1 (161 aa, no intron) covers part of MATalpha2; PAS_chr4_1001/XP_002494328.1 (103 aa) is the shifted-frame remainder of MATalpha1.
- **GS115 MATalpha1 is frameshifted.** One extra A at NC_012966.1:1649032 (context GAAACT[A]CCTCG, about codon 47 of 190) relative to CBS 7435. The GS115 ORF from the ATG is 50 aa (26% of CCA40251.1). A 41 bp probe shows the GS115 allele in GCF_000027005.1 and one other of the 188 Komagataella assemblies in the campaign library (GCA_029632965.1, "K. pastoris"), and the CBS 7435 allele in GCA_900235035.2 and 7 others (the other ~178 assemblies are too fragmented around the site). Whether this is a true GS115 lesion or a shared artefact was not checked (no reads). Recorded as `present: false` (no coordinates, no protein), the existing remnant convention, so that a 50-aa fragment is not searched as MATalpha1. This leaves the GS115 record with one core gene (MATalpha2).

### Genetic code

Standard (table 1): NCBI Taxonomy `GeneticCode` 1 for taxids 460519 (K. phaffii), 981350 (CBS 7435) and 644223 (GS115); the detect reports also resolve `genetic_code: 1`. Every CDS translates without internal stop under table 1 (validation re-translates each gene and matches the deposited protein at 100%).

### Strain labels

- CBS 7435 = NRRL Y-11430: stated in PMID 36520709 (Offei et al. 2022, abstract: "CBS7435 (NRRL Y-11430)", all common expression hosts derive from this background) and PMID 31654411 (Brady et al. 2019, "Y-11430 (CBS7435)"). PMID 38443579 (Claes et al. 2024) calls NRRL Y-11430 the basic strain from which all commercial hosts derive and says comparative sequencing indicates it derives from a K. phaffii type strain deposited in 1954. PMID 21575661 (Küberl 2011) calls CBS 7435 the parental strain of commonly used hosts.
- GS115: PMID 31654411 compares Y-11430 (CBS7435), GS115 and X-33 and reports that X-33 shares nearly all mutations and a defective HIS4 with GS115, i.e. GS115 is a derivative variant of the Y-11430 lineage with a his4 defect. Cregg et al. 1985 (PMID 3915774, abstract) describe a his4 (histidinol dehydrogenase) auxotrophic host.
- NOT VERIFIED: that GS115 was made by chemical mutagenesis of Y-11430 (not stated in any text I could open; Cregg 1985 full text had no content in PMC here and the abstract names neither GS115 nor the mutagen), the Küberl 2011 and De Schutter 2009 (PMID 19465926) full texts (not read), and the NCBI record of GS115 (BioSample SAMEA2272385 says only STRAIN GS115). The NRRL Y-11430 label is put on the CBS 7435 record only (`culture_collection_ids`) and none on GS115.

## Family decision: new scope-only locus `MATkp`

Protein evidence (pyhmmer 0.12.3, phmmer of the four CBS 7435 proteins against all 209 curated proteins in `db/Ascomycota` and the Dipodascales proteins; Pfam PF04769, PF00046, PF00505, PF09011 from InterPro and the curator MAT1-1-1/MAT1-2-1 profiles; full output reproducible with `work/prot_evidence.py` on the HPCC, not committed):

| Kp protein | Profile hits | Best curated homologs (bits) |
|---|---|---|
| MATalpha1 (190 aa) | Pfam MATalpha_HMGbox 37.5 (GA 28.1, E 3e-13); curator MAT1-1-1 30.5 (E 7e-11) | MTLalpha1 (C. albicans 49, D. hansenii 56, others to 65); MATALPHA1 (S. cerevisiae 41) |
| MATalpha2 (173 aa) | Homeodomain 27.0 (GA 22.9, E 7e-10) | MATALPHA2 S. cerevisiae 65.6, Z. rouxii 63.4, MTLalpha2 C. albicans 58.3, MATA2 (HMR) 67.9 |
| MATa1 (225 aa) | Homeodomain 51.7 (E 1e-17) | MATsc MATA1 52.5, 47.7; MTLA1 37-41; Yarrowia mata1 22 |
| MATa2 (178 aa) | HMG_box PF00505 13.6 (GA 25.4, E 2e-5, sub-GA); no curator MAT1-2-1 hit | MTLA2 C. albicans 38.3, D. hansenii 29.0, Geotrichum MATA 21.3 (MATsc has no matching a2 reference) |

Reading: all four genes are homologous to genes in BOTH `MATsc` (Saccharomycetaceae) and `MTL` (Serinales); the data do not make Komagataella closer to one family than the other (alpha1 is closer to MTLalpha1, MATa1 to MATsc MATA1, MATa2 only to MTLA2). Considerations that decided for a new family:

1. The architecture is neither family's. `MATsc` is the three-cassette HML/MAT/HMR system with HO; `MTL` is anchored on PAP1/OBP1/PIK1 inside the idiomorph (flank window 3 kb). Komagataella has two cassettes with SLA2/DIC1 outside (flank window 20 kb, as `MATsc` sla2).
2. `MATsc` has per-gene rules built for HML/HMR (`MATA2`/`MATALPHA2` not informative); `MTL` has Serinales-specific optional PAP/OBP/PIK genes and a 3 kb window. Extending either to Pichiales would require unscoped exceptions.
3. Scope-only with a disjoint genus scope (460517) cannot be reported twice. Extending `MATsc` (scope 4893 = Saccharomycetaceae) to Pichiales would break the meaning of "sc"; extending `MTL` (Serinales) would break its documented scope.

Decision: new `MATkp`. Genus scope, not Pichiales/Pichiaceae, because only Komagataella is curated (Ogataea, also Pichiaceae, has the same inversion system per Hanson 2014 but no record). `fallback_searchable: false`. No `fallback_subphylum`: the key exists only on `dipodascales-curation` (not on main's schema); add `fallback_subphylum: 147537` after that branch merges.

What I did NOT decide (left for the curator, flagged in `order.yml` comments):
- `idiomorph_class` is `unassigned` for both idiomorphs. MATalpha -> MAT1-1 is strongly supported (PF04769 above GA); MATa -> MAT1-2 is only weakly supported (HMG_box sub-GA, the curator profile does not hit).
- `gene_class` not set on the four core genes (the `order.schema.yaml` evidence bar forbids an implementer from inferring it). `sla2` carries `sla2_homolog` by precedent (`MATsc`/`MATyl`), which the curator never separately ruled on for this family.
- `two_idiomorphs_report` not enabled (curator: "test first").
- Whether `system` should be `pseudohomothallic` (used: K. phaffii is homothallic and switches, PMID 25349420/29061733/31756769) rather than `homothallic`; the db convention for S. pombe h90 is `pseudohomothallic`.

## Before / after detection

`matpredict detect --taxid <own taxid>`, `MATPREDICT_HTML=0`, `PYTHONPATH` = this worktree's `src` (v0.6.2 code), 52 runs on the HPCC (array of at most 10 tasks, jobs named `kom-detect`). Before: `MATPREDICT_DB_ROOT` = `origin/main` `db`. After: this branch's `db`. Genome FASTAs from the BFD library. CBS 7435 and GS115 are run with the own record excluded (`--exclude-records`), plus the self-hit and both-excluded variants.

The curator's list names nine Pichiales genomes though it says eight; all nine are included.

| Genome | Species | Before (main db) | After (new db) |
|---|---|---|---|
| GCA_023013455.1 | K. phaffii | none | MATkp MATa JAJJTB010000045.1:117-2287 high [MATa1, MATa2, sla2] |
| GCA_023013825.1 | K. phaffii | none | MATkp MATa JAJJTU010000041.1:117-2287 high [MATa1, MATa2, sla2] |
| GCA_023013855.1 | K. phaffii | none | MATkp MATa JAJJTV010000046.1:117-2287 high [MATa1, MATa2, sla2] |
| GCA_023014715.1 | K. phaffii | none | MATkp MATa JAJJUH010000040.1:117-2287 high [MATa1, MATa2, sla2] |
| GCA_023013895.1 | K. phaffii | none | MATkp MATa JAJJTX010000047.1:110-2099 high [MATa1, MATa2] |
| GCA_023014845.1 | K. phaffii | none | MATkp MATa JAJJUO010000047.1:117-2106 high [MATa1, MATa2] |
| GCA_030557195.1 | K. kurtzmanii | none | MATkp MATa JAKTVM010000015.1:1-3210 high [MATa1, MATa2, sla2, dic1] |
| GCA_030567415.1 | K. pastoris | MTL A JAJLSE010000039.1:2-2045 high [MTLA1, MTLA2] | MATkp MATa JAJLSE010000039.1:2-2138 high [MATa1, MATa2, sla2] |
| GCA_000182785.1 | K. pastoris | MTL A CABH01000012.1:64-1954 high [MTLA1, MTLA2] | MATkp MATa CABH01000012.1:55-2047 high [MATa1, MATa2] |
| GCA_900235035.2 CBS 7435, own MATa excluded | K. phaffii | MATsc MATa 2941-7997 med; MATsc MATalpha 136012-146986 med; MTL alpha 136012-146986 med | MATkp MATalpha LT962479.2:136012-141828 high [MATalpha2, sla2, dic1] (the MATa cassette is NOT called: no other MATa reference) |
| CBS 7435, self-hit | | | MATkp MATa 1771-7997 high [MATa1, MATa2, sla2, dic1]; MATkp MATalpha 136012-141828 high [MATalpha2, sla2, dic1] |
| GCA_000223565.1 CBS 7435 (FR839631.1), MATa excluded | K. phaffii | MATsc MATa 2836-7892 med; MATsc MATalpha 135891-146865 med; MTL alpha med | MATkp MATalpha 135891-141707 high [MATalpha2, sla2, dic1] |
| GCF_000027005.1 GS115, own MATalpha excluded | K. phaffii | MATsc MATalpha 1641728-1652702 med; MTL alpha med | none (the alpha cassette is referenced only by the GS115 record) |
| GS115, self-hit | | | MATkp MATalpha 1646885-1652702 high [MATalpha2, sla2, dic1] |
| both records excluded, either strain | | | none |
| GCA_001708105.1 | K. pastoris | MATsc MATalpha med; MTL alpha med; MATsc MATa high; MTL A high | MATkp MATalpha CP014587.1:1438280-1444097 high [MATalpha2, sla2, dic1]; MATkp MATa 1569316-1575518 high [MATa1, MATa2, sla2, dic1] |
| GCA_023213065.1 | K. ulmi | MATsc MATalpha med; MTL alpha med; MATsc MATa high | MATkp MATalpha 1439297-1445114 high; MATkp MATa 1569573-1575776 high |
| GCA_900010155.1 | K. pastoris | none | MATkp MATa FBUB01000071.1:33-2022 high [MATa1, MATa2] |
| Ogataea polymorpha GCA_019207825.1, O. parapolymorpha GCA_019207595.1, O. glucozyma GCA_003706155.2, P. kudriavzevii GCA_003033855.1, C. boidinii GCA_002008015.1, Ambrosiozyma GCA_001599075.1, Brettanomyces GCA_003705215.3 | other Pichiaceae | MATsc/MTL/MATyl/PM/MAT calls | identical calls (routing unchanged: phylum_fallback; MATkp is scope-only) |

Reading:
- The nine contigs that were uncalled or called only by Candida-record fragments (MTL A, two genes) are now called MATa at high confidence with MATa1 and MATa2 (and SLA2/DIC1 when the contig reaches them). They were not used for any record.
- On full genomes the call replaces three weak medium calls (MATsc/MTL fragments) at each cassette with one high-confidence MATkp call per cassette. Routing for Komagataella is `lineage` after (was `phylum_fallback`): the other Ascomycota families are no longer searched there. This is the intended effect of scope-only families but also means Komagataella genomes lose any non-MATkp call.
- Leave-one-out: excluding the own record recovers the MATalpha cassette of CBS 7435 from the GS115 record (different strain, 100% identical sequence), but the MATa cassette is not recoverable (no other MATa reference), and the GS115 alpha cassette is not recoverable from the CBS 7435 MATa record (no shared core gene; the SLA2/DIC1 flanks alone do not make a call). This is a recall gap in the alpha direction: the db has exactly one alpha reference, from a strain whose MATalpha1 is frameshifted, so MATalpha1 is not searchable at all.

## Tests

`pytest tests/db tests/detect` on the HPCC (test environment of mat-prof, job kom-tests): 1116 passed, 4 skipped, 0 failed, 240 s. A first run with `-x` stopped at `test_shipped_roster_flags_the_four_lineage_families` because `MATkp` is a new scope-only family; the test was updated and the full suite rerun.

## Questions for the curator

See the final report; the same list: (1) new `MATkp` vs extension of `MATsc`/`MTL`, (2) add a CBS 7435 MATalpha record (intact MATalpha1/MATalpha2, FR839631.1:135890-141707), (3) GS115 MATalpha1 handling, (4) idiomorph_class and gene_class values, (5) genus vs Pichiales scope and Ogataea, (6) `system` value, (7) two_idiomorphs_report, (8) `fallback_subphylum` after the Dipodascales merge, (9) the nine-vs-eight genomes discrepancy.

## Not verified

- Intron/exon structure of all genes beyond annotation agreement and in-silico translation; no RNA-seq reads re-examined.
- GS115 MATalpha1 frameshift: real lesion or assembly artefact.
- Whether GS115 lacks MATa/locus 1 or the assembly does.
- GS115 derivation by chemical mutagenesis of NRRL Y-11430 (no open text).
- Küberl 2011, De Schutter 2009, Valli 2016 full texts (not read; only abstracts), Hanson 2014 and Feng 2020 abstracts only (Heistinger 2018 full text was read). That the Valli 2016 curation specifically supports these gene models was not checked.
- `homothallic` status of CBS 7435 and GS115 specifically (papers cover the species and CBS 2612).
- Claim that the idiomorph MATa maps to MAT1-2 / HMG box.

## Addendum, 2026-10-18 (curator follow-up)

- Decisions applied: `MATkp` accepted as its own scope-only family; a CBS 7435 MATalpha record added (`460519_cbs-7435_MATkp_MATalpha`, FR839631.1:135891-141707, anchor LT962479.2:136012-141828 at 100.000% over 5,818 bp, all four genes 100%/100%, `taxonomy_current: true` on a compute node); GS115 MATalpha1 kept `present: false` (extra A at NC_012966.1:1649032; the `status: pseudogene` label waits for PR #101); `idiomorph_class` kept `unassigned`; `gene_class: alpha_box` set on MATalpha1 only (Pfam PF04769 37.5 bits, GA 28.1; curator MAT1-1-1 30.5), MATa genes unset. `origin/main` merged into the branch (clean).
- Strain evidence added: the UniProt organism line of the CBS 7435 proteins reads "Komagataella phaffii (strain ATCC 76273 / CBS 7435 / CECT 11047 / NRRL Y-11430 / Wegner 21-1)" (rest.uniprot.org, 2026-10-17), an openable source for CBS 7435 = NRRL Y-11430.
- Detection rerun with the three records (10 runs, `work/det2` on the HPCC; `--taxid` own, `MATPREDICT_HTML=0`):

| Genome | Exclusion | After (3 records) |
|---|---|---|
| CBS 7435 GCA_900235035.2 | none | MATkp MATa 1771-7997 high [MATa1, MATa2, sla2, dic1]; MATkp MATalpha 136012-141828 high [MATalpha1, MATalpha2, sla2, dic1] |
| CBS 7435 | CBS MATa record | MATalpha cassette called (MATalpha1, MATalpha2, sla2, dic1); the MATa cassette is not (no other MATa reference) |
| CBS 7435 | CBS MATalpha record | MATa called; MATalpha called from the GS115 record (MATalpha2, sla2, dic1; no MATalpha1: the GS115 record has none) |
| CBS 7435 | both CBS records | MATalpha called from GS115 (MATalpha2, sla2, dic1); MATa not called |
| GS115 GCF_000027005.1 | none or GS115 record | MATkp MATalpha 1646885-1652702 high [MATalpha1, MATalpha2, sla2, dic1] (the recall gap of the first run is closed: the CBS 7435 MATalpha1 finds the frameshifted GS115 copy) |
| GCA_001708105.1 (K. pastoris) | none | MATalpha 1438280-1444097 and MATa 1569316-1575518, both high, all four genes |
| GCA_023213065.1 (K. ulmi) | none | MATalpha 1439297-1445114 and MATa 1569573-1575776, both high, all four genes |
| GCA_023013455.1, GCA_030557195.1 | none | MATa high (unchanged from the first run) |

- The MATa-gene analysis is in `2026-10-18_komagataella-mata-analysis.md`.
- Tests after the merge and the third record: `tests/db` + `tests/detect` 1116 passed, 4 skipped, 1 failed in the full run. The failure, `tests/detect/test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`, is intermittent and is not caused by this branch: it fails in about 1 of 3 runs on an unmodified `origin/main` checkout as well (it uses a temporary db root), and passes when run alone. Not investigated further.
