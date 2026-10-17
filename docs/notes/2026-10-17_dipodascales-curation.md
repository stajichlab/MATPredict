# Dipodascales curation, 2026-10-17

Branch `dipodascales-curation`. Follows `2026-10-17_dipodascales-investigation.md` (branch `dipodascales-refs`). Nothing is merged. Curator decisions of 2026-10-17: extend `MATyl` to the Yarrowia genus, curate the Yarrowia MATB idiomorph, add Geotrichum, address the Dipodascales wrong-contig calls.

## What changed

| File | Change |
|---|---|
| `db/Ascomycota/order.yml` | `MATyl`: `taxonomic_scope` [4952] -> [4951]; genes `matb2`, `matb1` (core_MAT, idiomorph B). New family `MATdc` (scope [43987], idiomorphs A/B, genes MATA, MATB, sla2, optional apn2). |
| `db/_schema/metadata.schema.yaml` | Evidence `tier` enum [1,2] -> [1,2,3]. Tier 3 = unpublished annotated deposit with no paper. Needed because the curator asked for tier 3 on the Geotrichum records; the schema allowed only 1 and 2. |
| `db/Ascomycota/Saccharomycetales/4952_clib122_MATyl_B/` | New MATB record (placed beside the MATA record so both MATyl records share a directory). |
| `db/Ascomycota/Dipodascales/1173061_clib918_MATdc_A/` | New. |
| `db/Ascomycota/Dipodascales/1173061_cbs615-84_MATdc_B/` | New, deposit-only. |

Workflow run on the HPCC: `propose`, `validate`, `accept`, `build-gff` for each record. All validations pass (sequence match 100% / 100% coverage for every gene, accession resolved, taxonomy current).

## Verified against NCBI (E-utilities, 2026-10-17)

- Taxids: Yarrowia 4951 (genus), Y. lipolytica 4952 (species), Geotrichum 43987 (genus, lineage ...Dipodascales; Dipodascaceae). The Geotrichum deposits carry taxon 1173061 (*Geotrichum candidum*, species), used as the record taxid.
- AJ617307.1: matb2 CDS 259..1362 (CAE84425.1), matb1 2362..3270 (CAE84426.1), strain field "W29". PMIDs 14745027, 18461384, 32105315, 10503550 resolve to the titles claimed.
- CLIB122 RefSeq (GCF_000002525.2, chromosome C): YALI0_C07436g apn2 (XP_501555.1) 979054..980802, YALI0_C07458g matb2 (XP_501556.1) 981060..982163, YALI0_C07480g matb1 (XP_501557.1) 983163..984071, YALI0_C07502g sla2 (XP_501558.1, minus) 984210..987371. The earlier agent's matb2 position (981060-982160) is right. Its apn2 start (979564) was a BLAST-hit edge; the RefSeq gene starts at 979054 (582 aa; the MATA-record apn2 CAE84424.1 of 412 aa is contained in it).
- AJ617307.1 is 3353/3353 nucleotides identical to NC_006069.1:980802-984154. CAE84425.1 and CAE84426.1 are identical to XP_501556.1 and XP_501557.1.
- `NC_006069.1` is flagged `suppressed` by nuccore esummary (stale: "removed as a result of standard genome annotation processing"), although the assembly report lists it as current chromosome C. The record's segment therefore uses `CR382129.1`, the INSDC copy of the same chromosome (sequence identical across the segment), so the accession check passes. `assembly_location.contig` keeps `NC_006069.1`, the name in the assembly FASTA.
- G. candidum CLIB 918 assembly GCA_001402995.1 carries the INRA annotation: SLA2 `CDO51887.1` (BN980_GECA02s02540g, CCBN010000002.1:579055..582249), "mating type protein MATA" `CDO51888.1` (BN980_GECA02s02545g, 584221..585066, identical protein to deposit CCO73812.1), then APC5. The earlier agent's MATA end (585063) omitted the stop codon. No APN2 beside SLA2, so apn2 is not in the record.
- HF558448.1 (CLIB 918, 846 bp) and HF558449.1 (CBS 615.84, 1020 bp): unpublished, "Specialization of the cheese isolates of the species Geotrichum", INRA/MICALIS, 07-NOV-2012.

## Decisions and assumptions

1. **MATyl scope** 4952 -> 4951 (curator). Strict superset.
2. **MATB record** uses CLIB122 (same-strain assembly, `assembly_accession_basis: strain_match_verified`, 100% over the whole AJ617307.1 deposit) as primary evidence. Strain 22301-5 (CP096873.1) is mentioned as corroboration only. Evidence tiers: existence 1 (Butler 2004, Rosas-Quijano 2008, Han 2020), boundaries 2 (RefSeq annotation), idiomorph assignment 1 (functional). `coordinate_provenance: curator_derived` (not `published_explicit`: the coordinates come from the RefSeq annotation, not a paper).
3. **Geotrichum family**: own family `MATdc`, scope Geotrichum genus only (not the order, to avoid reporting Yarrowia twice). A = HMG-like, B = alpha1-like, by position and weak homology only; Pfam was not run. `idiomorph_class` is `unassigned`. Whether A/B are homologous to Yarrowia A/B is an OPEN QUESTION and is not asserted anywhere.
4. **CLIB 918 record (A)** is assembly-anchored (same strain as the deposit, deposit CDS = assembly CDS). **CBS 615.84 record (B)** is deposit-only (no assembly of that strain in the library), so it has no flank and is excluded from the coordinate benchmark. LMA-244 was deliberately NOT curated, so it stays a held-out test genome.
5. `mating_type.system: heterothallic` on both Geotrichum records is unverified (the schema has no "unknown" value). A review snippet calls G. candidum self-fertile. Stated in each record.
6. `sla2` is required in `MATdc`; `apn2` is optional (not adjacent in CLIB 918).
7. Segment type `insdc_nucleotide` (contig or WGS accession) was used for all three, so `validate` independently re-translates every gene from the recorded coordinates.

## Detection before and after

`matpredict detect --taxid <genus taxid>` with `MATPREDICT_DB_ROOT` = this branch's `db`, v0.6.2 code from this branch (46 genomes, HPCC array job). "Before" is the v0.6.0 campaign call (`results/2026-10-03_ascomycota_v060/`). Distance d is from the earlier agent's SLA2 tblastn hit (`classified.tsv`); d=0 means the locus contains the SLA2 hit. The after-runs used the genus taxid, so Yarrowia and Geotrichum route directly; v0.6.0 used the genome's own taxid (lineage routing, same families).

| Run | Taxon | Before (v0.6.0) | After |
|---|---|---|---|
| CLIB122 | Yarrowia lipolytica | uncalled (lineage); - | direct: MATyl B NC_006069.1:979054-987371 high [matb2|matb1|sla2|apn2] d=0 |
| CLIB122_selfexcl | Yarrowia lipolytica | uncalled (lineage); - | direct: none |
| 22301-5 | Yarrowia lipolytica | uncalled (direct); - | direct: MATyl B CP096873.1:1023291-1031608 high [matb2|matb1|sla2|apn2] d=0 |
| W29 | Yarrowia lipolytica | MATyl (direct); MATyl A LJBI01000003.1:964153-970900 d=0 | direct: MATyl A LJBI01000003.1:963643-970900 high [mata1|mata2|sla2|apn2] d=0 |
| W29_selfexcl | Yarrowia lipolytica | MATyl (direct); MATyl A LJBI01000003.1:964153-970900 d=0 | direct: none |
| Yyakushimensis | Yarrowia yakushimensis | uncalled (phylum_fallback); - | direct: MATyl B ULGW01000001.1:4453598-4461501 high [matb2|matb1|sla2|apn2] d=0 |
| CLIB918_selfexcl | Geotrichum candidum | MAT (phylum_fallback); MAT MAT1-2 CCBN010000009.1:933886-950991 d=other contig | direct: none |
| CLIB918 | Geotrichum candidum | MAT (phylum_fallback); MAT MAT1-2 CCBN010000009.1:933886-950991 d=other contig | direct: MATdc A CCBN010000002.1:579055-585063 high [MATA|sla2] d=0 |
| LMA244 | Geotrichum candidum | uncalled (phylum_fallback); - | direct: MATdc B QQZM01000017.1:282101-286546 high [MATB|sla2] d=0 |
| Gghanense | Geotrichum ghanense | uncalled (phylum_fallback); - | direct: MATdc B JAKKCX010000006.1:1653133-1675406 high [MATA|MATB|sla2] d=0 |
| out_Magnusiomyces | Magnusiomyces paraingens | uncalled (phylum_fallback); - | phylum_fallback: MATdc B NW_022474068.1:4134838-4141060 high [MATB|sla2] d=0 |
| out_Wickerhamiella | Wickerhamiella versatilis | uncalled (phylum_fallback); - | phylum_fallback: none |
| out_Starmerella | Starmerella bombicola | uncalled (phylum_fallback); - | phylum_fallback: none |
Self-excluded runs are not independent tests of a different idiomorph: A and B are different genes, so with the only A (or B) record withheld nothing is left to find. The independent tests are 22301-5 (B, found from the CLIB122 record), Y. yakushimensis (B, 82%), LMA-244 (B, found from the CBS 615.84 deposit at ~97%) and the other wrong-contig genomes below.

G. ghanense: both MATA-like and MATB-like genes sit within 22 kb of SLA2 and are reported as one locus labelled B (genes MATA, MATB, sla2). Not modelled; possible homothallic-like arrangement.

### All 34 wrong-contig Dipodascales calls (v0.6.0), re-run

| Outcome after | n |
|---|---|
| Call within 30 kb of SLA2, no wrong-contig call | 20 (Geotrichum 16, Yarrowia 2, Daia 1, Zygoascus 1) |
| Call within 30 kb of SLA2, a wrong-contig call also remains | 3 (Blastobotrys, Tardiomyces, Zygoascus) |
| No call, wrong-contig call remains | 9 (Sugiyamaella 4, Wickerhamiella 2, Nadsonia 1, Tardiomyces 1, Yarrowia 1) |
| No call at all | 2 (Geotrichum; SLA2 on 1-2 kb contigs) |

Totals: 23 anchored, 12 with a wrong-contig call remaining (11 `phylum_fallback` genomes plus the fragmented Yarrowia GCA_030581655.1, which has MATyl A and B calls on tiny contigs away from its SLA2 contig). Geotrichum: 16 anchored + 2 uncalled, no wrong call remains: Geotrichum now routes directly to `MATdc`, so the Pezizomycotina `MAT` family is no longer searched.

Per-genome rows (genome, genus, species, v0.6.0 call, routing now, calls now):

```
GCA_037042115.1_ASM3704211v1	Wickerhamiella	P97-06A	MAT MAT1-1	phylum_fallback	MAT MAT1-1 low d=other contig; MATyl B low d=other contig
GCA_054122405.1_ASM5412240v1	Geotrichum	candidum	MAT MAT1-2	direct	none
GCA_023627875.1_ASM2362787v1	Geotrichum	candidum	MAT MAT1-2	direct	none
GCA_030563525.1_ASM3056352v1	Nadsonia	commutata	MATtub MAT1-2	phylum_fallback	MATtub MAT1-2 high d=other contig
GCA_054122455.1_ASM5412245v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc B high d=0
GCA_003705705.2_ASM370570v2	Blastobotrys	raffinosifermentans	MAT MAT1-2	phylum_fallback	MATdc B high d=0; MATyl B low d=0; MAT MAT1-2 high d=other contig
GCA_054122495.1_ASM5412249v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc B high d=0
GCA_054124535.1_ASM5412453v1	Geotrichum	candidum	MATsc MATalpha	direct	MATdc B high d=0
GCA_024734315.1_ASM2473431v1	Tardiomyces	blankii	MATtub MAT1-2	phylum_fallback	MATyl B low d=0; MATtub MAT1-2 medium d=other contig
GCA_025605225.1_LMA-1146_V2	Geotrichum	bryndzae	MAT MAT1-2	direct	MATdc B high d=0
GCA_050930895.1_ASM5093089v1	Sugiyamaella	novakii	PM P	phylum_fallback	PM P medium d=other contig
GCA_054125585.1_ASM5412558v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
GCA_001402995.1_New2.3_08062011	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
GCA_030565365.1_ASM3056536v1	Zygoascus	flipseniorum	MAT MAT1-2	phylum_fallback	MATdc A high d=0
GCA_030581655.1_ASM3058165v1	Yarrowia	lipolytica	MATyl A	direct	MATyl A high d=other contig; MATyl B high d=other contig
GCA_054125645.1_ASM5412564v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc B high d=0
GCA_013365075.1_LMA-77_final	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
GCA_025234805.1_LMA-563_1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
GCA_030579895.1_ASM3057989v1	Sugiyamaella	grinbergsii	MAT MAT1-2	phylum_fallback	MAT MAT1-2 high d=other contig
GCA_054124675.1_ASM5412467v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc B high d=0
GCA_056776975.1_ASM5677697v1	Wickerhamiella	versatilis	MAT MAT1-1	phylum_fallback	MAT MAT1-1 medium d=other contig
GCA_030564625.1_ASM3056462v1	Zygoascus	hellenicus	MAT MAT1-2	phylum_fallback	MATdc A medium d=0; MAT MAT1-2 medium d=other contig
GCA_030566045.1_ASM3056604v1	Geotrichum	geniculatum	MAT MAT1-2	direct	MATdc B high d=0
GCA_030569095.1_ASM3056909v1	Daia	tannicola	MAT MAT1-2	phylum_fallback	MATdc A medium d=0
GCA_030573285.1_ASM3057328v1	Tardiomyces	blankii	MATtub MAT1-2	phylum_fallback	MAT MAT1-2 high d=other contig; MATtub MAT1-2 high d=other contig
GCA_054123975.1_ASM5412397v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc B high d=0
GCA_054125385.1_ASM5412538v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
GCA_014596635.1_LMA-40_V2	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
GCA_030575175.1_ASM3057517v1	Sugiyamaella	novakii	PM P	phylum_fallback	PM P medium d=other contig
GCA_030583505.1_ASM3058350v1	Sugiyamaella	castrensis	MAT MAT1-2	phylum_fallback	MAT MAT1-2 medium d=other contig
GCA_900518985.1_YAAL0	Yarrowia	alimentaria	MAT MAT1-2	direct	MATyl B high d=0
GCA_030560665.1_ASM3056066v1	Yarrowia	alimentaria	MAT MAT1-2	direct	MATyl B high d=0
GCA_043388505.1_ASM4338850v1	Geotrichum	albidum	MAT MAT1-2	direct	MATdc B high d=0
GCA_054124005.1_ASM5412400v1	Geotrichum	candidum	MAT MAT1-2	direct	MATdc A high d=0
```

## Routing for Dipodascales genera outside any family

The curator ruling of 2026-09-26 (`not_searched` when nothing matches) does NOT apply here: `route()` falls to `phylum_fallback` whenever the phylum has any family, and Ascomycota has ten. Measured: *Magnusiomyces paraingens*, *Wickerhamiella versatilis* and *Starmerella bombicola* all report `routing_mode: phylum_fallback` (both before and after this branch), and are searched against every Ascomycota family. 1,472 of 19,415 genomes in the v0.6.0 campaign were routed this way (Pichiales 575, Dipodascales 346, Phaffomycetales 187, ...).

New side effect of this branch: `MATdc` and `MATyl` are now also searched on phylum fallback, so they fire on out-of-scope genera. M. paraingens gets `MATdc B` high at SLA2 (d=0), Zygoascus and Daia get `MATdc A` (d=0), which look like true alpha1/HMG genes beside SLA2, but Wickerhamiella P97-06A, Tardiomyces blankii and Blastobotrys also get `MATyl B` low-confidence single-gene calls (matb1 + sla2 + apn2), which are probably weak homology.

## Proposed rule (not implemented; no `src/` change)

The cause of the 32 of 34 Pezizomycotina-family calls (MAT 27, MATtub 3, PM 2) is that families scoped to Pezizomycotina (`MAT`, `MATtub`) and Taphrinomycotina (`PM`) are fallback-searchable on Saccharomycotina genomes. Minimal rule, two lines of data and no code (the mechanism exists: `fallback_searchable`, curator ruling 2026-10-01, used in Basidiomycota):

```yaml
# db/Ascomycota/order.yml
- locus_name: "MAT"      # + fallback_searchable: false
- locus_name: "MATtub"   # + fallback_searchable: false
- locus_name: "PM"       # + fallback_searchable: false
```

Cost: Orbiliomycetes (67 genomes, 4 called by `MAT` on fallback) and the unclassified "NA"-class genomes (15 fallback `MAT` calls, 9 `MATtub`) lose that search, because some are Pezizomycotina outside the leotiomyceta and Pezizomycetes scopes; they would need their own scope entry. A fuller alternative is a subphylum guard inside the fallback branch of `route()` (keep only families whose scope shares the query's subphylum). Optionally also `fallback_searchable: false` on `MATdc` and `MATyl` to stop them firing out of scope (see question 3 in the report).

## Tests

HPCC pixi `test` environment (added with `pixi run -e test`; the default environment has no pytest): `python -m pytest tests/db tests/detect -q` on this branch: 1073 passed, 0 failed, 7 warnings. No collection error there (pygenomeviz is installed in that environment).

# Round 2 (curator answers of 2026-10-17)

## 1. Routing: subphylum guard (implemented)

`route()` (src/MATPredict/detect/family_registry.py), phylum-fallback branch only. New optional locus field `fallback_subphylum` (NCBI taxid; schema `order.schema.yaml` updated). In a fallback the genome's subphylum is whichever *declared* `fallback_subphylum` taxid appears in the lineage that routing has already fetched; a family is then searched only if it declares no subphylum or the genome's one.

Why this design (smallest clean one): the field is declared per family, not derived from `taxonomic_scope`, because scope taxids are of mixed rank (leotiomyceta, a class, a genus, four species) and deriving a subphylum would need one lineage lookup per family per run. The genome side costs no extra lookup and has no hard-coded subphylum table (the candidate set is the declared values). A genome whose lineage has none of the declared subphyla (unclassified Ascomycota, or a failed lineage lookup, which is already recorded in `routing_error`) is not restricted, exactly as before; so unclassified Pezizomycotina keep the `MAT` fallback. An undeclared family is never restricted. `explicit_phylum` (`--phylum`) and `exhaustive` are untouched.

Declared in `db/Ascomycota/order.yml`: `MAT`, `MATtub` = Pezizomycotina 147538; `PM`, `mat1`, `mat2`, `mat3` = Taphrinomycotina 451866; `MATsc`, `MATyl`, `MATdc`, `MTL` = Saccharomycotina 147537 (all three taxids confirmed by NCBI Taxonomy).

Tests (`tests/detect/test_family_registry.py`): a Saccharomycotina genome with no family gets no `MAT`/`MATtub`/`PM`; an Orbiliomycetes-style Pezizomycotina genome still gets `MAT` and `MATtub`; Taphrinomycotina gets `PM`; no-subphylum and failed-lookup genomes are unrestricted; explicit `--phylum` unchanged; real `order.yml`: Magnusiomyces falls back to none of MAT/MATtub/PM/mat1-3/MATyl/MATdc, Yarrowia routes to MATyl only.

Side effect to know: `fallback_searchable` also governs `--phylum` runs (`route()` filters on it there too), so `MATyl`/`MATdc` are not searched by `--phylum Ascomycota`, only by taxid routing or `--exhaustive`.

## 2. MATyl and MATdc are scope-only

`fallback_searchable: false` on both (the existing curator-2026-10-01 mechanism). They now fire only for genomes in the Yarrowia / Geotrichum genus.

## 3. A/B homology: profile results

Profiles: Pfam HMG_box PF00505.26 (GA 25.4), HMG_box_2 PF09011.17, MATalpha_HMGbox PF04769.19 (GA 28.1; this is the alpha1/MAT alpha 1 HMG-box family; InterPro API), plus the curator's MAT1-1-1 (alpha box) and MAT1-2-1 (HMG box) profiles from `origin/ascomycota-core-hmms` (`MAT_core`, no GA cutoff). pyhmmer 0.12.3, E<=10 reporting; protein lengths: Geotrichum MATA 281, MATB 339; Yarrowia mata1 176, mata2 291, matb1 302, matb2 367.

| Protein | Profile | domain bits | domain i-E | span on protein (aa) |
|---|---|---|---|---|
| Geotrichum MATA | curator MAT1-2-1 (HMG) | 39.4 | 1.8e-13 | 70-160 |
| Geotrichum MATA | MATalpha_HMGbox | 13.0 | 1.5e-05 | 80-126 (below GA) |
| Geotrichum MATB | MATalpha_HMGbox | 34.9 | 3.0e-12 | 78-146 (above GA 28.1) |
| Geotrichum MATB | curator MAT1-1-1 (alpha) | 22.6 | 2.7e-08 | 76-142 |
| Yarrowia matb1 | MATalpha_HMGbox | 14.4 | 5.8e-06 | 99-160 (below GA) |
| Yarrowia mata1, mata2, matb2 | all five profiles | no hit at E<=10 | | |
| any | HMG_box PF00505, HMG_box_2 PF09011 | no hit at E<=10 | | |

Protein-to-protein (phmmer, E<=100): Geotrichum MATA <-> Yarrowia mata2 (42 bits, E 2e-14); Geotrichum MATB <-> Yarrowia matb1 (26 bits, E 2e-9); no other pair, and in particular no A-versus-B cross hit.

What the scores show: Geotrichum MATA is HMG-box-like (curator HMG profile, strong; Pfam PF00505 itself does not hit, so "HMG-box-like" rests on the curator profile) and Geotrichum MATB is alpha-box-like (Pfam above its GA plus the curator alpha profile). Yarrowia mata2 corresponds to Geotrichum MATA and matb1 to Geotrichum MATB by sequence similarity, so the labels A and B line up across the two genera (A = HMG-like, B = alpha-like), with the caveat that the Yarrowia proteins carry no profile hit of their own for A (mata2's HMG box from the literature is NOT confirmed by these profiles) and only a sub-GA alpha-box hit for matb1. mata1 and matb2 have no detectable counterpart in Geotrichum. Labels were not renamed; `idiomorph_class` stays `unassigned` for both families (a mapping A -> MAT1-2, B -> MAT1-1 would be supported for Geotrichum; for Yarrowia only by the phmmer link). Recorded in the `MATdc` and `MATyl` comments of `order.yml`.

## 5. homothallic_candidate for MATdc, and G. ghanense literature

Mechanism: the existing different-domains path of `classify_locus` (curator ruling 2026-09-25): both core genes full-length models (>=50% of the reference protein, no cross-match), different `gene_class`, not overlapping, within 25 kb. `MATdc` genes now carry `gene_class` (`MATA` HMG_box, `MATB` alpha_box), supported by section 3. No new code; the report's `two_idiomorphs` machinery is not enabled for MATdc (`two_idiomorphs_report` not set). Unit test `tests/detect/test_homothallic_matdc.py` uses the G. ghanense GCA_022023575.1 pattern (MATB 1,653,133-1,653,648; MATA 1,667,948-1,669,624; sla2 1,672,203-1,675,406). Whether the real detect run emits the label for G. ghanense is in the retest section.

Literature on Geotrichum ghanense and Galactomyces (PubMed, NCBI nuccore/protein searches 2026-10-17; abstracts read, full texts not):
- G. ghanense: only PMID 19671722 (Int J Syst Evol Microbiol 2010, description from Ghanaian cocoa fermentations; type strain G6 = CBS 11010 = CCUG 56722). The abstract calls it "asporogenous" and says ascospores were not produced. No mating-type or sexual-state data. No GenBank MAT deposit for it (nuccore/protein searches for "Geotrichum ghanense" with mating/MAT: 0 records). Whether GCA_022023575.1 is CBS 11010 is from the earlier agent's strain match and I did not re-verify it.
- Genus Galactomyces/Geotrichum: PMID 5463862 (Science 1970, "Sexual reproduction on Geotrichum candidum": mating types reported); PMID 22494250 (Phytopathology 2012: 1:1 mating-type segregation in G. citri-aurantii and a mixed outcrossing/asexual system proposed for G. citri-aurantii and G. geotrichum); PMID 22798647 (IJSEM 2012: mating studies used to synonymise several Geotrichum species with Galactomyces candidus / pseudocandidus); the two unpublished deposits HF558448/HF558449 (G. candidum MATA, MATB).
- Not found: any curated example of a homothallic or both-idiomorph Geotrichum, any MAT locus characterisation for G. ghanense or other Geotrichum species beyond the deposits. bioRxiv was not searched with a tool I could verify (no usable bioRxiv search in this session), so absence there is unverified. The ghanense double-idiomorph pattern therefore has no external curated confirmation; it may equally be a fragmented/duplicated locus or a second paralog.

## 6. Directories and lineages

NCBI Taxonomy (efetch, 2026-10-17): Geotrichum = ...Dipodascomycetes; Dipodascales; Dipodascaceae; Geotrichum. Yarrowia lipolytica = ...Dipodascomycetes; Dipodascales; Dipodascales incertae sedis; Yarrowia. Both are Dipodascales, so unifying under `db/Ascomycota/Dipodascales/` is taxonomically logical, and was done: `4952_w29_MATyl_A` and `4952_clib122_MATyl_B` moved from `Saccharomycetales/` (git mv). The W29 record's lineage said `c__Saccharomycetes;o__Saccharomycetales;f__Dipodascaceae`; it is now `c__Dipodascomycetes;o__Dipodascales;f__Dipodascales_incertae_sedis` (`lineage_resolved_date` 2026-10-17). My first B record had the wrong family (Dipodascaceae); fixed.
Dependencies checked: nothing in src or tests keys on the directory name for these records (only fixture records in two tests use `Saccharomycetales`). `holdout.records_to_withhold` reads ranks from `taxonomy.lineage` prefixes (`o__`, `f__` ...): at ORDER radius a Yarrowia record now also withholds the Geotrichum records and vice versa (correct, they share Dipodascales), and at FAMILY radius the Yarrowia records are in their own family (`Dipodascales_incertae_sedis`) separate from Geotrichum (Dipodascaceae). `tests/db` and `tests/detect` pass (see the end).

## 7. AJ617307.1 / W29 strain label

- GenBank source features: both AJ617306.1 (MATA) and AJ617307.1 (MATB) say `/strain="W29"`, both from the Wolfe lab deposit linked to Butler et al. 2004 (PMID 14745027). The Butler full text could not be retrieved (PMC returned only metadata), so the paper's own statement of the source strain is unverified.
- blastn of both deposits against the W29 assemblies in the BFD library: GCA_001307795.1 (LJBI01000003.1: 99.74% over 8,175 bp to AJ617306.1), GCF_001761485.1 (NC_090772.1: 99.69%), GCA_003054345.1 (CP028450.1: 99.7% in three HSPs, one 5 kb gap-split alignment). None of the three has any hit to AJ617307.1. Controls: E122 MATA GCA_059623445.1 carries MATA (99.74%); CLIB122 GCF_000002525.2 and 22301-5 GCA_023272835.1 carry AJ617307.1 at 100.000% over all 3,353 bp.
- So all three W29 assemblies carry MATA. The existing MATA record (W29, GCA_001307795.1) is consistent. The W29 label on AJ617307.1 is the odd one: that sequence is a MATB allele (identical to CLIB122's) and is not in W29. No allele label was changed. Both records' notes now state this; the W29 record's note also records the two further W29 assemblies as sequence matches.

## 4. Accepted without action

Evidence tier 3 enum extension and `mating_type.system: heterothallic` for the Geotrichum records: accepted by the curator (2026-10-17), no change.

## Full retest (own-taxid `detect`, 1,524 genomes)

Set: all 398 Dipodascales genomes, plus every genome the v0.6.0 campaign routed by `phylum_fallback` (1,472: Pichiales, Phaffomycetales, Ascoideales, Saccharomycodales, Lipomycetales, Saccharomycetales, Trigonopsidales, Sporopachydermiales, Alloascoideales, Orbiliales 67, 39 of unassigned class), 1,524 distinct in all (the 34 wrong-contig genomes are inside the Dipodascales set). Command: `matpredict detect --genome ... --taxid <the genome's own taxid>` with this branch's `db` and `src` (HPCC, SLURM array of 10 genomes per task, MATPREDICT_HTML=0). Taxids: column 6 of `assembly_summary_GCA.txt`/`_GCF.txt`; four accessions are absent from those files (GCF_000002525.2 used 284591, GCF_014490615.1 and the Y. lipolytica ones 4952, GCF_001983325.1 4909, GCF_000441935.1 430498 from NCBI Taxonomy).

Deviations from the requested job spec: partition `stajichlab` and 3 GB per task instead of `batch` and 16 GB: with `batch` every array sat at `AssocGrpMemLimit` for over an hour (the account's memory cap was used by other arrays). Measured peak RSS of `detect` on these genomes was under 0.5 GB (sacct, 46 test jobs). Concurrency was capped at 40. All 1,524 runs exited 0.

Self-hit handling: the curated records come from CLIB122 (`4952_clib122_MATyl_B`), CLIB 918 (`1173061_clib918_MATdc_A`) and, as region evidence only, 22301-5. CLIB122 and CLIB 918 were run with `--exclude-records` of their own record and give no call (expected: each is the only record of its idiomorph; these two rows are self-excluded, not tests). W29 (existing MATA record) was run with its record present: a self-hit. 22301-5 (`GCA_023272835.1`, MATB high) is not a self-hit by record but its MATB locus is 100% identical to the CLIB122 reference, so it is not an independent test either. Independent: Y. yakushimensis (MATB, 82% matb2) and G. candidum LMA-244 (MATB, the CBS 615.84 deposit at ~97%).

IMPORTANT caveat on "before": v0.6.0 vs this run differ in code as well as in the routing change (the branch contains main up to #87, v0.6.2), so not every difference below is caused by the guard. To separate the two, 15 genomes with changed non-Dipodascales calls were also run on main (624048d, no guard) with their own taxid: see "Attribution".

### Before (v0.6.0) vs after, by order


| group | n | called before | called after | not_searched after | newly called | calls lost | calls >30kb from SLA2 before | after |
|---|---|---|---|---|---|---|---|---|
| Pichiales | 575 | 384 | 381 | 0 | 10 | 13 | 0 | 0 |
| Dipodascales | 398 | 114 | 157 | 0 | 80 | 37 | 34 | 3 |
| Phaffomycetales | 187 | 172 | 170 | 0 | 0 | 2 | 0 | 0 |
| Ascoideales | 80 | 48 | 44 | 0 | 0 | 4 | 0 | 0 |
| Saccharomycodales | 71 | 16 | 16 | 0 | 0 | 0 | 0 | 0 |
| Orbiliales | 67 | 6 | 4 | 0 | 0 | 2 | 0 | 0 |
| Lipomycetales | 59 | 56 | 4 | 0 | 0 | 52 | 0 | 0 |
| (none) | 39 | 20 | 21 | 0 | 1 | 0 | 0 | 0 |
| Saccharomycetales | 24 | 20 | 19 | 0 | 0 | 1 | 0 | 0 |
| Trigonopsidales | 16 | 5 | 0 | 0 | 0 | 5 | 0 | 0 |
| Sporopachydermiales | 4 | 4 | 4 | 0 | 0 | 0 | 0 | 0 |
| Alloascoideales | 3 | 3 | 3 | 0 | 0 | 0 | 0 | 0 |
| Symbiotaphrinales | 1 | 1 | 1 | 0 | 0 | 0 | 0 | 0 |

"Called" = at least one call. "Calls >30 kb from SLA2" counts calls on a different contig from, or more than 30 kb from, the earlier agent's SLA2 tblastn hit (`results/2026-10-17_dipodascales/` method); it can only be computed for the Dipodascales genomes (389 of 398 have an SLA2 hit), so it is 0 by construction elsewhere. `not_searched` after: 0 in every order (no genome in this set reached `not_searched`; routing is `lineage` for Yarrowia/Geotrichum, `phylum_fallback` otherwise).

### Dipodascales by genus

| group | n | called before | called after | not_searched after | newly called | calls lost | calls >30kb from SLA2 before | after |
|---|---|---|---|---|---|---|---|---|
| Geotrichum | 91 | 21 | 69 | 0 | 51 | 3 | 18 | 0 |
| Yarrowia | 80 | 49 | 78 | 0 | 29 | 0 | 3 | 3 |
| Wickerhamiella | 58 | 2 | 0 | 0 | 0 | 2 | 2 | 0 |
| Blastobotrys | 36 | 14 | 6 | 0 | 0 | 8 | 1 | 0 |
| Starmerella | 35 | 1 | 0 | 0 | 0 | 1 | 0 | 0 |
| Sugiyamaella | 24 | 9 | 0 | 0 | 0 | 9 | 4 | 0 |
| Magnusiomyces | 15 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Entelexis | 13 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Zygoascus | 8 | 3 | 0 | 0 | 0 | 3 | 2 | 0 |
| Nadsonia | 5 | 1 | 0 | 0 | 0 | 1 | 1 | 0 |
| None | 5 | 3 | 1 | 0 | 0 | 2 | 0 | 0 |
| Groenewaldozyma | 4 | 4 | 3 | 0 | 0 | 1 | 0 | 0 |
| Tardiomyces | 4 | 3 | 0 | 0 | 0 | 3 | 2 | 0 |
| Spencermartinsiella | 4 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Diddensiella | 3 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Casaregolazyma | 3 | 1 | 0 | 0 | 0 | 1 | 0 | 0 |
| Limtongella | 2 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Dengshuqunia | 2 | 2 | 0 | 0 | 0 | 2 | 0 | 0 |
| Daia | 2 | 1 | 0 | 0 | 0 | 1 | 1 | 0 |
| Deakozyma | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Middelhovenomyces | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Grinbergsozyma | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Crinitomyces | 1 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Totals over the 1,524 genomes: called 849 before, 824 after; 91 newly called, 116 genomes lost their call(s).

### Newly called (91 genomes)

| Group | Family | Idiomorph | Confidence | Genomes |
|---|---|---|---|---|
| Geotrichum | MATdc | B | high | 38 |
| Geotrichum | MATdc | A | high | 13 |
| Yarrowia | MATyl | B | high | 27 |
| Yarrowia | MATyl | A | high | 2 |
| Pichiales | MTL | A (idiomorph_gene_only) | high | 9 |
| Pichiales | MATsc | MATalpha (idiomorph_gene_only) | medium | 1 |
| unclassified | MAT | MAT1-2 (partial_locus) | low | 1 |

Geotrichum: 21 -> 69 called of 91; G. ghanense GCA_022023575.1 now reports `MATdc B ... homothallic_candidate` (JAKKCX010000006.1:1653133-1675406, genes MATA, MATB, sla2, confidence high): the new rule fires on real data. 22 of 91 Geotrichum genomes are still uncalled (contig-edge SLA2 hits and fragmented assemblies, not examined genome by genome).

### Calls lost (116 genomes), and what they were

Attributable to the guard (the lost call is from a Pezizomycotina `MAT`/`MATtub`, Taphrinomycotina `PM`, or out-of-scope `MATyl` family on a Saccharomycotina genome):

| Order | Lost genomes | Lost calls of which anchored (core gene plus SLA2/APN2 in the same locus, or within 30 kb of SLA2 for Dipodascales) |
|---|---|---|
| Dipodascales | 37 | 16 wrong-contig (the goal) + 25 anchored calls in 20 genomes (Blastobotrys 7, Sugiyamaella 5, Dengshuqunia 2, Zygoascus, Starmerella, Groenewaldozyma, Casaregolazyma, 2 unclassified) |
| Lipomycetales | 52 | 96 calls with core + SLA2/APN2 in the locus (e.g. `MAT1-2-1|SLA2|APN2`, `MAT1-1-1|MAT1-2-1|APN2|SLA2`), 4 without a flank gene |
| Trigonopsidales | 5 | 6 calls with core + flank |
| Pichiales / Ascoideales / Phaffomycetales / Saccharomycetales / Orbiliales | 22 | 3 with core + flank, 20 without a flank gene |

The Dipodascales row reproduces the goal: wrong-contig genomes 34 -> 2 (both fragmented Yarrowia: GCA_030581655.1 and GCA_059050905.1 have MATyl calls on tiny contigs that are not the SLA2 contig). All 18 wrong-contig Geotrichum calls are gone.

THE COST (needs a curator decision): the guard also removes anchored calls that look real. 52 of 59 Lipomycetales genomes (Lipomyces) lose a `MAT`/`MATtub` call whose locus holds a core gene plus SLA2 and APN2, and 20 Dipodascales genomes (Blastobotrys, Sugiyamaella, ...) lose a call within 30 kb of SLA2 (the alpha1-like and HMG genes the investigation found at SLA2). These came from the Pezizomycotina family searching on yeasts via the fallback; whether those loci are real MAT genes is the question the investigation raised (SLA2-anchored alpha1-like hits, median 33% identity, Pfam not run on them). They are lost because no Saccharomycotina family covers these genera, not because they were shown wrong.

### Attribution (code drift versus the guard)

15 genomes whose non-Dipodascales call changed were re-run on main 624048d (no guard, same own taxid): 3 of the 10 new Pichiales `MTL`/`MATsc` calls also appear on main (so they come from code changes since v0.6.0, not the guard); 7 Pichiales `MTL A` `idiomorph_gene_only` calls at contig starts (positions 119-2,107) appear ONLY on the branch (a side effect of the narrowed family set that I did not explain; these are low-evidence, no flank gene, and should be treated as suspect). On main the Orbiliales MATsc calls (2) and one Pichiales MATsc `MATa` partial call are present, and on the branch absent (the Orbiliales ones are the guard working as intended; Orbiliales still get `MAT`/`MATtub`). `GCF_038497685.1` (Lipomyces suomiensis) has `MAT`/`MATtub`/`MATsc MATalpha` calls on main, none on the branch (the MATsc call is lost with the others).

Routing check from `families_attempted`: Saccharomycotina genomes (Pichiales, Phaffomycetales, Ascoideales, Saccharomycodales, Lipomycetales, Saccharomycetales, Trigonopsidales, Sporopachydermiales, Alloascoideales, Dipodascales other) searched only `MATsc` and `MTL`; Orbiliales and Symbiotaphrinales searched `MAT` and `MATtub`; Yarrowia only `MATyl`, Geotrichum only `MATdc`; 39 genomes of unassigned class (no declared subphylum in the lineage) keep the unrestricted fallback, as designed (27 of them search all eight non-scope-only families). Dipodascales genera outside Yarrowia/Geotrichum (Magnusiomyces, Wickerhamiella, Starmerella...) are still `phylum_fallback`, never `not_searched`: they are searched against `MATsc` and `MTL` and nothing else.

## Tests (final)

HPCC pixi `test` env, `python -m pytest tests/db tests/detect -q`: 1,082 passed, 1 failed. The failure, `tests/detect/test_pipeline.py::test_fast_path_missing_gene_is_polished_in_the_existing_clusters_window`, asserts the call order `["exonerate", "miniprot"]` of two polishers run concurrently; it passed in the first full run of this branch (1,073 passed) and passes 5/5 locally, so it is an order-dependent flake unrelated to these changes. One failure that was mine (`test_fallback_searchable.py` pinned the flagged families to the four Basidiomycota ones) was fixed to include `MATyl` and `MATdc`. New tests: 7 in `test_family_registry.py` (guard), 3 in `test_homothallic_matdc.py`.

# Round 3: Pichiales `MTL A` calls and the two Orbiliales losses (analysis only)

No `src/` or `db/` change. Jobs `dipo-pich-*` (partition stajichlab: the `preempt` account was refused for my user).

## 1. The nine Pichiales `MTL A` idiomorph_gene_only calls

All nine are *Komagataella* (six K. phaffii, one K. kurtzmanii, two K. pastoris). Seven appear only on this branch; two (GCA_000182785.1, GCA_030567415.1, both K. pastoris) also appear on main 624048d. (The tenth new Pichiales call, Brettanomyces bruxellensis `MATsc` MATalpha, is not an MTL call and was not analysed.)

| Genome | Taxon | Contig (length) | Hit span | Dist to contig start / end | Contigs / N50 / size | Identity (coverage) of K. phaffii deposit proteins on the contig | Placement in K. phaffii CBS 7435 |
|---|---|---|---|---|---|---|---|
| GCA_000182785.1 | K. pastoris | CABH01000012.1 (2082 bp) | 64-1954 | 63/128 | 125 / 169 kb / 9.3 Mb | MATa1 92.889% cov 100; MATa2 89.326% cov 100; SLA2-2 cov -%; DIC1 cov -% | CBS 7435 LT962479.2 4892-2845 85.106% |
| GCA_023013455.1 | K. phaffii | JAJJTB010000045.1 (2287 bp) | 210-2097 | 209/190 | 48 / 559 kb / 9.4 Mb | MATa1 100.000% cov 100; MATa2 100.000% cov 100; SLA2-2 cov 5%; DIC1 cov -% | CBS 7435 LT962479.2 2732-5018 99.213% |
| GCA_023013825.1 | K. phaffii | JAJJTU010000041.1 (2287 bp) | 210-2097 | 209/190 | 44 / 909 kb / 9.3 Mb | MATa1 100.000% cov 100; MATa2 100.000% cov 100; SLA2-2 cov 5%; DIC1 cov -% | CBS 7435 LT962479.2 2732-5018 99.213% |
| GCA_023013855.1 | K. phaffii | JAJJTV010000046.1 (2287 bp) | 210-2097 | 209/190 | 47 / 525 kb / 9.3 Mb | MATa1 100.000% cov 100; MATa2 100.000% cov 100; SLA2-2 cov 5%; DIC1 cov -% | CBS 7435 LT962479.2 2732-5018 99.213% |
| GCA_023013895.1 | K. phaffii | JAJJTX010000047.1 (2215 bp) | 119-2006 | 118/209 | 48 / 602 kb / 9.3 Mb | MATa1 100.000% cov 100; MATa2 100.000% cov 100; SLA2-2 cov 3%; DIC1 cov -% | CBS 7435 LT962479.2 4946-2732 100.000% |
| GCA_023014715.1 | K. phaffii | JAJJUH010000040.1 (2287 bp) | 210-2097 | 209/190 | 42 / 525 kb / 9.3 Mb | MATa1 100.000% cov 100; MATa2 100.000% cov 100; SLA2-2 cov 5%; DIC1 cov -% | CBS 7435 LT962479.2 2732-5018 99.213% |
| GCA_023014845.1 | K. phaffii | JAJJUO010000047.1 (2215 bp) | 210-2097 | 209/118 | 47 / 701 kb / 9.3 Mb | MATa1 100.000% cov 100; MATa2 100.000% cov 100; SLA2-2 cov 3%; DIC1 cov -% | CBS 7435 LT962479.2 2732-4946 100.000% |
| GCA_030557195.1 | K. kurtzmanii | JAKTVM010000015.1 (111107 bp) | 155-2107 | 154/109000 | 36 / 591 kb / 9.2 Mb | MATa1 99.556% cov 100; MATa2 97.753% cov 100; SLA2-2 cov 4%; DIC1 cov 100% | CBS 7435 LT962479.2 140791-249721 98.016% |
| GCA_030567415.1 | K. pastoris | JAJLSE010000039.1 (2262 bp) | 155-2045 | 154/217 | 40 / 1004 kb / 9.4 Mb | MATa1 92.889% cov 100; MATa2 89.326% cov 100; SLA2-2 cov 4%; DIC1 cov -% | CBS 7435 LT962479.2 4983-2845 85.549% |
Method: blastn of each hit contig against CBS 7435 (GCA_900235035.2, LT962479.2) and GS115 (GCF_000027005.1, NC_012966.1); tblastn of the earlier agent's `Kphaffii_FR839631` (DIC1-1 CCA40181.1, MATa2 SCV12319.1, MATa1 CCA40182.1, SLA2-2 CCA40183.1) and `Ohpoly_AJ617305` proteins onto each contig and onto the reference regions. Note these deposit proteins are not curated records in `db/`.

What the detect call matched: the `MTL` family's references `498019_b11205_MTL_A` MTLA1 (C. auris, 39-42% identity, 190-250 bp) and `5476_sc5314_MTL_A` MTLA2 (C. albicans, 33-34%, 320 bp), both `polished_single` fragments of 191-321 bp, both genes within 2 kb on the same 2.1-2.3 kb contig; no SLA2/APN2 or OBP1/PIK1/PAP1 neighbour was found by detect.

What the contigs are (verified by blast): 8 of the 9 hit contigs are 2,082-2,287 bp and 99-100% identical (K. phaffii) or 85% (K. pastoris) over their whole length to CBS 7435 LT962479.2:2,732-5,018, the region between DIC1 and SLA2-2 at the start of that chromosome. There, in CBS 7435 by tblastn: DIC1-1 1,771-2,685, MATa2 3,381-2,848 (minus), MATa1 4,163-4,837, SLA2-2 7,997-4,863 (minus): the K. phaffii MATa locus, DIC1 - MATa2 - MATa1 - SLA2. The contigs carry MATa1 and MATa2 FULL LENGTH (100% coverage, 100% identity to the deposit in K. phaffii; 89-93% in K. pastoris; 98-100% in the K. kurtzmanii contig), and only the first 3-5% of SLA2-2 at the contig edge; DIC1 is outside the contig in 8 of 9. The kurtzmanii contig (111 kb) is the complete case: SLA2-2 (partial, contig start) - MATa1 - MATa2 - DIC1 (100% coverage), then 109 kb of the chromosome. The same chromosome also has the other cassette (CBS 7435 LT962479.2:136,012-141,828: SLA2-2, matalpha1/matalpha2 by Ohpoly proteins at 44-47%, DIC1) and GS115 NC_012966.1 carries the alpha configuration at the DIC1 - matalpha2 - matalpha1 - SLA2-2 site (1,646,885-1,652,702), so K. phaffii has both configurations in its genome, consistent with a switching or two-cassette arrangement. These contigs are the MATa state; I did not check the biology further.

Verdicts per genome: all nine are the REAL K. phaffii-group MATa locus (MATa1 and MATa2 genes), not paralogs: the MATa1/MATa2 sequences match the deposit set at 89-100% and sit in the DIC1/SLA2 synteny of CBS 7435. "Contig-edge fragment" applies only to the flanks: the contig is a short MAT-only contig produced by assembly breaks at the repeats around the locus (distance to contig end 118-217 bp = SLA2-2 truncated), so the locus is genuine but the call is flank-poor, which is why it is `idiomorph_gene_only` and not `mat_locus`. The MTL family call is right in substance but is made by weak cross-lineage hits (33-42% to Candida references) because no Pichiaceae/Komagataella MAT record exists.

### Why seven are branch-only

Not the guard's family narrowing itself and not a family-specific suppression: a polish-cap artefact. Same genome, same taxid, current main versus this branch, `--evidence-diagnostics`:
- On main the same 2-gene MTL cluster (MTLA1+MTLA2) is formed and admitted, then reported as suppressed `modelled_gene_bar` because 0 of its genes were polished (`polish_capped: true` in the diagnostics row for this cluster in 7/7 genomes; admitted MTL clusters in the genome 18-24 on main vs 15-18 on the branch, of which 12-18 vs 9-12 capped).
- On the branch, with only MATsc and MTL searched, the cluster is within the per-family cap (`max_polished_clusters_per_family`, 6) and is polished (2 `polished_single` models), so it is reported.
- Mechanism: clusters are formed jointly from hits of all searched families, so the number and ranking of admitted `MTL` clusters (and so who falls under the cap of 6) depends on which other families are searched. Nine families on main give 5 competing clusters on this contig (MAT 4 genes, MATsc, MATyl, MTL, MATtub), which changes the ranking. The two genomes that are also called on main are the ones where the cluster stayed within the cap.
- So the guard changed the result only indirectly (fewer families, different cluster set); the call itself is the correct detection and the main-branch silence is a cap-induced false negative. The v0.6.0 campaign (9 families) was subject to the same artefact, which is why the Komagataella MATa locus was "uncalled" in 100+ genomes. I did not run an uncapped main to confirm (a run with a larger cap would prove it directly); evidence is the diagnostics flag.

Separation from code drift: the 15-genome main re-run (previous section) shows 3 genomes (GCA_000182785.1, GCA_001719535.1, GCA_030567415.1) called on main as well as on the branch (drift since v0.6.0: v0.6.0 left them uncalled), and the other 7 only on the branch (cap artefact via family set).

### Proposed rule (not implemented)

Order the per-family polish cap so a cluster with two distinct core genes of one idiomorph (or any anchored cluster) is never capped out, or rank by independent evidence (gene count, best identity) before cluster position; then the result no longer depends on which families are searched. Separately, a Komagataella (Pichiaceae) MATa/MATalpha record would turn these from weak MTL cross-lineage calls into proper family calls.

## 2. The two Orbiliales losses

- GCA_047671915.1 (Hyalorbilia sp. HsImV27) and GCA_047716685.1 (Hyalorbilia oviparasitica ASM4771668v1), both class Orbiliomycetes, both `phylum_fallback` in v0.6.0.
- v0.6.0 call: `Ascomycota:MATsc`, idiomorph `undetermined`, `idiomorph_gene_only`, genes MATALPHA2|MATA2, a single 321-bp hit (JBAMKL010000648.1:5,396-5,716 and JBAMJO010000854.1:10,847-11,167).
- Now: routing `phylum_fallback`, families searched `MAT` and `MATtub` only; no call. The same genomes on main (9 families) still give the MATsc call.
- Cause: the guard (`MATsc` is declared Saccharomycotina, so it is no longer searched on an Orbiliomycetes genome). Not the scope-only change (MATsc is not scope-only) and not code drift (main reproduces the call). A Pezizomycotina class still gets the `MAT`/`MATtub` fallback: 4 other Orbiliales genomes (Dactylellina cionopaga, Arthrobotrys iridis, two Hyalorbilia) keep their `MAT` calls, including a second H. oviparasitica genome (GCA_047716655.1, MAT1-2-1|APN2|COX13).
- Verdict: the lost calls were a single 321-bp cross-lineage HMG-like fragment found by a Saccharomycetaceae family in a filamentous fungus; no flank gene; I read them as spurious, not real MAT loci (unverified: the fragments were not inspected further). Neither of the two genomes gets a `MAT` call now, as before.

## 3. Questions for the curator

1. Accept the nine Komagataella calls as real (MATa1/MATa2) and curate a K. phaffii MATa (and MATalpha) record so they route to a proper Pichiaceae/Komagataella family instead of weak MTL fallback hits? CBS 7435 LT962479.2:2,732-5,018 and the DIC1-MATa2-MATa1-SLA2-2 layout are the obvious source (the deposit `FR839631` is not yet in `db/`).
2. Fix the polish cap so that the result does not depend on which other families are searched (rank by evidence or exempt two-core-gene clusters)? This affects every family and needs the single-gene rule evaluation the curator already wants before the guard merges.
3. Should a call made only by weak cross-lineage hits (33-42%) on a contig-edge fragment be `high` confidence? These nine are `high` with `idiomorph_gene_only`; consider a lower tier until a same-genus record exists.
4. The kurtzmanii contig carries the whole MATa locus with SLA2 and DIC1 flanks yet is `idiomorph_gene_only` with no flank gene because the flank references (SLA2/DIC1 of the Serinales MTL records) do not match: a Komagataella record would fix that too.
5. Not verified: that the CBS 7435 and GS115 regions are the two mating-type configurations (layout from tblastn only; annotation and literature not checked); that K. pastoris at 85-93% identity is the same locus rather than a diverged copy (synteny supports it); and the uncapped-main confirmation of the cap explanation.
