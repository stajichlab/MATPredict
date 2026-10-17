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
