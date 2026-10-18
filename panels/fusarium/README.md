# Fusarium MAT idiomorph panels for `matpredict reads-type`

Fourteen new panels (`<panel>_MAT1-1.fasta`, `<panel>_MAT1-2.fasta`, one record each) plus the existing
`../fusarium_oxysporum_fola/` panel (not copied here). Built 2026-10-17 from NCBI assemblies; validated on
56 read tests (50 runs) streamed from ENA. Companion files: `validation.tsv` (every read test),
`species_panel_map.tsv` (species table below, with notes), `insilico_assembly_typing.tsv` (panel x assembled locus typing).

```
matpredict reads-type --idiomorph MAT1-1=panels/fusarium/pseudograminearum_MAT1-1.fasta \
                      --idiomorph MAT1-2=panels/fusarium/pseudograminearum_MAT1-2.fasta --reads R1.fq.gz R2.fq.gz --out out.tsv
```

## Main findings

1. **One panel does not cover Fusarium; roughly one panel per species or tight species group is needed.** At k=31 a
   panel types a strain only if >=~50% of its idiomorph k-mers occur in the strain. MAT1-1-1 protein identity to the
   *F. fujikuroi* MO44 query (db 5127_mo44) is 100% in *F. verticillioides*, 84% in *F. oxysporum*, 52-54% in the
   sambucinum clade and *F. graminearum*, 36-38% in the solani complex (tblastn of the five repo-db proteins against 366
   assemblies, `insilico_assembly_typing.tsv`). Fraction of one locus's 31-mers present in another species' locus
   (5 kb window, best pair): FGSC species 0.64-0.73 among themselves; *F. verticillioides* vs *F. musae* MAT1-1 0.84; *F. circinatum* vs
   *F. subglutinans* MAT1-1 0.62; *F. oxysporum* vs *F. odoratissimum* 0.83-1.0; *F. verticillioides* vs *F. andiyazi* MAT1-1 0.47;
   *F. fujikuroi* vs *F. proliferatum* MAT1-2 0.44; *F. pseudograminearum* vs *F. culmorum* MAT1-2 0.43; *F. sambucinum* vs *F. venenatum* MAT1-1
   0.33; solani complex: 32 MAT1-1 loci fall in 11 clusters and 18 MAT1-2 loci in 7 clusters at >=0.5.
   MAT1-1 and MAT1-2 diverge at different rates: the *F. culmorum* MAT1-1 panel record also types *F. pseudograminearum* and *F. cerealis*
   MAT1-1 (breadth 0.56-0.74 on reads) but their MAT1-2 is not typed (breadth 0.36 -> no call).
2. **Where each panel stops** (reads, `validation.tsv`): the Fg panel gives `low_depth` on *F. pseudograminearum* reads (shared depth 0.01); the
   Fps panel gives `low_depth` on *F. culmorum* MAT1-2 (0.36 breadth); the sambucinum panel on *F. venenatum* (0.30); the solani panel on a *F. falciforme* MAT1-2 allele
   and on solani IISc-1 (different MAT1-1 lineage). Worked across species: Fv panel -> *F. musae* (MAT1-1, 0.85), proliferatum panel -> *F. annulatum*
   (MAT1-2, 0.67), circinatum panel -> *F. temperatum* (MAT1-2, 0.59), fola panel -> *F. odoratissimum*, *purpurascens*, *tardichlamydosporum*, *veterinarium*.
3. ***F. graminearum* and the FGSC are homothallic and a correct call is `both`.** All assembled FGSC loci (graminearum PH-1, asiaticum x3,
   boothii x8, gerlachii x2, meridionale x2, aethiopicum x2, louisianense x2) carry MAT1-1-1/-2/-3 and MAT1-2-1 in one locus; I found no heterothallic FGSC member
   among the species in selected.tsv. The `graminearum` panel therefore cuts the two idiomorph records from PH-1 side by side (MAT1-2 record 265500-268500, MAT1-1 record 267990-274400,
   511 bp overlap = the shared control). Reads: PH-1 (x2), CML3066, NRRL 29169, *F. boothii* ET-2022-34, *F. aethiopicum* ET-2022-49 all `both`.
   Population-level departures from `both` (a lost half of the locus) were not searched for.
4. **The brief's Fol4287 = MAT1-2 is wrong for the reference genome and the reads:** RefSeq GCF_000149955.1 carries MAT1-1-1/-2/-3 at NC_030993.1:3023257-3028476 (99.6% identical to the Fola MAT1-1 record)
   and reads SRR3139042 type as MAT1-1 with the fola panel. Other brief examples: *F. fujikuroi* IMI58289 is MAT1-2 (HF679028.1, 96% to AF100926 MAT1-2-1) and *F. verticillioides* 7600 is MAT1-1
   (CP114046.1/CM000583.1); no WGS reads for either were available in ENA under those strain names, so they serve as panel sources/in silico references only.

## Panel design

Each record = the idiomorph "core" (largest gap of the colinear BLASTN chain between the MAT1-1 and MAT1-2 source loci) plus up to 600 bp of flank on each side, but a flank
is kept only if the two source loci align there at >=97% identity, so that flank k-mers are identical in both records and become the shared single-copy control
(`Panel`: 319-847 shared k-mers per panel). Source loci come from different strains or species lineages (graminearum: one locus), flank identity was 98.1-99.8% except sacchari (left flank 92%, dropped), solani (left 93%, dropped) and vanettenii (right 85%, dropped). Coordinates are 1-based inclusive on the forward strand of the named
sequence. Records contain no non-ACGT. The `reads-type` min-unique-run filter (100) removes the few SNP-created flank k-mers. Every record was re-checked with tblastn against the repo-db MAT proteins:
MAT1-1 records contain MAT1-1-1/-2/-3, MAT1-2 records MAT1-2-1 (sub-35%-identity HMG-box-like hits to the other idiomorph's queries are noise; graminearum MAT1-1 record/pseudograminearum/culmorum
MAT1-1 records show a 81-aa 32%-identity MAT1-2-1-like hit, as in the PH-1 locus itself).

| panel | MAT1-1 record (assembly, strain) | bp | MAT1-2 record (assembly, strain) | bp | unique k-mers 11 / 12 | shared |
|---|---|--:|---|--:|--:|--:|
| culmorum (F. culmorum) | `JBIQGB010000002.1:1585671-1591557` (GCA_052570695.1 336336) | 5887 | `JBIQGA010000006.1:1586012-1590766` (GCA_052570705.1 309752) | 4755 | 4915 / 3782 | 648 |
| falciforme (F. falciforme) | `JBMISE010000001.1:333765-343604` (GCA_049863005.1 UA-600A) | 9840 | `JAOQAU010000021.1:605589-613426` (GCA_027945495.1 A01-2) | 7838 | 8681 / 6686 | 585 |
| fujikuroi (F. fujikuroi) | `CP023106.1:1145947-1151465` (GCA_009663115.1 I1.3) | 5519 | `HF679028.1:1115055-1119902` (GCA_900079805.1 IMI 58289) | 4848 | 4582 / 3903 | 655 |
| graminearum (F. graminearum) | `NC_026475.1:267990-274400` (GCF_000240135.3 PH-1) | 6411 | `NC_026475.1:265500-268500` (GCF_000240135.3 PH-1) | 3001 | 5877 / 2489 | 481 |
| poae (F. poae) | `CM197424.1:3006269-3012462` (GCA_061603655.1 DAOMC 252228) | 6194 | `CM033786.1:3467058-3472362` (GCF_019609905.1 DAOMC 252244) | 5305 | 5181 / 4292 | 791 |
| proliferatum (F. proliferatum) | `JBVTKN010000006.1:1024266-1029801` (GCA_059629225.1 66420) | 5536 | `JBVTKM010000006.1:1018806-1023654` (GCA_059629185.1 66442) | 4849 | 4489 / 3802 | 641 |
| pseudograminearum (F. pseudograminearum) | `CM003199.1:3030886-3036777` (GCF_000303195.2 CS3096) | 5892 | `CM004557.1:3006703-3011430` (GCA_001703955.2 RBG5266) | 4728 | 4940 / 3777 | 743 |
| sacchari (F. sacchari) | `JALNLS010000010.1:1162081-1167008` (GCA_048414715.1 CNO1) | 4928 | `JADANP010000006.1:287877-292090` (GCA_017165645.1 FS66) | 4214 | 4455 / 3741 | 381 |
| sambucinum (F. sambucinum) | `CP167250.1:6122684-6128889` (GCA_041870755.1 hop_citra) | 6206 | `JBJXRM010000002.1:6013001-6018204` (GCA_047301595.1 GUCCf09) | 5204 | 5234 / 4234 | 847 |
| solani (F. solani) | `JAMBZA010000017.1:1517673-1526510` (GCA_023522795.1 SB1) | 8838 | `JBHMJZ010000983.1:42964-49890` (GCA_045529305.1 CBS 140079) | 6927 | 8381 / 6469 | 319 |
| vanettenii (F. vanettenii) | `GG698897.1:2934967-2940328` (GCA_000151355.1 77-13-4) | 5362 | `CM049509.1:3403522-3407895` (GCA_025433535.2 Fs6.1) | 4374 | 4768 / 3781 | 371 |
| verticillioides (F. verticillioides) | `CP114046.1:999334-1005139` (GCA_027571605.1 7600) | 5806 | `CM132614.1:1076369-1081392` (GCA_053843395.1 FvLNF15-11) | 5024 | 4831 / 4049 | 646 |
| ../fusarium_oxysporum_fola (existing) | see its README | 5220 | see its README | 5103 | 4717 / 4619 | 362 |

Missing pieces: no panel has a sequence for lineages I could not assemble (see species table). For a species whose assemblies show only one idiomorph, a
reads-type panel cannot be built (it needs both idiomorphs plus shared flank), e.g. *F. andiyazi*, *F. venenatum*, *F. avenaceum*, *F. keratoplasticum*, *F. euwallaceae*, *F. nygamai*
(only MAT1-1 found among 1-16 assemblies each). *F. sacchari* carries two divergent MAT1-1 alleles (CNO1 in the panel, PB09 held out); the sacchari panel types PB09 only at breadth 0.33.

## Species in selected.tsv -> panel

Status key: `reads-validated` = at least one run of this species typed correctly with this panel against an independent truth (same-strain assembly) and all such tests passed;
`assembly-supported` = every assembled locus of the species (both idiomorphs, or `both` for homothallic) is typed correctly in silico (150-bp windows, 10-bp step, same `reads-type` call logic) but no reads were tested;
`partial` = some assembled loci typed, or only one idiomorph tested (the other idiomorph unverified); `reads x/y pass` prefixes a species with a failed test;
`inferred (untested)` = no assembly of the species; panel suggested only from complex membership or a doubtful species label (cf./sp.); `no panel` = no panel types it (no assembly,
single idiomorph assembled, or a divergent lineage). Run counts are rows of selected.tsv (960). Totals: reads-validated 414 runs (19 species), assembly-supported 118 (7), partial or reads-with-failures 104 (13), inferred 49 (9), no panel 275 (79). A panel being listed does not guarantee a call for every strain of the species: allele lineages
inside a species can fall outside it (F. solani, F. sacchari, F. proliferatum 1 of 18 assembled loci, F. circinatum 2 of 17).

| species | group | runs | panel | status |
|---|---|--:|---|---|
| F. fujikuroi | fujikuroi_complex | 50 | fujikuroi | assembly-supported |
| F. proliferatum | fujikuroi_complex | 50 | proliferatum | reads 3/4 pass; partial |
| F. graminearum | graminearum_complex | 50 | graminearum | reads-validated |
| F. culmorum | other | 50 | culmorum | reads-validated |
| F. poae | other | 50 | poae | reads-validated |
| F. pseudograminearum | other | 50 | pseudograminearum | reads-validated |
| F. sambucinum | other | 50 | sambucinum | assembly-supported |
| F. oxysporum | oxysporum | 50 | ../fusarium_oxysporum_fola | reads-validated |
| F. andiyazi | fujikuroi_complex | 36 | no panel | no panel |
| F. sacchari | fujikuroi_complex | 34 | sacchari | reads-validated (partial: 5/6 assembled loci typed) |
| F. sp. | solani | 31 | no panel | no panel |
| F. venenatum | other | 27 | no panel | no panel |
| F. falciforme | solani | 25 | falciforme | reads-validated |
| F. cf. solani | solani | 22 | solani | inferred (untested) |
| F. verticillioides | fujikuroi_complex | 21 | verticillioides | reads-validated |
| F. avenaceum | other | 21 | no panel | no panel |
| F. odoratissimum | oxysporum | 19 | ../fusarium_oxysporum_fola | reads-validated |
| F. solani | solani | 18 | solani | reads 1/2 pass; partial |
| F. keratoplasticum | solani | 17 | no panel | no panel |
| F. lateritium | other | 16 | sambucinum | partial |
| F. circinatum | fujikuroi_complex | 13 | circinatum | reads-validated (partial: 15/17 assembled loci typed) |
| F. petroliphilum | solani | 13 | no panel | no panel |
| F. xylarioides | fujikuroi_complex | 11 | xylarioides | reads-validated |
| F. cf. falciforme | solani | 10 | falciforme | inferred (untested) |
| F. boothii | graminearum_complex | 9 | graminearum | reads-validated |
| F. purpurascens | oxysporum | 9 | ../fusarium_oxysporum_fola | reads-validated (partial: 1/1 assembled loci typed) |
| F. euwallaceae | solani | 8 | no panel | no panel |
| F. annulatum | fujikuroi_complex | 7 | proliferatum | reads-validated (partial: 2/2 assembled loci typed) |
| F. meridionale | graminearum_complex | 7 | graminearum | assembly-supported |
| F. tardichlamydosporum | oxysporum | 7 | ../fusarium_oxysporum_fola | reads-validated (partial: 2/2 assembled loci typed) |
| F. phialophorum | oxysporum | 6 | ../fusarium_oxysporum_fola | partial |
| F. kuroshium | solani | 6 | no panel | no panel |
| F. nygamai | fujikuroi_complex | 5 | no panel | no panel |
| F. asiaticum | graminearum_complex | 5 | graminearum | assembly-supported |
| F. duoseptatum | oxysporum | 5 | ../fusarium_oxysporum_fola | inferred (untested) |
| F. cf. verticillioides | fujikuroi_complex | 4 | verticillioides | inferred (untested) |
| F. brachygibbosum | other | 4 | no panel | no panel |
| F. equiseti | other | 4 | no panel | no panel |
| F. ambrosium | solani | 4 | no panel | no panel |
| F. virguliforme | solani | 4 | no panel | no panel |
| F. gerlachii | graminearum_complex | 3 | graminearum | assembly-supported |
| F. cf. hostae | other | 3 | no panel | no panel |
| F. commune | other | 3 | no panel | no panel |
| F. musae | other | 3 | verticillioides | reads-validated (partial: 1/1 assembled loci typed) |
| F. nematophilum | other | 3 | no panel | no panel |
| F. grosmichelii | oxysporum | 3 | ../fusarium_oxysporum_fola | inferred (untested) |
| F. veterinarium | oxysporum | 3 | ../fusarium_oxysporum_fola | reads-validated |
| F. tucumaniae | solani | 3 | no panel | no panel |
| F. vanettenii | solani | 3 | vanettenii | partial |
| F. chuoi | fujikuroi_complex | 2 | no panel | no panel |
| F. subglutinans | fujikuroi_complex | 2 | verticillioides | partial |
| F. temperatum | fujikuroi_complex | 2 | verticillioides | partial |
| F. aethiopicum | graminearum_complex | 2 | graminearum | reads-validated |
| F. louisianense | graminearum_complex | 2 | graminearum | assembly-supported |
| F. acuminatum | other | 2 | no panel | no panel |
| F. citri | other | 2 | no panel | no panel |
| F. clavum | other | 2 | no panel | no panel |
| F. delphinoides | other | 2 | no panel | no panel |
| F. dimerum | other | 2 | no panel | no panel |
| F. hostae | other | 2 | no panel | no panel |
| F. irregulare | other | 2 | no panel | no panel |
| F. langsethiae | other | 2 | no panel | no panel |
| F. redolens | other | 2 | no panel | no panel |
| F. sporotrichioides | other | 2 | no panel | no panel |
| F. tricinctum | other | 2 | no panel | no panel |
| F. nirenbergiae | oxysporum | 2 | ../fusarium_oxysporum_fola | partial |
| F. tardicrescens | oxysporum | 2 | ../fusarium_oxysporum_fola | inferred (untested) |
| F. awan | solani | 2 | no panel | no panel |
| F. duplospermum | solani | 2 | no panel | no panel |
| F. mekan | solani | 2 | no panel | no panel |
| F. obliquiseptatum | solani | 2 | no panel | no panel |
| F. papillatum | solani | 2 | no panel | no panel |
| F. pseudensiforme | solani | 2 | no panel | no panel |
| F. rekanum | solani | 2 | no panel | no panel |
| F. tuaranense | solani | 2 | no panel | no panel |
| F. variasi | solani | 2 | no panel | no panel |
| F. acutatum | fujikuroi_complex | 1 | no panel | no panel |
| F. ananatum | fujikuroi_complex | 1 | no panel | no panel |
| F. anthophilum | fujikuroi_complex | 1 | circinatum | partial |
| F. bactridioides | fujikuroi_complex | 1 | circinatum | partial |
| F. begoniae | fujikuroi_complex | 1 | no panel | no panel |
| F. brevicatenulatum | fujikuroi_complex | 1 | no panel | no panel |
| F. cf. nygamai | fujikuroi_complex | 1 | no panel | no panel |
| F. concentricum | fujikuroi_complex | 1 | proliferatum | partial |
| F. denticulatum | fujikuroi_complex | 1 | no panel | no panel |
| F. ficicrescens | fujikuroi_complex | 1 | no panel | no panel |
| F. fracticaudum | fujikuroi_complex | 1 | circinatum | partial |
| F. globosum | fujikuroi_complex | 1 | no panel | no panel |
| F. guttiforme | fujikuroi_complex | 1 | verticillioides | partial |
| F. napiforme | fujikuroi_complex | 1 | no panel | no panel |
| F. pseudoanthophilum | fujikuroi_complex | 1 | no panel | no panel |
| F. pseudocircinatum | fujikuroi_complex | 1 | no panel | no panel |
| F. pseudonygamai | fujikuroi_complex | 1 | no panel | no panel |
| F. ramigenum | fujikuroi_complex | 1 | no panel | no panel |
| F. udum | fujikuroi_complex | 1 | xylarioides | assembly-supported |
| F. azukiicola | other | 1 | no panel | no panel |
| F. cerealis | other | 1 | pseudograminearum | reads-validated (partial: 2/2 assembled loci typed) |
| F. cf. aywerte | other | 1 | no panel | no panel |
| F. chlamydosporum | other | 1 | no panel | no panel |
| F. coffeatum | other | 1 | no panel | no panel |
| F. incarnatum | other | 1 | no panel | no panel |
| F. kistleri | other | 1 | no panel | no panel |
| F. rubellum | other | 1 | no panel | no panel |
| F. sinense | other | 1 | no panel | no panel |
| F. terricola | other | 1 | no panel | no panel |
| F. torreyae | other | 1 | no panel | no panel |
| F. cugenangense | oxysporum | 1 | ../fusarium_oxysporum_fola | inferred (untested) |
| F. foetens | oxysporum | 1 | ../fusarium_oxysporum_fola | inferred (untested) |
| F. kalimantanense | oxysporum | 1 | ../fusarium_oxysporum_fola | inferred (untested) |
| F. akasia | solani | 1 | no panel | no panel |
| F. bataticola | solani | 1 | no panel | no panel |
| F. brasiliense | solani | 1 | no panel | no panel |
| F. cucurbiticola | solani | 1 | no panel | no panel |
| F. cuneirostrum | solani | 1 | no panel | no panel |
| F. drepaniforme | solani | 1 | no panel | no panel |
| F. floridanum | solani | 1 | no panel | no panel |
| F. lichenicola | solani | 1 | no panel | no panel |
| F. metavorans | solani | 1 | no panel | no panel |
| F. neocosmosporiellum | solani | 1 | no panel | no panel |
| F. oligoseptatum | solani | 1 | no panel | no panel |
| F. phaseoli | solani | 1 | no panel | no panel |
| F. piperis | solani | 1 | no panel | no panel |
| F. suttonianum | solani | 1 | no panel | no panel |
| F. tenuicristatum | solani | 1 | no panel | no panel |
| F. tonkinense | solani | 1 | no panel | no panel |
| F. waltergamsii | solani | 1 | no panel | no panel |
| F. warna | solani | 1 | no panel | no panel |

## Validation (details in `validation.tsv`)

56 tests on 50 runs, first 0.4-2 M read pairs of each (60-90 MB per mate from `ftp.sra.ebi.ac.uk` via ranged HTTPS; 150-350 MB for old HiSeq runs), typed with the unmodified `reads-type` (k=31,
min-unique-run 100). **48 PASS, 8 NO_CALL, 0 wrong calls.** Of the 8 no-calls, 6 are cross-species tests expected to fail (listed in finding 2), 1 is a within-species failure (solani IISc-1, lineage not in panel), and 1 is the unresolved FPI024 conflict.
Expected alleles come from the assembly of the same strain (published or SRA-linked, MAT genes found with tblastn against the repo db proteins, `*` below) or, for ERR1288718 and SRR27384772 (*F. graminearum*),
from the species-level fact that all assembled *F. graminearum* loci are homothallic. Panel-source strains were also tested (marked in `truth_basis`), so those passes are not held out. Held-out passes: culmorum 270555; poae Fp133 and Fp023; pseudograminearum CS3270 and Fp22-2;
verticillioides BRIP53590 and BRIP14953; proliferatum 31915 and 31866; sacchari PB09 and GXUF-1; xylarioides IMI 392254; falciforme A02, UASB and CF00175; circinatum FL72_t3 and FK870_A2; and every fola test (the panel is from Fola).

Reads-type needs shared (flank) depth >= 1.0. With 0.5-0.7 M pairs of 150-bp reads shared depth was 0.9-1.4; ~1.5-2 M pairs gave 2-7. Plan for >=1.5 M pairs per strain (older 100-bp HiSeq runs need more).
Breadth on correct calls: 0.51-1.00 (lowest: fola on MRL8996 0.51, F. boothii/aethiopicum both-call 0.61-0.65, cerealis on Fps panel 0.56).

## Unverified or unresolved

- **SRR31408079 (*F. proliferatum* FPI024)**: assembly GCA_058129285.1 says MAT1-2, reads say something MAT1-1-like at shared depth 0.9 (no call). A run/assembly mismatch or a divergent MAT1-1 allele; not resolved.
- No truth beyond the same-strain assembly (no literature or phenotype check) for any run; the species-level `both` truth for ERR1288718 and SRR27384772 is not strain-specific.
- No read test at all for: fujikuroi panel (no same-strain ENA reads found), sambucinum panel, vanettenii panel, xylarioides MAT1-1 outside the source strain, the FGSC species asiaticum/gerlachii/meridionale/louisianense (assembly-supported only), *F. udum*, *F. nirenbergiae*, *F. phialophorum*.
- *F. poae* MAT1-1 and *F. culmorum*/*F. proliferatum* second alleles: in silico only (no run with a known MAT1-1 poae genotype was found).
- Species with several assemblies were not all scanned (first 1-16 each; some NCBI downloads failed: e.g. GCA_058139895.1, GCA_059629165.1, GCA_013618375.1, GCA_019055245.1, GCA_019425555.1), so rarer MAT alleles are missed (e.g. second *F. sacchari* MAT1-1 allele).
- Run labels (species) are the submitters'; *F. solani*/FSSC species labels in particular did not match MAT locus lineages in assemblies (e.g. solani and falciforme MAT1-2 alleles cluster together).
- Scratch download was ~7 GB of compressed reads (more than the brief's few GB); all deleted after typing.

## Provenance

All sequence is from NCBI Datasets genome assemblies (accessions in the panel table). Candidate loci were found by tblastn of the repo-db MAT proteins (5127_mo44, 5127_mo45, 5518_3639, 59765_880621a-1, 61372_suf959, i.e. GenBank AF100925.1, AF100926.1, AF318048.1, AB011379.2, AB011378.1)
against 366 Fusarium assemblies (up to 1-16 per species), loci classed MAT1-1 if a MAT1-1-1 hit >=220 aa at >=35% identity, MAT1-2 if a MAT1-2-1 hit >=85 aa at >=60% or >=150 aa at >=38%. Reference reads (ENA) were streamed, not stored.
