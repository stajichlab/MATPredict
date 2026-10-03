# LCG distant-genus placements: overrides recorded, for manual review (2026-10-02)

Status: decided (overrides recorded); the curator will review the trees and
taxa by hand later. This report lists what to look at.

## Question
13 LCG genomes sit next to a genus far from their file-name genus in the
nf_phyling BUSCO trees (B12). Does their rDNA agree with the tree? Should
MATPredict record them as taxon overrides?

## Ruling
J. Stajich, 2026-10-02: record all 13 as overrides, and keep this report as a
review entry for a later manual check of the trees and taxa.
- 12 rows: genus rank, status `confirmed` (tree + rDNA agree).
- Circinella_muscae_NRRL_1360: status `unconfirmed`, not renamed. Its nuclear
  markers say Benjaminiella; its rDNA says Circinella. Possible mixed sample.
- All 13: `exclude from species-level scoring`.
- File: `db/taxon_overrides.tsv` (now 18 rows).

## Data
- Held-out set: LCG. No training, curation or classifier use.
- Trees (FastTree, SH-like local supports, 941 taxa):
  - `results/2026-10-01_lcg_name_check/nf_out/protein/buildtree/fungi_odb12/fasttree/protein-lcg_namecheck_v1-taxa_941.fungi_odb12.fasttree.support.treefile`
  - `results/2026-10-01_lcg_name_check/nf_out/protein/buildtree/mucoromycota_odb12/fasttree/protein-lcg_namecheck_v1-taxa_941.mucoromycota_odb12.fasttree.support.treefile`
- Tip names: `LCG__<genome>`; references are `JENA__` and `BFD__` tips.
  Taxon table: `results/2026-10-01_lcg_name_check/taxa.tsv`.
- Flags and supports per genome: `results/2026-10-01_lcg_name_check/flags.tsv`.
- rDNA: `results/2026-10-01_lcg_name_check/its_check/NOTE.md`,
  `its_verdicts.tsv`, `regions.*.fa`.
- IQ-TREE subtree (325 tips, LG+F+R8, UFBoot 1000):
  `results/2026-10-01_lcg_name_check/iqtree_subtree/`.
  - mucoromycota_odb12: job 29344888, running at the time of writing.
  - fungi_odb12: job 29344889 cancelled by the curator after 10.4 h at search
    iteration 20 (projected more than 35 h to finish). The FastTree tree is the
    fungi_odb12 result. The last best ML tree (no supports) is at
    `iqtree_subtree/fungi_odb12/subtree.fungi_odb12.treefile`.

## Results
| Genome | Tree genus (nearest ref clade) | SH-like fungi / mucoro | rDNA | Status |
|---|---|---|---|---|
| Phycomyces_blakesleeanus_NRRL_1555 | Syncephalastrum (S. contaminatum) | 0.957 / 1.000 | ITS2 + LSU: Syncephalastrum; no full ITS | confirmed |
| Phycomyces_blakesleeanus_NRRL_1556 | Gilbertella (G. persicaria) | 0.999 / 1.000 | ITS 100% G. persicaria (type) | confirmed |
| Phycomyces_nitens_NRRL_2700 | Gilbertella (G. persicaria) | 0.999 / 1.000 | ITS 100% G. persicaria (type); = NRRL 1556 | confirmed |
| Pilaira_anomala_NRRL_2527 | Syncephalastrum (S. contaminatum) | 1.000 / 0.999 | ITS2 + LSU: Syncephalastrum; no full ITS | confirmed |
| Pilaira_anomala_RSA_1997_Plus | Cunninghamella (C. blakesleeana) | 1.000 / 1.000 | ITS 99.6% C. echinulata | confirmed; same-strain JGI genome is P. anomala |
| Thamnidium_elegans_NRRL_2467 | Syncephalastrum (S. racemosum) | 1.000 / 0.997 | ITS 100% S. racemosum; 2nd ITS copy unmatched | confirmed; same-strain JGI genome is T. elegans |
| Syncephalastrum_racemosum_NRRL_1506 | Phycomyces (P. blakesleeanus) | **0.565** / 1.000 | ITS 100% P. blakesleeanus (type) | confirmed; weak fungi_odb12 support |
| Syncephalastrum_racemosum_NRRL_1623 | Pilaira (P. anomala) | 1.000 / 1.000 | ITS 100% P. anomala | confirmed; tree reference has Mucor rDNA |
| Syncephalastrum_sp._NRRL_1485 | Dichotomocladium (D. elegans) | 1.000 / 1.000 | LSU 100% D. robustum (type); no full ITS | confirmed |
| Rhizomucor_pusillus_NRRL_2543 | Phycomyces (P. nitens) | 1.000 / 1.000 | ITS 100% P. nitens; 2nd ITS copy unmatched | confirmed |
| Cunninghamella_japonica_NRRL_2464 | Actinomucor (A. elegans) | 1.000 / 1.000 | ITS 100% A. elegans | confirmed |
| Circinella_tenella_NRRL_A-23557 | Lichtheimia (10-tip clade) | 1.000 / 1.000 | ITS2 + LSU: L. brasiliensis (type); no full ITS | confirmed |
| Circinella_muscae_NRRL_1360 | Benjaminiella (B. poitrasii) | 1.000 / 1.000 | ITS 99.7% C. muscae | **unconfirmed**: rDNA and markers disagree |

Topology check on the interim IQ-TREE best trees (no supports, not final):
the same tree genus for all 13. All FastTree flags: mucoromycota_odb12 101 of
101 still flagged with the same genus; fungi_odb12 101 of 102. The exception is
Pirella circinans var. volvogradensis RSA 2566, which became consistent in the
fungi_odb12 interim tree.

### Effect on held-out scoring
`results/2026-10-01_heldout_rerun/analyze.py` rerun (backups `*.pre_b12_its`):
- LCG species summarised 200 -> 199; excluded genomes 5 -> 18.
- Circinella tenella leaves the species table (its only genome is excluded).
- Lower genome counts for C. muscae, C. japonica, P. blakesleeanus, P. nitens,
  P. anomala, R. pusillus, S. racemosum and Syncephalastrum sp.
- No other line of `summary.txt` changed.

## Manual review checklist (curator, later)
1. Open each genome's tip in both FastTree trees. Check the sister clade and the
   SH-like value against the table above.
2. Syncephalastrum_racemosum_NRRL_1506: fungi_odb12 support is 0.565. Check
   whether the mucoromycota_odb12 UFBoot result (when it finishes) supports it.
3. Circinella_muscae_NRRL_1360: decide between mixed sample, misnamed culture,
   or an unassembled Benjaminiella rDNA. Read coverage would separate these.
4. Thamnidium_elegans_NRRL_2467 and Rhizomucor_pusillus_NRRL_2543: each has a
   second ITS copy that matches nothing above 89%. Not tested.
5. Syncephalastrum_racemosum_NRRL_1623: the tree reference it pairs with (Jena
   Pilaira anomala CBS 695.68) has Mucor saturninus rDNA. Check that reference.
6. The 2 same-strain mismatches suggest LCG sample swaps or contamination. The
   JGI genomes (Pilano1, Thaele1) of the same strain ids have the expected rDNA.
7. Reference genomes with foreign rDNA: see the side note in
   `analysis/open-questions.md`.

## Limits
- rDNA in these assemblies is collapsed to 1-2 copies; 4 of 13 have no full ITS.
- rDNA shows what rDNA is in the assembly, not which organism the nuclear
  genome is from (NRRL 1360 shows they can differ).
- Genus rank only. Species names from ITS hits are in the basis column, not
  recorded as identities.
