# 39-bit flank-carried floor and fragment classifier, implemented (2026-09-27)

Code: branch flank-bitscore, commits 51ec961 (bitscore floor) and 8d80bed
(fragment classifier), pushed onto polish-scope-cuts. Frozen worktree
run-8d80bed. Tests: 809 pass (was 790).

## Rulings implemented
1. flank_carried: strongest core hit = highest bitscore; kept at low only if
   >= 39 bits (FLANK_CARRIED_MIN_BITSCORE; per-family
   `flank_carried_min_bitscore`, no family overrides) and within the window.
   GeneEvidence and the report now carry each gene's best bitscore.
2. When no core protein is modelled, the classifier scores the tblastn HSP
   translations (3 frames, genome code), min_margin kept; the verdict reports
   `classifier_input: hsp_fragment` (else `model`, or `mixed`).

## Zygo 23 (both inputs)
scaffold 23/23 locus, 23/23 idiomorph; contig 23/23, 23/23.

## Mucoromycota scan, 8d80bed vs 4457c3a (293 genomes)
245 -> 249 calls: 243 same, 4 new, 2 idiomorph changes, 0 lost
(compare_vs_4457c3a.tsv).
- New (the flank floor): 4 Umbelopsis loci, all `undetermined`, low,
  partial_locus -- U. vinacea WA0000051536 (margin 13.7), U. vinacea
  gzUmbVina2 (9.4), U. sp. WA50703 x2 (16.8). All lean sexM (Minus) but under
  the 25-bit min_margin. U. nana is still withheld (33.5 bits < 39).
- Changed (the fragment classifier): Umbelopsis sp. M5902 Minus -> Plus
  (Plus 103.7 vs Minus 27.9); Mucor hiemalis gzMucHiem1 Plus -> Minus
  (Minus 97.2 vs Plus 50.9). Both as predicted.
- 9 calls are now fragment-decided (fragment_labels.txt): the 6 above plus 3
  U. isabellina, unchanged Plus (margins 68.5-92.2).

## Audit replay on the real code (replay_bitscore.txt)
Audit's own tblastn bitscores, run-8d80bed function and db.
- Ascomycota: 39-bit keeps 21/89 (E rule 9/89) -- matches the evaluation's
  prediction (21/89). Nothing kept under E is withheld now.
- Mucoromycota group: 12/19 (E rule 8/19); +4 = the Umbelopsis loci above.

## The 12 newly kept Ascomycota calls
All core genes unmodelled; best in-locus bitscore 39.3-46.6.
| Genome | Species (order) | Core hit (bits, id) | Gene order | Verdict |
|---|---|---|---|---|
| GCA_031142775.1 | Microdochium paspali (Xylariales) | MAT1-2-1 46.6, 24% | COX13-APN2-[MAT]-SLA2, 1.7-1.9 kb from APN2 | real (synteny); weak core, rank 2 genome-wide |
| GCA_060307455.1 | Xylariales sp. XT01 | MAT1-2-1 39.3, 25% | COX13-APN2-[MAT]-SLA2 | real (synteny); weak core |
| GCA_001600575.1 | Didymobotryum rigidum (Xylariales) | MAT1-1-2 44.7, 36% | SLA2-[MAT1-1-2, -1-3]-APN2-COX13 | real (audit's named case) |
| GCA_030574475.1 | Trigonopsis californica | MAT1-1-1 44.7, 37% | APN2-[MAT1-1-1]-SLA2, rank 1 | real (synteny) |
| GCA_031125595.1 | Saccharomycopsis capsularis (Ascoideales) | MAT1-1-1 42.0, 43% | SLA2-[MAT1-1-1] (0.7 kb)-APN2, rank 1 | real (synteny) |
| GCA_030564325.1 | Blastobotrys capitulatus (Dipodascales) | MAT1-1-1 40.8, 26% | APN2-[MAT1-1-1]-SLA2, rank 1 | real (synteny); weak core |
| GCA_037042115.1 | Wickerhamiella sp. (Dipodascales) | MAT1-1-1 39.7, 34% | APN2-[MAT1-1-1]-SLA2, rank 1 | real (synteny) |
| GCA_001939105.2 | Sugiyamaella xylanicola (Dipodascales) | mata1 41.2, 45% | sla2-[mata1]-apn2 (MATyl), rank 1 | real (synteny) |
| GCA_012184355.1 | Dactylellina cionopaga (Orbiliales) | MAT1-1-1 40.8, 32% | [MAT1-1-1]-APN2-COX13, rank 1 | real (Orbiliales pattern) |
| GCA_047671935.1 | Hyalorbilia sp. (Orbiliales), MAT | MAT1-2-1 42.4, 45% | SLA2-[MAT1-2-1]-APN2-COX13, rank 4 | real (synteny) |
| GCA_047671935.1 | same genome, MATtub | MAT1-2-1 42.4, 42% | same locus | duplicate of the MAT call |
| GCA_047716655.1 | Hyalorbilia oviparasitica (Orbiliales) | MAT1-2-1 42.4, 45% | COX13-APN2-[MAT1-2-1] (4.4 kb), rank 4 | real (Orbiliales pattern) |

11 distinct loci, all with the MAT gene between or beside the expected
anchors (SLA2/APN2/COX13); the 12th is the same Hyalorbilia locus called by a
second family (the cross-family duplicate issue already noted in the audit).
No noise found. Limit: "real" rests on synteny; three core hits are weak
(24-26% identity).

## The two Basidiomycota HD calls (not re-run; full-run reports)
- Pleurotus tuoliensis GCA_003243755.1, QCWT01000017.1:489817-497904: HD1
  42.7 bits (29%) with MIP1 and beta_fg in the same cluster -> real (both
  Agaricomycete HD anchors).
- Coriolopsis trogii GCA_007896425.1, SZVL01000001.1:2429139-2438918: HD1+HD2
  48.1 bits (28%) with MIP1 and beta_fg -> real.

## Limits
The floor rests on 8 real loci; the 4 new Umbelopsis calls are undetermined,
not labelled; the Ascomycota verdicts rest on synteny, not truth.
