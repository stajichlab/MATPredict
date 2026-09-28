# Mucoromycota curation + secondary-call guard (2026-09-27)

Branch curation-umbelopsis (local), rebased on polish-scope-cuts 44567fc.
Run: frozen run-2e9aa97 vs split-locus run 4174440 (results/2026-09-27_split_locus/).

## Commits
- a2f1da9 classifier rebuilt after rebase (LOO 87/87, worst 13.0 bits)
- bbcb524 guard: withhold a call the classifier ran on and left `undetermined`
  when the same family has a determined call in the genome (not homothallic).
  Replay (replay_guard.txt): a broad enum-family guard would have withheld 8
  Saccharomyces silent-cassette calls; this narrow form withholds 0 outside
  Mucoromycota (Cryptococcus 243, Dothideomycetes 99, Serinales 2368, cap6 567).
- 2e9aa97 tier-2 S. racemosum NRRL 2496 MAT Plus record (13706_nrrl-2496_MAT_Plus);
  classifier rebuilt (LOO 88/88, worst 13.5 bits); B1.9 logged.

## Lichtheimiaceae / Syncephalastraceae (revised rule: margin > 25 + Mucorales order)
- S. racemosum NRRL 2496: unannotated sexP ORF MCGN01000004.1:1,755,515-1,756,447 (+),
  310 aa, HMG box aa 112-175, classifier sexP 179.1 vs sexM 38.9 (+140), 69 bp
  upstream of rnhA ORY97819.1 (sexP->rnhA, Mucorales orientation); tptA/algA/glrA
  elsewhere (partial order). Same HMG 100% in B6101 (called Plus beside rnhA).
  Final ML tree: called S. racemosum loci outside sexP clade (tree too weak).
- Lichtheimiaceae: FAIL gene-order test. L. ornata GCF_029851405.1 strongest
  sexP-scoring HMG XP_058338267.1 (+103.9) and R. pusillus CBS 183.67
  XP_069244530.1 (+74.8) sit on scaffolds with no tptA/rnhA/algA/glrA within
  30 kb; no HMG near their rnhA or tptA. No record built.

## Scan results (293 genomes; Zygo 23/23 locus+idiomorph on scaffolds and contigs)
Genomes with any call 247 -> 247. 205 genomes margin-only changes.
| group | calls | genomes gained/lost | guard-withheld |
| Umbelopsidaceae | 12 -> 15 | 2 / 0 | 6 |
| Lichtheimiaceae | 3 -> 2 | 0 / 1 | 0 |
| Syncephalastraceae | 4 -> 4 | 0 / 0 | 0 |
| other | 242 -> 241 | 0 / 1 | 11 |

Guard withheld 17 second calls: 6 Umbelopsis, 5 Apophysomyces, 4 R. arrhizus
(new calls enabled by the records), Benjaminiella poitrasii, M. ardhlaengiktus.
One spurious-looking second call survives: U. ramanniana AG NW_026252095.1
Minus/medium, margin 25.6 (just above the 25-bit floor).

Umbelopsis: U. vinacea x2, WA50703 x2, U. nana -> Minus (medium/high);
U. isabellina B7317 and MPG-14A low -> high; WA0000067209 and M5902 low -> medium;
AD052 lone call undetermined/medium (24.6, kept: nothing to defer to).

Side effects still present:
- gzUmbRama1 and U. ramanniana AG high -> medium (distant rnhA hit in cluster).
- Mucor griseocyanus Minus -> undetermined (21.7; classifier retrain).
- R. pusillus FCH_5_7 and Absidia glauca AG_v1 lose their Minus/medium partial
  calls: best cluster now has 0 modelled genes. Not diagnosed.
- R. stolonifer PRFJ01/PRFJ02 gain Minus/medium partial_locus (26.1). Unverified.
- S. racemosum NRRL 2496 itself is still called only at the weak glrA/rnhA
  cluster on MCGN01000001.1 (undetermined 11.8); the record's own locus on
  MCGN01000004.1 is not called. Not diagnosed.
- S. monosporum x2 margin 10.3 -> 8.0 (still undetermined).
