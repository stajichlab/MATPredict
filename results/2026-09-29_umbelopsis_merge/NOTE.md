# curation-umbelopsis: deterministic rebuild and merge into PR #9

## Question
After shipping the deterministic PR #9 classifier (8eff88e), does
curation-umbelopsis, rebuilt deterministically with its three records, lose
calls outside the Lichtheimiaceae known gap? Curator merge rule (2026-09-29):
merge if all losses fall in that gap or in genomes shown to be misidentified
Circinella.

## Data and code version
- Branch rebased onto polish-scope-cuts 8eff88e; the guard and the old branch
  classifier commits dropped or superseded. Backup: backup/curation-umbelopsis-19e1d12.
- One deterministic full rebuild including the records 41833 (Umbelopsis Plus),
  44442 (Umbelopsis Minus), 13706 (S. racemosum NRRL 2496): commit f59353c.
  Frozen worktree run-f59353c.
- Baseline: the shipped classifier runs (results/2026-09-29_aligner_default/,
  byte-identical HMMs to 8eff88e).
- Inputs: Mucoromycota 293 (early_diverging list, taxid routing), LCG 621
  (--phylum Mucoromycota), Jena 64 scaffolds, Zygo 23 both inputs.

## Method
compare.py (calls matched by contig overlap; label/confidence changes; lost
calls with withholding reasons); scripts/check_record_selfcall.py; direct
pyhmmer scoring; blastp of annotated LCG flank proteins.

## Results
- Classifier: sexP 78 (27 genera), sexM 10 (8 genera); LOO 88/88, worst
  correct margin 13.6 bits; gate 98.4 bits (9/189 paralog negatives, 76/88
  held-out MAT at or above); P1 paralog class kept. Tests: 938 passed.
- Zygo 23: scaffold 23/23 locus and idiomorph; contig 23/23.
- Record self-call: 41833 Plus/high, 44442 Minus/high, 13706 Plus/medium.
- Umbelopsidaceae: 12/14 -> 13/14 called; 4 undetermined -> Minus/high
  (U. vinacea x2, WA50703 x2); U. nana gained (Minus/high); 3 U. isabellina
  and M5902 raised to high/medium.
- Mucoromycota 293: called 253 -> 253; lost 1, gained 1, changed 9.
  Label drop: M. griseocyanus CBS 116.08 Minus/high -> undetermined
  (margin 26.3 -> 17.8, under the 25-bit floor; the extra sexM training
  record lowers its sexM score). Not a lost call.
- LCG: 536 -> 535 called; gained 3 (M. ramannianus A-21216 Minus/high,
  Mucor sp. NRRL 1454 Minus/high, S. racemosum NRRL 2495 Minus/medium);
  lost 5 (below). U. ovata NRRL 13127T undetermined -> Minus/high.
- Jena: 61 -> 61; CBS206_69 undetermined -> Minus/low.
- Clean strain labels (17; disputed and misidentified excluded): 10 agree /
  3 disagree / 4 uncalled in both runs.

### Lost calls and causes
- Circinella minor GCA_016758965.1, LCG C. minor NRRL 1365 and CBS 143.56,
  C. umbellata NRRL 2417: withheld by the MAT-gene gate. Cause is the sexP
  MODEL, not HMM drift: the same exon translated scores 111.8 bits (shipped
  HMM) and 113.7 (branch HMM), but detect's branch model is a different
  protein (loser identity 32.8% cov 40 vs 31.6% cov 32), built after the new
  S. racemosum sexP reference changed polishing; it scores 86.4 < 98.4.
  Circinella is Lichtheimiaceae in BFD: inside the known gap.
- Mucor pusillus NRRL A-13674 (= Rhizomucor pusillus, Lichtheimiaceae):
  its sexP+sexM cluster on scaffold_161 is now `polish_capped: true` (the
  per-family cap of 6 skipped it; the new Umbelopsis flank references add
  competing clusters), so 0 genes are modelled and the modelled-gene bar
  withholds it; its fragment score is Plus 224.8. Inside the known gap.
- Rhizopus microsporus NRRL A-17693: a misidentified Circinella minor.
  tptA and glrA are 100% identical to C. minor NRRL 1365 and rnhA 97.5%;
  against R. microsporus NRRL 5546, rnhA 29.2% and glrA 75.7%. Added to
  results/2026-09-29_strain_labels_and_absidia/misidentified_strains.tsv.

### Merge decision
All five losses are Lichtheimiaceae or a misidentified Circinella, so the
curator's rule is met and the branch was fast-forwarded into polish-scope-cuts.

## Limits
- The Circinella model change and the M. pusillus polish-cap skip are
  side effects of adding references; they are general mechanisms that could
  affect other genera as more records are added.
- M. griseocyanus lost its label (not its call); the merge rule does not
  cover label drops.

## Curator decisions
Made: merge rule (2026-09-29). Open: whether the polish cap should protect a
cluster whose fragment score is strong (M. pusillus, 224.8); whether reference
additions that change polishing need a guard (Circinella).

## Files
compare.py, compare_output.txt, changes.tsv, umbelopsidaceae.tsv,
selfcall.tsv, zygo23_score.txt, reports_all.tar.zst, jobs.txt, run_*.slurm.
