# Validation of review findings F3 and F4 (2026-09-28)

Read-only. Code at PR #9 head ad9be39 (polish-scope-cuts worktree).
Scripts: f3_fragment_loo.py, f4_analyze.py, scan.slurm (reuses
results/2026-09-27_pheromone_positional/scan_genome.py). Tables: f3_scores.tsv,
f3_summary.txt, f3_absolute.txt, f4_summary.tsv, out_neg/.

## F4: strict-CAAX finder, false-positive rate

Rule scored: >= 1 strict-CAAX ORF (C[VI][IV][AVMG], 20-130 codons, in-frame
Met) within +-10 kb (column T_10kb).

| Set | Flagged / total | Rate |
|---|---|---|
| Curated mating receptors (6 ground-truth genomes) | 6/9 | 0.67 |
| Other STE3 copies, same genomes | 0/25 | 0.00 (95% upper 0.137) |
| STE3 loci, 15 Agaricales CAAX-panel genomes | 30/142 | 0.21 |
| Random windows, same 15 genomes | 298/13,000 | 0.023 |
| Random windows, 58 Basidiomycota genomes (earlier set) | 1,437/50,000 | 0.029 |
| STE3 loci, rust (7 new + 4 earlier genomes) | 4/39 | 0.10 |
| Random windows, rust | 90/8,000 | 0.011 |
| STE3-like loci, Pezizomycotina (7 of 30 genomes had loci) | 1/10 | 0.10 |
| Random windows, Pezizomycotina | 155/7,000 | 0.022 |

- Chance component in the Agaricales panel: 142 x 0.023 = 3.3 expected flags
  among the 30 flagged receptor loci (about 11%). Scaled to the 118 gained
  calls: about 13 false admissions if non-mating receptor loci flag at the
  random-window rate.
- Upper bound: the only direct measure of non-mating receptor loci is 0/25,
  whose 95% upper limit is 13.7%, i.e. up to ~19 of 142 loci, or ~55% of the
  30 flags. The data cannot exclude a false-admission share near the review's
  ~38% estimate.
- Rust receptor loci flag at 10% vs 1.1% random; the rule does not recover
  known rust mating receptors (earlier exploration), so these flags are
  unexplained. Pezizomycotina is a weak negative control: its Ste3-like loci
  were mostly not found by the Basidiomycota query set (23/30 genomes gave no
  locus), and Pezizomycotina a-factor-like pheromones also end in CAAX.

## F3: classifier margin on fragments (leave-one-genus-out)

108 positive proteins (85 training + 23 Zygo), 30 genera; HMMs rebuilt without
each genus. Inputs: full protein, HMG box (PF00505 envelope +-5 aa), and 3
random 50-90 aa windows spanning the HMG-box midpoint. Negatives: 189
non-locus HMG paralogs (FastTree, status nonlocus, clade other_HMG), scored
with the shipped HMMs.

- sexM vs sexP typing on fragments is reliable: 0 wrong calls out of 540
  positive inputs (every margin > 0). At 25 bits: full 106/108, HMG box
  104/108, windows 294/324 called; below 25 they become undetermined.
- The margin does NOT separate MAT HMG genes from paralogs. At 25 bits,
  46/189 paralog full proteins, 37/189 paralog HMG boxes and 63/558 paralog
  windows exceed the floor. Raising the floor to 50 bits removes paralogs
  (2/189, 1/189, 0/558) but loses correct calls (96/108, 82/108, 204/324).
- Absolute best-HMM score separates better on full proteins (>=100 bits:
  positives 96/108, paralogs 9/189) but not on windows (>=80: 230/324 vs
  25/558; >=100: 132/324 vs 1/558).

## Conclusions

- F3: keep min_margin for TYPING (25 is conservative; no wrong-direction call
  was seen at any floor). It must not be read as evidence that the gene is a
  MAT gene. MAT-vs-paralog needs a separate criterion: locus context (flank
  genes) or an absolute-score floor on modelled proteins (~100 bits); on
  fragments neither margin nor absolute score is clean.
- F4: chance alone predicts ~11% false admissions; the 0/25 non-mating set is
  too small to exclude ~50%. Label CAAX-dependent calls unverified until a
  larger labelled set of non-mating STE3 loci (>=100) narrows the bound.

## Limits

- Ground truth: 9 mating and 25 non-mating receptor loci in 6 genomes.
- Pezizomycotina control is weak (see above). Rust flags are not labelled.
- Fragment windows are simulated, not real tblastn HSPs; the paralog set is
  the non-locus HMG copies in the scanned Mucoromycota-group genomes.
- rust.txt was rewritten once by an earlier selection run before the rust job
  started; the file on disk lists the 7 genomes actually scanned.
