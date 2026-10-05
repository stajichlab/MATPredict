# Basidiomycota PR test: the 1,360 CAAX-only receptor calls
Status: open (interim negative set; curator review pending)

## Question
Do the 1,360 Basidiomycota PR calls admitted only by the strict-CAAX scan
(`verification: unverified`) reflect real signal? Two tests: (1) the strict
CAAX rule on curated B-locus receptors and on non-mating STE3 copies; (2) how
many of the 1,360 calls the rule would give by chance.

## Data and code version
- Code: `origin/main` 8f07a22 (PR #29). The test uses the standalone scan
  `results/2026-10-05_caax_receptor_test/scan_genome.py` (copy of
  `results/2026-09-27_pheromone_positional/scan_genome.py`, only the STE3 query
  path changed). It uses the same strict motif `C[VI][IV][AVMG]` within 10 kb as
  `src/MATPredict/detect/caax.py`, but it finds receptors by miniprot of 1,094
  STE3 queries, not by the detect receptor hits.
- The 1,360 calls: `results/2026-10-03_basidiomycota_v060/loci.tsv` (1,472 PR
  loci; 1,360 `unverified`; 1,027 genomes; detection pass: 1,345 strict, 15
  relaxed). Counts by order: `unverified_calls_by_order.tsv`.
- Panel: 6 genomes of the 2026-09-27 study (34 STE3-like loci) plus 4 new
  Agaricomycete genomes (Trametes, Grifola, Russula, Heterobasidion; 30 loci).
- Sample: 502 genomes, 40 orders, up to 30 per order, drawn at random from all
  genomes of each order regardless of call status (`make_sample.py`, seed
  20261005). Genomes of 500 Mb or more were left out.

## Method
1. Find STE3-like loci per genome (miniprot). Flag a locus if a strict-CAAX ORF
   (Met 20-130 codons upstream of the stop) lies within 10 kb.
2. Panel: label loci mating if they overlap the curated record's receptor
   cluster; all other STE3-like loci are "other" (assumed non-mating). Rates
   with Wilson 95% intervals (`evaluate_panel.py`).
3. Sample: per genome, measure p = share of 1,000 random 20 kb windows with a
   strict-CAAX ORF (windows away from STE3 loci). Expected flagged loci by
   chance = n_loci x p. Per order: observed O, expected E, excess share
   (O-E)/O, bootstrap over genomes (`analyze_sample.py`).
4. Apply each order's E/O to its unverified calls (approximation: a call is a
   merged locus, a flagged locus is a miniprot locus).

## Results
### 1. Labelled panel (`panel_summary.txt`, `panel_loci.tsv`)
| Group | Flagged / n | Rate | Wilson 95% |
|---|---|---|---|
| Mating, independent records (6 genomes, all lineages) | 6/9 | 66.7% | 35.4-87.9 |
| Mating, independent, Agaricomycetes | 5/5 | 100% | 56.6-100 |
| Mating, CAAX-selected records (not independent) | 8/9 | 88.9% | 56.5-98.0 |
| Other, all lineages | 1/46 | 2.2% | 0.4-11.3 |
| Other, Agaricomycetes | 1/33 | 3.0% | 0.5-15.3 |
| Other, old panel (2026-09-27) | 0/25 | 0% | 0-13.3 |
| Other, new genomes | 1/21 | 4.8% | 0.8-22.7 |

- The one flagged "other" copy is in Russula nobilis (1 of 8 non-record copies).
- The 95% upper limit on the false-positive rate fell from 13.7% (0/25) to
  11.3% (1/46). The curator's bar of 100 labelled non-mating loci is not met.
- A first run labelled Trametes and Heterobasidion receptors as "other"
  (6/26 flagged). Cause: records use GenBank contig names, the GCF genomes use
  RefSeq names. The record sequence matches the genome exactly (record base 1 =
  NW_007360328.1:1,556,459 and NW_009258203.1:825,882). Fixed in
  `evaluate_panel.py`.

### 2. Chance model on 502 sampled genomes (`sample_summary.txt`, `sample_per_genome.tsv`)
- STE3-like loci 3,594; flagged 588; expected by chance 87.8; ratio 6.70;
  excess share (O-E)/O = 0.85, bootstrap 95% 0.82-0.87.
- Per order (excess share, bootstrap 95%, genomes sampled):
  Agaricales 0.89 (0.87-0.91, 30); Polyporales 0.91 (0.88-0.94, 30);
  Cantharellales 0.81 (0.72-0.88, 30); Boletales 0.79 (0.53-0.93, 30);
  Russulales 0.74 (0.58-0.85, 30); Trichosporonales 0.83 (0.71-0.89, 30);
  Tilletiales 0.97 (0.94-0.99, 30).
- Weak or no excess: Hymenochaetales O/E 1.64, excess 0.39 (0.08-0.70);
  Phallales O/E 0.80 (5 genomes); Agaricostilbales 0/20 flagged.
- Applying each order's chance share to the unverified calls: about 179 of
  1,359 calls (13%) are at chance level, about 1,180 are beyond chance
  (`sample_summary.txt`, last table). One call (Agaricostilbales) has no
  estimate. Agaricales: 702 calls, chance share 0.11 (76 calls). Hymenochaetales:
  15 calls, chance share 0.61.
- Genome level: 289 of 502 sampled genomes have a flagged locus; 234 of them
  carry an unverified v0.6.0 PR call, 55 do not. 213 have no flagged locus; 6 of
  them carry a call (the detect receptor hits differ from the miniprot loci).

### 3. Precursor homology (Hx) added to the CAAX rule (`hx_summary.txt`, `analyze_hx.py`)
Hx = tblastn hit (E <= 1) of the 31 curated precursors within 10 kb. Sample
statistics leave out 4 genomes of genera that supplied precursors (Cryptococcus,
Schizophyllum, Coprinopsis, Ustilaginales genera); 498 genomes, 3,560 loci.

| Rule | Mating, independent (n=9) | Other (n=46) | Sample observed / chance | Excess share |
|---|---|---|---|---|
| T (CAAX) | 6/9 (35-88%) | 1/46 (0.4-11.3%) | 578 / 86.9 | 0.85 (0.82-0.87) |
| Hx | 3/9 (12-65%) | 0/46 (0-7.7%) | 194 / 61.0 | 0.69 (0.61-0.74) |
| T and Hx | 3/9 (12-65%) | 0/46 (0-7.7%) | 149 / 1.8 | 0.99 (0.98-0.99) |

- T and Hx together is very specific but finds only about a third of the
  independent mating receptors. It is a high-confidence tier, not a replacement.
- Hx is sparse where precursors diverge: Cantharellales T 49 / Hx 0; Boletales
  23 / 1; Tilletiales 21 / 0; Auriculariales 45 / 2. Calls in those orders rest on
  CAAX alone. Hx is strong in Polyporales (T&Hx 45 vs 0.2 expected) and
  Cystofilobasidiales (25 vs 0.0).
- I did not apply T and Hx to the 1,360 calls: that needs the scan on all 1,027
  genomes (about 10 core-hours).

### 4. Receptor-protein pilot, 6 Agaricomycete species (`pilot_summary.txt`, `pilot_scores_A.tsv`)
Leave-one-species-out: for each held-out species, rank its STE3-like loci by
similarity or HMM score built from the other species' mating receptors (plus
curated REF receptors). 47 loci, 14 mating (Coprinopsis 3, Grifola 3,
Heterobasidion 3, Schizophyllum 2, Trametes 2, Russula 1). Variant A uses
Agaricomycete REFs only; variant B all 34 REFs.

| Score | Pooled AUC (A) | Bootstrap 95% | Pooled AUC (B) | Mean per-species AUC (A) |
|---|---|---|---|---|
| Nearest-neighbour similarity (mating minus other) | 0.59 | 0.42-0.76 | 0.59 | 0.70 |
| Profile HMM from mating receptors | 0.70 | 0.54-0.85 | 0.73 | 0.64 |
| Strict CAAX flag (reference) | 0.95 | 0.86-1.00 | 0.95 | 0.96 |

- Similarity does not separate mating from non-mating copies (interval includes
  0.5). The HMM shows a modest signal (mean score 429 mating vs 303 other).
- The top-ranked locus in a held-out species is a mating receptor in 2/6
  species by similarity and 1/6 by HMM. Neither picks the mating copy.
- The CAAX AUC is not an independent result: 4 of the 6 species' receptors were
  chosen from CAAX evidence.
- This agrees with the 2026-09-27 receptor tree: Agaricomycete mating receptors
  are not monophyletic.

## What changed in detection
Nothing. This is a measurement only.

## Limits
- "Beyond chance" means CAAX ORFs sit near STE3-like loci more often than in
  random windows. It does not show a locus is a mating receptor. The panel
  false-positive rate (1/46) suggests most flagged loci are mating, but the
  interval is wide (upper 11.3%).
- The "other" label is an assumption. C. cinerea has about 14 receptors in the
  genome (see the B43 record), so some "other" copies may be mating. This makes
  the false-positive rate conservative, not optimistic.
- The four new records were curated from CAAX-positional evidence. They are
  excluded from the sensitivity estimate (independent 6/9 only).
- Sensitivity rests on 9 independent receptors in 6 genomes. The scan misses
  Cryptococcus (precursor 45 kb away) and both Rhodotorula (CTxA motif).
- The chance model uses random windows away from STE3 loci. If STE3 loci sit in
  regions with a different short-ORF density, the expected count is biased.
  Not measured.
- Tandem receptor copies are not independent (Sebacinales: 489 loci in 12
  genomes). Bootstrap over genomes handles this only partly. Intervals for
  orders with fewer than 10 sampled genomes are unreliable.
- The sample covers 502 genomes, not all 1,027 called genomes. Genomes of 500
  Mb or more were left out.
- The step 4 extrapolation mixes loci and calls; treat 179 as approximate.
- Pilot (section 4): 6 species, 14 mating receptors; the "other" label is an
  assumption; intervals are wide. It tests only one training design.
- No promoter or motif analysis was done (no curated promoter set).

## Curator decisions
Made: 2026-10-04, interim negative set from the in-repo STE3 copies until a
literature set exists; pass criterion = hit rate on curated B-locus receptors
and false-positive rate on non-mating STE3, each with Wilson 95% intervals.
Curator view (2026-10-05, not a ruling): PR calls will be a weaker confirmation
of mating relatedness and may belong to a broader receptor class. The work may
focus on (a) detection at all, (b) detection in context around known genes and
genomic motifs or promoters, (c) protein alignments of candidates to classify
related groups. HMMs from true positives may not be useful; tried in a small
pilot (section 4). Approved: try the precursor-homology signal and a small
receptor-protein pilot.
Open: whether any order may lose the `unverified` label. Candidate evidence is
excess share with a lower bound of 0.7 or more (Agaricales, Polyporales,
Cantharellales, Trichosporonales, Tilletiales, Cystofilobasidiales,
Auriculariales). This is a proposal, not a ruling. I suggest Hymenochaetales stays
unverified. Needed: a literature set of non-mating STE3 loci (>= 100).

## Files
- `results/2026-10-05_caax_receptor_test/`: `panel_summary.txt`,
  `panel_loci.tsv`, `sample_summary.txt`, `sample_per_genome.tsv`,
  `sample_genomes.tsv`, `unverified_calls_by_order.tsv`, scripts
  (`scan_genome.py`, `evaluate_panel.py`, `make_sample.py`, `analyze_sample.py`,
  `run_new_genomes.sh`, `run_sample.sh`), raw tables
  `out_new.tar.zst`, `out_sample_loci_random.tar.zst`.
- Hx and pilot: `hx_summary.txt`, `analyze_hx.py`, `pilot_extract.py`,
  `pilot_test.py`, `run_pilot.sh`, `pilot_summary.txt`, `pilot_scores_A.tsv`,
  `pilot_scores_B.tsv`, `pilot_loci.tsv`, `pilot_proteins.faa` (SLURM job 29406597).
- SLURM jobs 29404658 (4 genomes, about 35 s each) and 29404666 (502 genomes,
  11 tasks, 12-23 min each, 0 failures).
- Earlier work: `analysis/2026-09-27_receptor-loci.md`.
