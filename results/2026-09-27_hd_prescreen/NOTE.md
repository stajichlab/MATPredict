# HD prescreen before polishing (redHD, Rhodotorula) -- exploration

Question (curator, 2026-09-27): can an HMM screen redHD candidate clusters
before polishing, to cut the redHD cost without losing HD calls?

Data: run-02434ee reports, 221 genomes. 14,266 redHD candidate clusters: 222
positives (overlap the called HD locus), 608 admitted negatives (would be
polished), 13,436 non-admitted negatives (never polished). Record-strain
genomes excluded from evaluation (217 pos, 590 admitted neg left).

## Where the time goes (cProfile, median genome GCA_019059545.1, 4 threads)

| step | redPR only (ac33880) | + redHD (02434ee) | added |
|---|---|---|---|
| total | 10.5 s | 38.7 s | +28.2 s |
| tblastn localization | 1.4 s | 22.6 s | +21.2 s |
| polishing (exonerate) | 5.6 s | 11.1 s | +5.5 s |

A pre-polish filter can remove at most the +5.5 s (about 20% of the added cost).

## Separation (leave-one-species-out HD1/HD2 HMMs, 3-5 redHD records each, scored on tblastn HSP translations)

- HMM score: positives 159.5-540.8 (median 237.8); admitted negatives mostly
  no hit, max 251.1; threshold at the worst positive removes 589/590 admitted
  negatives.
- Plain tblastn best bitscore separates as well: positives 148-894, admitted
  negatives 27.7-150.0; threshold 148 removes 589/590.
- The one negative above threshold: Sporobolomyces pararoseus GCA_010758995.1
  (HMM 251, blast 150), uncalled HD-like cluster -- possibly a real HD locus
  that failed the bar; not checked.
- Admitted clusters per genome: median 3 -> 1 after the screen.

## Projection

Removing ~2 of 3 polished redHD clusters per genome saves roughly 2/3 of the
+5.5 s polishing, about 3.7 s of 38.7 s (~10%). Prescreen cost: ~20 ms per
cluster (HMM on existing HSPs; no extra search). The +21 s tblastn
localization (HD queries hitting genome-wide homeodomains; ~60 non-admitted
clusters per genome) is untouched by any pre-polish screen.

Files: collect.py, extract_hsps.py, score.py, summarize.py, clusters.tsv,
clusters_sampled.tsv, hsps.tsv, scores.tsv, summary.txt.
