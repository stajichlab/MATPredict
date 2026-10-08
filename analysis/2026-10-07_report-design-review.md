# Per-genome HTML/PDF report (`matpredict report genome`): build and design review
Status: open (implemented on branch `report-genome`; awaiting curator review)

## Question
Can one `detect` run be turned into a report that a biologist reads in five seconds and a curator can audit, that
works offline, on screen (desktop, phone, light, dark) and prints to PDF? Spec:
`docs/superpowers/specs/2026-10-07-genome-report-and-reads-intake-design.md`, Part 1.

## Data and code version
- Code: branch `report-genome`, new package `src/MATPredict/report/` (`html.py`, `svg.py`, `pdf.py`, `cli.py`) and
  `src/MATPredict/detect/provenance.py` (the new `run` block of `detection_report.yaml`).
- Inputs: one real report (`results/2026-10-06_fola_50a`, old format, no `run` block) and six synthetic reports
  (`tests/report/make_fixtures.py`): Mucoromycota with the HMM classifier, two idiomorphs (unlinked), not searched,
  nothing called with an assembly gap and withheld loci, Basidiomycota HD with subloci plus an unconfirmed PR call in
  a receptor array, and hostile strings (`<script>` sample name, `<img onerror>` contig).
- Engines: Chrome 1xx headless (macOS) and WeasyPrint 70.0 (Pango from Homebrew).
- Reviewer: an independent agent briefed as a senior web, information and scientific-visualization designer; review
  only, no edits. Three rounds on rendered HTML, screenshots and PDFs from both engines.

## Method
1. Render every fixture to HTML, a Chrome PDF and a WeasyPrint PDF; screenshots at 1200 px (light, dark) and inside a
   true 390 px iframe (headless Chrome will not size a window below about 500 px, so the first phone shots were wrong).
2. Reviewer lists findings P0 (must fix), P1, P2, each with a concrete fix; I apply them; the reviewer verifies
   against new renders.

## Results
Round 1 (first version): "not yet high quality". 6 P0, 12 P1, 6 P2. The P0s:
- two-idiomorph result read as two equal answers (chips first, flag below, cause list repeated);
- an old report showed the folder name ("fixtures") as the sample and the renderer's version as the run's;
- the figure did not implement its own encoding (no dashed outline for models that differ, legend key matching
  nothing, hit-only genes reduced to grey hatching, locus span drawn like a gene);
- figure unreadable on a phone and one table overflowing its card;
- ids broken mid-token (`578113_sxlc146_MAT_M / AT1-2`), print fonts not inherited in tables;
- wrong statements: a PR call with no flanking gene labelled "Complete locus: MAT gene(s) with their flanking genes",
  a provenance sentence pointing at citations that are not on the page, an unreadable CAAX warning.
Also found by me before review: WeasyPrint drew the figure black (it ignores page CSS inside SVG) and collapsed
the CSS-grid fact lists; every figure reused `id="hatch"`; the sample name in the `@page` header could close the
`<style>` element; Chrome on macOS sometimes did not exit after writing the PDF, or exited before the file landed.

Round 2: "close". Verified most fixes; one new P0 that the rewrite introduced: with one idiomorph candidate the
page said "No other idiomorph scored" and hid `idiomorph_margin` (46.7), while the table below showed the set-aside
MAT1-1-3 hit. In the pipeline (`pipeline.py`, idiomorph margin block) that margin is the narrowest overlap
resolution in identity points and it caps the confidence tier, so the page contradicted itself. The synthetic
fixtures had scores that did not match their own bitscores, so the tests could not catch it. Six P1s (phone figure
still too small, gene names breaking at hyphens, two-idiomorph wording, filled legend keys, raw strings, WeasyPrint
running header naming the wrong locus) and eight P2s.

Round 3: "high quality on screen in every case"; overall once the printed region string is fixed. Remaining then:
P1 region string cut off in print (nowrap plus a scroll box; fixed: wraps in print), P1 two legend keys looked the
same (fixed: grey fill with dark dashed outline vs pale fill with dotted outline), P2 phone figure opened at its
left end (fixed: scrolls to the first core gene on load), P2 blank quarter of page 1 in print (open: WeasyPrint
still moves the score box whole; probably its flex fragmentation). Final fixes checked by me on the renders, not
by a fourth review round.

What the report now does:
- Result first, in plain words, for each outcome: "Mating type MAT1-2, high confidence" with span, contig and gene
  order; "Two idiomorphs found (Plus and Minus): needs review" with supported and not-ruled-out causes; "No MAT locus
  called" with a next step; "Not searched" with the way out.
- One card per locus: SVG figure (role by fill and lightness, core drawn taller; dashed dark outline where exonerate
  and miniprot disagree, with the alternate model under it; pale dashed for hit only; strand arrows, exons; "called
  locus" bracket; contig end cap or distance), legend with only the keys used, flags, facts, samtools region string,
  idiomorph evidence with the margin named for what it measures, gene table (model, identity with coverage, E-value
  with bits, position, exons, reference record; cross-hit rows marked "not counted").
- What was searched, withheld candidates (collapsed on screen, open in print), provenance, glossary.
- Print: page header (sample, current locus in WeasyPrint, title), page numbers, no breaks inside a locus summary,
  light colours in print from a dark screen.
- PDF: 3 s per report with Chrome (about 0.5 MB, fonts embedded), WeasyPrint about 0.1 MB. 0 failures in 12 repeated
  Chrome renders after the wait fix.
- Tests: 30 in `tests/report/` and `tests/detect/test_provenance.py`; full suite on macOS: the 13 failures present
  before the change (missing Linux binaries, viz extras) and no new ones.

## Limits
- The fixtures are synthetic except the Fola 50a report; the Basidiomycota, Mucoromycota and no-call layouts have
  not been seen on real campaign reports (those live on the HPCC). Render a sample of real reports before release.
- The reviewer is a model briefed as an expert, not a person; a curator's read is still needed.
- Not done: a figure for an assembly gap at an uncalled family's locus (text only); a reads-type report (spec Part 2).
- WeasyPrint in the image is checked only by the Docker CI smoke test (no Docker daemon on the build machine);
  the pixi environment pins 69.x while the local renders used 70.0.
- `not_detected` reasons are the pipeline's own sentences, shown as written (written for curators, not plain words).

## Curator decisions
Decided 2026-10-07: (1) `detect` writes report.html by default, opt out with `--no-html` / `MATPREDICT_HTML=0`
(set in the batch scripts). (2) PDF engine: WeasyPrint in the pixi environment and image (25 packages, 10.3 MB
download, against about 300 MB for chrome-headless-shell 154 and its libraries); revisit if insufficient. Open: (3) Withheld loci collapsed (current) or hidden by default?

## Files
- Code: `src/MATPredict/report/`, `src/MATPredict/detect/provenance.py`; tests `tests/report/`,
  `tests/detect/test_provenance.py`.
- Renders and screenshots of all three rounds: session scratchpad only (not tracked); regenerate with
  `matpredict report genome --run tests/report/fixtures/<case>.yaml --out X.html --pdf X.pdf`.
