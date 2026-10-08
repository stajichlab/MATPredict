# MATPredict: per-genome HTML/PDF report, and reads intake for the on-demand service

Date: 2026-10-07
Status: draft for review. Part 1 phases 1-2 (provenance block, `report genome`, `--pdf`, design review) are
implemented on branch `report-genome` (2026-10-07); Part 2 is not implemented.
Builds on `2026-10-04-packaging-service-and-reports-design.md` (workstreams B
and C) and `2026-10-05-reads-typing-and-population-locus-discovery-design.md`
(Tool A). That document's choice of service host (Galaxy, NRP web app,
Nextflow) is not repeated here; this one specifies the two pieces every host
option needs first: a report a person can read, and a way to get reads plus an
organism into `reads-type` or `detect`.

## Context

The curator asked (2026-10-07) how a small web service could take a genome, or
raw reads (SRA run or upload) plus the organism, and return a mating-type
report; whether it belongs in Galaxy; whether NRP (nrp.ai) Kubernetes pods
can run it; and what other generic hosts exist. The answer, in short:

- Every host option runs the same container (`ghcr.io/stajichlab/matpredict`,
  offline taxonomy, embedded database), so the host choice is a deployment
  decision made later (2026-10-04 spec, B2).
- What no host provides, and all of them need, is (1) a readable per-genome
  report (README goal 4; today the output is YAML and GFF3 only), and (2) a
  reads intake: organism to taxid, SRA to reads, taxid to a `reads-type`
  panel, and a rule for when reads must be assembled and sent to `detect`
  instead.

## Measured facts this design depends on

| Fact | Value | Source |
|---|---|---|
| `detect` wall time, one core | 1.5 min (Phycomyces) to about 1 h (>500 Mb genomes) | README; 2026-10-04 spec |
| `detect` memory | about 0.5-1 GB per genome | 2026-10-04 spec |
| `reads-type` wall time | about 1.5 min per strain at 8M reads, pure Python, 1 CPU | `analysis/2026-10-06_fola-reads-type.md` |
| `reads-type` agreement | 145/148 Fola strains with samtools breadth calls; 3 refusals (`low_depth`), 0 opposite calls | same |
| `reads-type` divergence limit | a SNP every 40 bp gives `none`; every 60 bp is called correctly (synthetic) | same |
| `reads-type` panel | user supplies `--idiomorph NAME=FASTA`; no builder from `db/` | `src/MATPredict/reads/cli.py` |
| Report provenance in `detection_report.yaml` | none: no input file name, taxid, organism, MATPredict version, database version, date or parameters | `src/MATPredict/detect/report.py:write_detection_report` |
| Database release stamp | not implemented (`curate-db release` raises NotImplementedError) | `src/MATPredict/db/cli.py` |
| Assembly path under-reports hybrids | A. fumigatus putative hybrids: MAT1-1 on a 2.4 kb contig withheld or unreported | 2026-10-05 spec, results 2026-10-05 |
| Drawing code | pyGenomeViz for curated records only; no drawing of a `detect` call | `src/MATPredict/db/draw.py` |

## Part 1: `matpredict report` (per-genome HTML, printable to PDF)

### 1.1 Goals
1. One self-contained `.html` file per run: no network, no external CSS, JS or
   fonts; opens from disk, e-mail or a Galaxy history; about 50-300 kB.
2. A first screen that answers "what mating type is this genome, and how sure
   are we" in plain words, for a reader who does not know the pipeline.
3. Every number behind the answer is on the page, so an expert can audit the
   call without opening the YAML.
4. Prints to PDF cleanly from any modern browser ("Save as PDF"), and from the
   command line (`--pdf`) without a browser.
5. Reads the existing `detection_report.yaml` and `detected_loci.gff3`; old
   reports (no provenance block) still render, with "not recorded" in place of
   the missing fields.

### 1.2 Command
```
matpredict report genome --run DIR [--out FILE.html] [--pdf FILE.pdf]
                         [--sample NAME] [--title TEXT]
matpredict report reads  --tsv reads_type.tsv [--sample NAME] [--out ...] [--pdf ...]   # Part 2
```
`--run` is a `detect --out-dir`. `--sample` overrides the sample name when the
report has no provenance. `report aggregate` (2026-10-04 spec, D) is not part
of this phase.

Decided (curator, 2026-10-07): `detect` writes `report.html` on every run by default; `--no-html` or
`MATPREDICT_HTML=0` opts out (batch scripts use the variable, which older worktrees ignore); `--pdf` adds the PDF.
A report error never fails the run.

### 1.3 Provenance block (small change to `detect`)
Add a `run` mapping at the top of `detection_report.yaml`, written
unconditionally:

```yaml
run:
  matpredict_version: 0.6.1
  database: {root: /app/db, content_sha256: 3f2a..., records: 150}
  taxonomy_source: local NCBI taxonomy snapshot 2026-10-01
  genome: {file: genome.fna, sha256: ..., contigs: 812, length_bp: 48213377, n50: 1203321}
  proteins: null
  taxid: 5507
  organism: Fusarium oxysporum
  lineage: [Fungi, Dikarya, Ascomycota, Pezizomycotina, Sordariomycetes, Hypocreales, Nectriaceae, Fusarium]
  parameters: {phylum: null, exhaustive: false, genetic_code: null, max_polished_clusters_per_family: 6, ...}
  started: 2026-10-07T14:03:11Z
  wall_seconds: 162
```
- `content_sha256` is a hash over the sorted database files `detect` reads,
  until `curate-db release` exists. It lets two reports say whether they used
  the same database.
- The genome sha256 and the assembly statistics (contigs, length, N50) cost
  one extra pass over the FASTA; the statistics put "uncalled" in context (a
  fragmented assembly is the most common reason for a lost locus).
- Additive: no existing key changes, so the regression check and rollout
  scripts are unaffected.

### 1.4 Page structure
In order, top to bottom:

1. **Header.** Sample name, organism (italic binomial) and taxid, short
   lineage, date, MATPredict version, database hash, taxonomy snapshot.
2. **Result.** One sentence per searched family and one badge, for example
   "MAT1-2, high confidence, 1 locus on NODE_86 (12.1 kb), complete
   (MAT1-2-1 between SLA2 and APN2)". Plain-language sentences for the other
   outcomes:
   - nothing called: the reason (no hits, withheld by a named gate, assembly
     gap at the expected locus) and "an uncalled genome is weak evidence that
     the locus is absent";
   - `not_searched`: why, and the `--phylum` / `--exhaustive` way out;
   - two idiomorphs: the arrangement and every possible cause (homothallism,
     mixed culture, merged diploid assembly), never a verdict;
   - `zygosity: unknown`: the assembly may have collapsed a heterozygous locus;
   - routing error, or a call outside the genome's phylum (`verification`).
3. **Per-locus card** (one per `detected` entry):
   - gene-order figure (1.5);
   - key facts: family, idiomorph and margin, confidence, locus class,
     detection pass, coordinates and span, distance to contig ends, merged
     subloci, receptor array (PR), flags;
   - idiomorph evidence: candidate scores, classifier scores and margin
     (Mucoromycota), resolved cross-hits (winner and loser identity and
     coverage) in a small table;
   - gene table: gene, role, coordinates, strand, exons, identity, coverage,
     method, model status (with exonerate/miniprot disagreement shown),
     reference record, e-value, bitscore; missing and unsearchable genes as
     rows marked "not found" / "not searchable".
4. **Search summary.** Routing mode and what it means, genetic code, families
   attempted, `not_detected` with reason and best fraction found.
5. **Withheld loci** (collapsed by default): counts per gate, then a table
   with location, genes and the reason, with a one-line explanation of each
   gate (`modelled_gene_bar`, `below_fraction_floor`, `mat_gene_gate`,
   `paralog_class`, `flank_carried`).
6. **How to read this report.** Glossary of confidence tiers, locus classes,
   model statuses, routing modes; the limits that apply to every report
   (2026-10-04 spec, "Risks and limits").
7. **Provenance.** The `run` block in full, the output files and how to cite.

### 1.5 Gene-order figure
- Inline SVG, written by MATPredict (no matplotlib or pyGenomeViz), so it is
  vector in the PDF, small, and needs no extra dependency in `detect`.
- One horizontal track at genome scale over the locus +/- a margin: genes as
  arrows in their strand direction, exons as filled blocks joined by intron
  lines; the label is the gene name; fill encodes the role (core MAT,
  conserved flank, variable flank), and the outline style encodes the model
  status (solid = polished and agreeing, dashed = the two tools disagree,
  dotted = hit only, not modelled). Colour is never the only code: role is
  also in the label row and the legend.
- A second thin track shows the alternate model (miniprot) where it differs.
- The contig drawn as a line with its ends marked when they fall inside the
  window ("contig end 41 kb" otherwise); assembly gaps as hatched blocks.
- Scale bar in kb; coordinates in the caption.
- Width set by the viewBox, so the figure scales to the screen and to the
  printed page.

### 1.6 Visual and print design
- Light, document-style page (a report, not an app); a dark scheme for the
  screen under `prefers-color-scheme: dark`, always printed light.
- System font stack, tabular numerals in tables, a maximum line length near
  75 characters for prose, tables full-width.
- Colour tokens as CSS custom properties; the role palette chosen for colour
  vision deficiency and checked for contrast (WCAG AA for text).
- Print: `@page` size A4/Letter-neutral margins, page header (sample) and
  footer (page numbers) where the engine supports margin boxes; no break
  inside a locus card, figure or table row; long tables repeat their header;
  collapsed sections are opened before printing (a small `beforeprint`
  handler; without JS everything prints open because the print stylesheet
  does not hide content); links print with their target where useful.
- No horizontal scroll at phone width: wide tables scroll inside their own
  container on screen and wrap in print.

### 1.7 PDF
1. In a browser: "Print, Save as PDF" (the print stylesheet does the layout).
   A "Save as PDF" button in the header calls `window.print()`.
2. On the command line, `--pdf`: WeasyPrint when it is installed (pure
   Python, about 50 MB with Pango; an optional pixi feature `report-pdf`, and
   in the container), otherwise a headless Chrome or Chromium on PATH
   (`--headless --print-to-pdf --no-pdf-header-footer`). Neither found: a
   clear error naming both.
3. The CSS uses only features both engines support (flexbox, no CSS grid
   subtleties, no JS-dependent layout), and the test suite renders one PDF
   with whichever engine is present.

### 1.8 Implementation
- New package `src/MATPredict/report/`: `model.py` (YAML + GFF3 into
  dataclasses with defaults for old reports), `svg.py` (figure),
  `html.py` (page, escaping every string), `pdf.py` (engine selection),
  `cli.py`.
- No new run-time dependency for HTML (standard library `html` and string
  templates). Jinja2 is not used, to keep `detect`'s environment unchanged.
- Tests (`tests/report/`): fixture reports covering one call, nothing called,
  `not_searched`, two idiomorphs, a Mucoromycota classifier call, a merged
  Basidiomycota A/B call with subloci, a PR call with a receptor array, an old
  report without `run`; checks: renders without error, every gene of the
  YAML appears, escaped input (a contig named `<script>`), parses as HTML,
  one PDF render when an engine is present.

### 1.9 Review
Done 2026-10-07: two rounds, findings and changes in `analysis/2026-10-07_report-design-review.md`.
Before merge, a web- and data-design review of rendered examples (screen at
phone and desktop width, and the PDF): hierarchy and first-screen answer,
figure legibility, table density, colour and contrast, print pagination. The
findings and what was changed are recorded in `analysis/`.

## Part 2: reads intake (reads plus organism to a MAT call)

### 2.1 Inputs
- Reads: one or more SRA run accessions (SRR/ERR/DRR), or uploaded FASTQ
  (plain, gz, zst; single or paired).
- Organism: a taxid, or a scientific name resolved against the offline
  taxonomy table (exact match first, then a case-insensitive match; an
  ambiguous name is an error that lists the candidates).

### 2.2 Getting reads from SRA
- Ask the ENA portal API (`filereport?result=read_run&fields=fastq_ftp,read_count,base_count,library_layout,instrument_platform`)
  for the run's FASTQ URLs. ENA mirrors SRA runs as gzipped FASTQ over HTTPS,
  and a stream can be cut after N reads, so a typing job downloads about
  1-2 GB instead of the whole run.
- Fallback when ENA has no FASTQ for the run: `fasterq-dump` (sra-tools, an
  optional pixi feature) with a read limit.
- Refuse, with a reason: long-read-only runs for the k-mer path (to test:
  ONT/PacBio error rates break exact 31-mers), amplicon or RNA-seq libraries
  (`library_strategy` other than WGS), runs whose `scientific_name` disagrees
  with the given organism at genus level (a warning, not a refusal, because
  SRA names are often wrong).

### 2.3 Panel builder (new: `matpredict reads-panel`)
From a taxid, build the `reads-type` panel from `db/`:
1. Accepted records of the same species with at least two idiomorphs of one
   family; each idiomorph sequence is the region between the inner flank
   ends from `locus.gbk`, plus the flank ends that become the shared control.
2. If the species has none: same genus, marked `panel_level: genus`.
3. The panel is cached per (taxid, database hash) and its record ids and
   level are written into the reads report.

### 2.4 Route
| Situation | Route | Expected cost |
|---|---|---|
| Species-level panel exists (Fola, A. fumigatus, Rhizopus, Clavispora, ...) | `reads-type` (k-mers) | minutes, under 1 GB |
| Genus-level panel only | `reads-type`, report marked lower confidence; `none` triggers the assembly route | minutes |
| No panel; Basidiomycota HD/PR (multiallelic, divergent); or `none`/`low_depth` that the user wants pursued | assemble, then `detect` | SPAdes or MEGAHIT on about 40x: 16-32 GB, 1-3 h (to measure) |
| Later (not built): targeted assembly | DIAMOND blastx of reads against MAT and flank proteins, assemble only recruited reads, `detect` | to measure; much cheaper than a full assembly |

Both routes can run for one sample when the user asks; the report then shows
the two calls side by side, because they fail differently: assembly loses
short or heterozygous loci (A. fumigatus hybrids), k-mers lose divergent alleles.

### 2.5 Reads report
`reads-type` gains `--json` (one record per sample: call, breadth, depth and
expected breadth per idiomorph, shared depth, reads used, flags, panel record
ids and level). `report reads` renders it with the same page frame as Part 1:
the result sentence, a breadth-versus-expected chart per idiomorph, depth
ratio, and the panel provenance. "both" is always shown with the depth ratio
and the possible causes (heterokaryon, diploid, mixed or contaminated library).

### 2.6 Before a public service offers reads typing
1. The simulated mixes (50:50, 80:20, 95:5) the curator asked for
   (2026-10-05) are run and the `both` thresholds set from them.
2. One independent test set (not the Fola cohort the thresholds were set on),
   for example the LCG strains with assembly calls.
3. The panel builder reproduces the Fola panel from the database records
   (AB011379.2, AB011378.1) and gives the same 148 calls.

## Phased plan

| Phase | Content | Depends on |
|---|---|---|
| 1 | `run` provenance block in `detect`; `report genome` HTML + SVG figure; tests | - |
| 2 | `--pdf` (WeasyPrint / headless Chrome); design review of rendered examples; fixes | 1 |
| 3 | `reads-type --json`, `report reads`; `reads-panel` builder; ENA streaming input | 1 |
| 4 | Simulated-mix calibration and an independent test of `reads-type` | 3 |
| 5 | Wire both into the chosen host (Nextflow first; then Galaxy or the NRP app, 2026-10-04 spec B2) | 2, 3 |

## Questions for the reviewer
1. ~~`detect --html` default~~: decided 2026-10-07, on by default everywhere, opt out.
2. ~~PDF engine~~: decided 2026-10-07, WeasyPrint in the pixi environment (10.3 MB download against about 300 MB
   for Chromium); revisit if its output proves insufficient.
3. ~~Withheld loci~~: decided 2026-10-07, always-present section, table collapsed on screen and expanded in
   print.
4. Genus-level panels: offer them at all in the service, or species-level only
   until a divergence test exists?
5. Assembly route in the service: offer it (hours of compute per sample) or
   leave it to Nextflow/HPC users?
