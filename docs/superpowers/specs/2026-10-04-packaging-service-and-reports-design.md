# MATPredict: packaging, on-demand service, and reports

Date: 2026-10-04
Status: draft for review. Nothing here is implemented unless it says so.

## Context

MATPredict v0.6.x detects MAT loci in a genome FASTA and writes
`detection_report.yaml`, `detected_loci.gff3`, `detected_loci.fasta` and an
optional evidence log. It now runs at scale: about 23,600 BFD genomes on
2026-10-03/04 (Ascomycota 19,415, Basidiomycota 3,275, Mucoromycotina 974).

The curator asked for four things (2026-10-04):

1. Run from a container. Can the database be embedded? Is a conda package
   needed in addition to pixi?
2. Run on demand: a user uploads a genome (or a folder of genomes) and gets a
   report back, for example on a Kubernetes pod through nrp.ai.
3. A richer per-locus report: gene-order plot with evidence tracks, optional
   ab initio gene prediction in the locus window, a species-context panel
   (tree of close relatives with a MAT ideogram per tip), and gene trees of the
   key MAT genes.
4. Run statistics: a per-genome report and an aggregate report for a
   collection.

This document sorts these into four workstreams, states what is measured and
what is not, gives a recommendation for each, and lists the questions a
reviewer should settle.

## Measured facts this design depends on

| Fact | Value | Source |
|---|---|---|
| Curated database, apparent size | 3.8 MB in 619 files (`du` shows 122 MB: block overhead) | `find db -type f -printf %s` |
| `.faa`, `.hmm`, `.yml`/`.yaml`, `.tsv` files (`.gbk` 2.3 MB, `.gff3` 0.06 MB not counted) | 1.5 MB; which files `detect` reads at run time is not checked | same |
| Per-genome wall time, one core | median 101 s (<50 Mb) to 716 s (>500 Mb, Basidiomycota); Ascomycota median 162 s, max 3,095 s; Mucoromycotina median 183 s | campaign runs 2026-10-03 |
| Memory | 32 parallel genomes peaked at 16.7 GB RSS (about 0.5 GB each) | `sacct` MaxRSS, pilot jobs |
| Network use | NCBI E-utilities for taxid lineage and genetic code; none with `--phylum` and `--genetic-code` | `db/taxonomy.py`, `db/ncbi_client.py` |
| NCBI response cache | 1.7 GB on /bigdata after the campaigns | `.matpredict_cache` |
| Container | Dockerfile with pixi build stage and Ubuntu 24.04 runtime; embeds `src/` and `db/`; CI build and smoke test pass (PR #12, not yet merged) | `.github/workflows/docker.yml` |
| Image size | not measured | to do |
| Ab initio predictor in the environment | Augustus is a dependency but no code calls it; Helixer left out for packaging reasons (TensorFlow pins, Python <= 3.11 wheels) | `pixi.toml` comments |
| Drawing code | pyGenomeViz locus and synteny drawing for curated records (`curate-db draw-locus`, `draw-synteny`); campaign figures are scripts in `results/` | `src/MATPredict/db/draw.py` |

## Workstream A: packaging and distribution

### A1. Container
- It already works (PR #12): the image carries the pixi environment, the
  package and `db/`, with `MATPREDICT_DB_ROOT=/app/db`. It is published to GHCR
  on a `v*` tag.
- Recommendation: keep the database embedded. At 3.8 MB it costs nothing, and
  the classifier HMMs and family rules are tied to the code version (a
  database from another version can change calls). One image = one code +
  database version = reproducible reports.
- Gaps to close:
  - Offline taxonomy. A service should not depend on NCBI E-utilities per job.
    `detect` needs three things per taxid: the lineage (taxids up to the root),
    the phylum name, and the nuclear genetic code. All three are in
    `nodes.dmp` (parent, rank, genetic code id) and `names.dmp` (scientific
    name). Measured on the HPCC copy (`/srv/projects/db/taxonomy`, May 2026):

    | Option | Size |
    |---|---|
    | Full taxdump (`taxdump.tar.gz`) | 74 MB compressed; nodes + names 504 MB unpacked |
    | Slim table, all taxa: taxid, parent, rank, gencode + scientific names (2.83 M taxa) | 25 MB zstd |
    | Slim table, Fungi only (taxid 4751 subtree, 219,293 taxa) | 2.2 MB zstd |

    Recommendation: embed the slim all-taxa table (25 MB) plus `merged.dmp`
    (1.9 MB raw; old taxids that NCBI merged) in the image, and add a local
    taxonomy backend that reads it, so a job makes no network call. All taxa,
    not Fungi only, so a wrong or non-fungal taxid gets a clear answer rather
    than "unknown". Rules:
    - Every report states the taxonomy snapshot date.
    - A taxid newer than the snapshot: use E-utilities only if network use is
      allowed (CLI default); in the service, report it as unresolved and ask
      for `--phylum` / `--genetic-code`.
    - A maintainer command rebuilds the table from a new taxdump at release.
    Alternatives: the full taxdump with `taxonkit` (already a dependency;
    larger image), or require `--phylum` and `--genetic-code` from every user.
    Compressed reading (measured 2026-10-04):
    - `taxonkit` v0.20.0 reads gzip and zstd content (it detects the format),
      but only under the plain names (`nodes.dmp`, ...); `nodes.dmp.gz` or
      `.zst` names give "taxonomy data not found". A symlink `nodes.dmp ->
      nodes.dmp.zst` works. Compressed nodes + names + merged + delnodes: 67.4 MB
      gzip, 65.8 MB zstd. Each call loads the whole taxonomy: about 2 s and
      440 MB RAM, so batch taxids into one call.
    - Python 3.14 (our environment) has `gzip` and `compression.zstd` in the
      standard library, so a MATPredict reader of the slim `.zst` table needs no
      new dependency.
  - Measure the image size and cold-start time.
  - Apptainer/Singularity: build from the GHCR image for HPC use (README has the
    command; test on UCR HPCC).

### A2. Conda package (bioconda)
- pixi installs from conda channels, so a bioconda recipe does not conflict with
  pixi; pixi can consume it.
- Reasons to add one:
  - BioContainers builds a Docker and Singularity image automatically from every
    bioconda release.
  - Galaxy tool wrappers and nf-core modules expect a bioconda package plus a
    BioContainers image.
  - `conda install matpredict` for users who do not use pixi.
- Cost: a recipe (`noarch: python`, run deps = the binaries in `pixi.toml`),
  and bioconda review. The database ships inside the package (small, version-
  locked).
- Recommendation: yes, after v0.6.1, as the base for workstream B options that
  use Galaxy or Nextflow.

### A3. The "MAT atlas" (new data product)
Workstreams C and D need reference context: which close relatives carry which
loci, their gene sequences, and a tree. The campaign runs already contain this
for about 23,600 genomes. Package it as a versioned data release, separate from
the code image:
- per genome: taxonomy, calls (family, idiomorph, confidence, coordinates,
  locus size, flank status);
- per locus: protein sequences of polished core and flank genes; the locus
  region sequence (locus +/- 10 kb) for drawing;
- per gene family (sexP/sexM, MAT1-1-1, MAT1-2-1, HD1, HD2, STE3, flanks):
  reference alignment and tree;
- a species tree or taxonomy cladogram over atlas genomes.
- Size: not measured. Rough estimate: proteins well under 100 MB; region
  sequences about 28,000 loci x 20 kb, a few hundred MB compressed. Measure
  before choosing where to host it (Zenodo, GHCR OCI artifact, or an NRP S3
  bucket).
- Held-out sets (LCG, Jena) stay out of the atlas used for any build or score.
  Whether they may appear as context tips is a curator decision.

## Workstream B: on-demand service

### B1. What a job needs
One genome = one single-core process, about 0.5-1 GB RAM, 2-12 minutes
typical, up to about 1 h for >500 Mb genomes. Inputs: a FASTA (plain or gz),
optional taxid or phylum, optional genetic code. Outputs: a folder of about
1-10 MB per genome. This fits a Kubernetes Job (one pod per genome or per
small batch) well.

### B2. Options

| Option | Upload/UI | Compute | Effort | Notes |
|---|---|---|---|---|
| B-a. Small web app on NRP | FastAPI (or similar) + upload page | one k8s Job per genome in an NRP namespace; results to S3; report page | medium-high | We own auth, quotas, storage clean-up, abuse limits. NRP specifics (job API for a service account, public ingress, S3, usage policy for a public service) are **not verified** here and must be checked with NRP. |
| B-b. Galaxy tool | Galaxy upload and history; collections for folders | Galaxy job runner (usegalaxy.* or a Galaxy on NRP; Galaxy has a k8s runner) | low-medium once bioconda exists | Users already know Galaxy; collections give the aggregate mode; HTML report shows in the history. Needs A2. |
| B-c. Nextflow pipeline (nf-core style) | none (CLI, or Seqera Platform launchpad) | any executor incl. k8s, SLURM, AWS Batch | medium | Best for "folder of genomes" and for HPC users; Seqera Platform can give a launch form. |
| B-d. CLI + container only | none | user's machine or HPC | done (A1) | Baseline. |

Recommendation:
1. B-d now (exists).
2. B-c next: a small Nextflow pipeline (scatter genomes, `detect`, per-genome
   report, aggregate report). It is also the engine behind B-a or B-b.
3. Then pick B-b or B-a for the upload interface. B-b gives the most for the
   least work if Galaxy is acceptable. B-a only if a branded upload page on NRP
   is a requirement; first confirm with NRP that a public on-demand service is
   allowed and how jobs are launched.

### B3. Service-level rules (any option)
- Pin one image version per deployment; print it in every report.
- No NCBI calls per job (A1 offline taxonomy). No e-mail in the image.
- Input limits: maximum genome size and count per submission; reject non-FASTA.
- Results expire after N days; uploaded genomes are not kept.

## Workstream C: per-genome report

One self-contained HTML file per genome (plus the existing YAML/GFF3/FASTA),
with these panels. Each panel lists the tools and inputs it used.

### C1. Locus summary table (no new compute)
Per called locus: family, idiomorph and margin, confidence, locus class,
contig and coordinates, locus size (inner flank ends, when both flanks are
polished), genes found (core/flank/optional), gene-model status per gene
(`polished_agree`/`disagree`/`single`/`unpolished`), assembly context (contig
length, distance to contig end, gaps near the locus), flags (two idiomorphs,
homothallic candidate, flank-carried, unverified). Also: routing mode, genetic
code, families searched, and withheld loci with their reasons.

### C2. Gene-order plot with evidence tracks (low compute)
For each locus, the locus +/- a window, drawn with pyGenomeViz:
- track 1: called genes (core, flank, optional), coloured by role;
- track 2: tblastn hits (raw localisation evidence);
- track 3: exonerate models; track 4: miniprot models (shows where they
  disagree);
- track 5: the input annotation, if the user gives a GFF3;
- track 6 (optional, C3): ab initio genes.
All of this except C3 is already in the evidence the pipeline holds.

### C3. Ab initio gene prediction in the locus window (optional)
Purpose: show every gene in the window, not only the genes we search for.

| Tool | In env | Speed on a 20-60 kb window | Notes |
|---|---|---|---|
| Augustus | yes | seconds | needs a species model; no Mucoromycota-specific model is known to us; choose the nearest available model per phylum (to test) |
| Helixer | no | not measured (GPU helps; CPU possible on small windows) | lineage models (fungi); TensorFlow stack does not fit the main environment; run as a separate container step |
| miniprot with a broad protein set (for example SwissProt fungi) | yes (miniprot) | seconds-minutes | homology evidence for non-MAT genes |

Recommendation: Augustus in v1 (already installed; measure accuracy on curated
loci, where we know the true models). Helixer as an optional containerised step
in the Nextflow pipeline (B-c), not in the core package.

### C4. Species-context ideogram (needs A3)
Left: a tree of the target plus 5-10 close relatives from the atlas. Right: one
row per tip, a simple gene-order ideogram of each MAT locus (genes as arrows,
coloured by role; idiomorph label), aligned on the core gene.
- Choosing relatives: nearest atlas genomes by NCBI taxonomy (same genus, then
  family, then order), prefer curated-record genomes and reference assemblies,
  at most one per species, both idiomorphs if present.
- Tree: v1 = taxonomy cladogram (no branch lengths; cheap, no compute). v2 =
  phylogenomic placement (BUSCO markers + a reference tree; tens of minutes per
  genome) only if the reviewer judges it worth the cost.
- This is the "Mucoromycotina synteny figure" from the campaign, made
  per-target and automatic.

### C5. Gene trees of key MAT genes (needs A3)
For each core gene the target carries (sexP/sexM together with a non-MAT HMG
outgroup; MAT1-1-1/MAT1-2-1; HD1/HD2; STE3) and optionally one flank gene
(for example tptA):
- add the target's protein(s) to the atlas reference alignment for that family
  (MAFFT `--add`), then a fast tree (FastTree or IQ-TREE `-fast`) on the target
  plus the relatives from C4, or a placement (EPA-ng) on the full reference tree;
- mark the target tips; show support values.
- Expected run time: seconds to a few minutes for 10-50 sequences (not measured).
- Caution from the campaign: on the full-length Mucoromycotina tree the sexP and
  sexM clades were clean (0 tips of the other type) but the non-MAT outgroup was
  not one clade, so rooting needs care. Use midpoint or a fixed outgroup set
  chosen per family and say which in the figure legend.

## Workstream D: aggregate report (collection of genomes)

One HTML + TSV set for a run of N genomes:
- run summary: genomes in, reports, failures (with the last error line), time
  per genome, image/code version;
- call table: one row per locus (the `loci.tsv` the campaign scripts write);
  one row per genome (`genomes.tsv`);
- call rate by taxon rank (order/class), routing mode and family; with and
  without unverified-only calls (for example Basidiomycota PR);
- idiomorph balance per taxon (MAT1-1 vs MAT1-2, Plus vs Minus), two-idiomorph
  and homothallic-candidate lists;
- locus-size distribution per family/genus and idiomorph;
- reasons for no call (bar withheld, nothing localised, assembly gap, not
  searched);
- optional: the C4 ideogram for the whole collection (one row per genome) when
  N is small enough to draw (for example <= 100).
- A MultiQC custom-content file, so the aggregate drops into existing QC
  reports.

The campaign scripts (`analyze_full.py`, `summarise.py`, `analyze_asco.py`,
`collect.py`) already compute most of this; D is mainly moving that logic into
`matpredict report aggregate` with tests.

## Proposed command surface (for review)

```
matpredict detect ...                       # unchanged
matpredict report genome  --run DIR [--atlas ATLAS] [--annotation GFF3] [--abinitio augustus|none]
matpredict report aggregate --runs DIR... --out DIR [--taxonomy samples.tsv]
matpredict atlas build    --runs DIR... --out ATLAS   # maintainers only
```

## Phased plan

| Phase | Content | Depends on |
|---|---|---|
| 0 | Merge #11-#13, tag v0.6.1; measure image size; offline taxonomy (taxonkit + taxdump) | - |
| 1 | `report genome`: C1 + C2 (no new compute); `report aggregate`: D from the campaign scripts | 0 |
| 2 | Nextflow pipeline (B-c): scatter, detect, report genome, report aggregate; containers per step | 1 |
| 3 | Atlas v1 (A3) from the 2026-10 campaigns; C4 with taxonomy cladogram; C5 with MAFFT --add + FastTree | 1 |
| 4 | C3 Augustus in the window (with an accuracy check on curated loci); Helixer as an optional pipeline step | 2 |
| 5 | Bioconda recipe (A2); then Galaxy tool (B-b) or NRP web front end (B-a) | 0, 2 |

## Questions for the reviewer

1. Embed the database in the image and package (recommended), or ship it
   separately?
2. Offline taxonomy: embed the slim all-taxa table (25 MB, recommended), the full
   taxdump with taxonkit, or require phylum and genetic code from the user?
3. Upload interface: Galaxy, a custom NRP web app, or Nextflow + Seqera? Is a
   public on-demand service allowed on NRP, and under which account?
4. Ab initio: is Augustus with the nearest species model good enough for the
   evidence track, or is Helixer required from the start?
5. Species context: taxonomy cladogram (cheap) or phylogenomic placement
   (accurate, costly)? How many relatives per panel?
6. Gene trees: which genes per lineage (core genes only, or one flank gene)?
   Placement on a fixed reference tree, or a fresh small tree per report?
7. May held-out LCG/Jena genomes appear in the atlas as context tips (never in
   builds or scores)?
8. Should unverified-only calls (for example Basidiomycota PR via strict CAAX)
   be shown in reports by default, or behind a switch?
9. Where to host the atlas (Zenodo, GHCR OCI artifact, NRP S3)?

## Risks and limits to state in every report

- Lineages with no curated record are searched against the whole phylum; family
  labels there are not reliable (for example MATsc and MTL in one Alaninales
  genome), and a missing call is a reference gap, not evidence of absence.
- Gene models: exonerate and miniprot disagree on boundaries often in
  under-curated lineages (21 genes in 15 Alaninales genomes).
- The tool is not a mating-type assay: two idiomorphs in one assembly can be a
  homothallic species, a mixed culture or a merged diploid assembly.

## Side note (2026-10-05): reads-based typing and population locus discovery
Two further lines of work are specified separately:
`docs/superpowers/specs/2026-10-05-reads-typing-and-population-locus-discovery-design.md`
(idiomorph typing from unassembled reads against a species' known idiomorphs;
and discovery of MAT-like loci from presence/absence and depth patterns across a
population of strains, with Batrachochytrium and Xylariales as targets). They
would be developed on their own branch and are not part of the phases above.
