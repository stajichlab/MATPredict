# MATPredict

Find and type the mating-type (MAT) loci of fungi in genome assemblies.

MATPredict searches a genome assembly with a curated reference database of MAT
loci. It reports where each locus is, which genes it carries, and which
idiomorph or mating type it is (for example Plus/Minus, MAT1-1/MAT1-2, HD/PR
alleles). Every reference record traces to a publication or a validated
deposit, and every call reports the evidence behind it.

Current release: [`v0.6.1`](https://github.com/stajichlab/MATPredict/releases/tag/v0.6.1).
Changes: [CHANGELOG.md](CHANGELOG.md).

## Contents

- [Project goals](#project-goals)
- [Workflow](#workflow)
- [Supported lineages](#supported-lineages)
- [Curation process](#curation-process)
- [Repository layout](#repository-layout)
- [Quickstart](#quickstart)
- [Usage](#usage)
- [Validation and known limits](#validation-and-known-limits)
- [Authors](#authors)
- [Citation](#citation)
- [License](#license)

## Project goals

1. **Curated reference loci.** Develop a curated set of model MAT loci for
   different taxonomic groups, one record per strain and idiomorph, with
   provenance for every gene.
2. **Detection.** Use translated searches to find the key MAT genes and their
   flanking genes in new genome assemblies, annotated or not.
3. **Typing.** Classify each locus by the mating-type naming scheme of its group
   (for example Plus/Minus, MAT1-1/MAT1-2, HD/PR alleles), and report when the
   evidence is not enough to decide.
4. **Standard reports.** Extract the locus region and its genes, and report them
   in a standard format with figures and tables. Report uncalled genomes with a
   reason, and two idiomorphs in one genome as a flag for review.
5. **Synteny.** Visualize the synteny of MAT regions across sets of species.
6. **Scale.** Run on thousands of genomes (for example the full BFD fungal
   genome collection) on a SLURM cluster.

## Workflow

```
                 genome assembly (FASTA)  [+ proteins]  [+ NCBI taxid]
                                  │
                                  ▼
 ┌──────────────────────────────────────────────────────────────────────┐
 │ 1. ROUTE      taxid lineage ──► curated families in db/<Phylum>/      │
 │               order.yml. No match ──► phylum fallback, or             │
 │               "not_searched". --phylum forces one phylum.             │
 ├──────────────────────────────────────────────────────────────────────┤
 │ 2. SEARCH     tblastn (genome) / diamond (proteome) with the          │
 │               reference proteins: core MAT genes + flanking genes     │
 ├──────────────────────────────────────────────────────────────────────┤
 │ 3. CLUSTER    hits ──► candidate loci per contig (lineage gap rules)  │
 ├──────────────────────────────────────────────────────────────────────┤
 │ 4. POLISH     exonerate / miniprot gene models on a sliced window     │
 │               (at most 6 loci per family)                             │
 ├──────────────────────────────────────────────────────────────────────┤
 │ 5. GATE       modelled-gene bar, fraction floor, MAT-gene gate,       │
 │               flank-carried rule, paralog class, CAAX check           │
 ├──────────────────────────────────────────────────────────────────────┤
 │ 6. TYPE       idiomorph from gene content, or profile-HMM classifier  │
 │               (Mucoromycota sexP/sexM); margin and confidence tier    │
 └──────────────────────────────────────────────────────────────────────┘
                                  │
                                  ▼
     detection_report.yaml   detected_loci.gff3   [detected_loci.fasta]
     (calls, evidence, withheld loci, uncalled reason, two-idiomorph flag)
```

The reference database feeds steps 1, 2 and 6:

```
 literature / GenBank deposit
        │  propose (PMID/DOI + source text required)
        ▼
 db/candidates/ ──► automated validation ──► curator sign-off ──► db/<Phylum>/<Order>/
                    (accession resolves,       (accept or          locus.gbk, locus.gff3,
                     translation matches,       reject with         proteins.faa,
                     taxonomy current)          a reason)           metadata.yaml
                                                                         │
                                         classifiers, regression check ◄─┘
```

## Supported lineages

The curated families are defined in each phylum's `order.yml`.

| Phylum | Families (`locus_name`) | Records | Definition |
|---|---|---|---|
| Ascomycota | MAT, PM, mat1, mat2, mat3, MATsc, MATyl, MTL, MATtub | 57 | [db/Ascomycota/order.yml](db/Ascomycota/order.yml) |
| Basidiomycota | HD, PR, MAT, Aalpha, Abeta, Balpha, Bbeta, aLocus, bLocus, rustHD, redPR, redHD, wallMAT | 73 | [db/Basidiomycota/order.yml](db/Basidiomycota/order.yml) |
| Mucoromycota | MAT (sexP / sexM, with an idiomorph classifier) | 20 | [db/Mucoromycota/order.yml](db/Mucoromycota/order.yml) |

Other phyla (for example Mortierellomycota, Kickxellomycota, Chytridiomycota)
have no curated family. By default their genomes are reported as
`not_searched`. Use `--exhaustive` to search them with every family.

## Curation process

Each reference record is a folder `db/<Phylum>/<Order>/<record_id>/` with four
files. The record id is `<taxid>_<strain>_<locus>_<idiomorph>`.

| File | Content |
|---|---|
| `metadata.yaml` | Taxonomy, organism, strain, mating type, locus coordinates, genes, evidence (tier, PMID/DOI, source text), validation result, curation history. Schema: [db/_schema/metadata.schema.yaml](db/_schema/metadata.schema.yaml) |
| `locus.gbk` | The locus as GenBank |
| `locus.gff3` | Gene features of the locus |
| `proteins.faa` | Reference proteins, tagged with gene name and role (`core_MAT`, `flanking_conserved`, `flanking_variable`) |

Rules:

1. **Evidence.** Tier 1: a published and experimentally validated locus. Tier 2:
   genome-derived coordinates, admitted only with extra checks (domain evidence,
   alignment to curated orthologs, or a gene tree) and a record of how the
   locus was derived. A source never comes from inference alone.
2. **Proposal.** A candidate record goes to [db/candidates/](db/candidates/) with
   a resolvable PMID/DOI or accession and the exact source text.
3. **Automated validation.** The accession must resolve, the translation must
   match the cited protein, and the taxonomy must be current.
4. **Curator sign-off.** Only the curator moves a record to `accepted`. Rejected
   candidates are kept with a reason.
5. **Regression check.** A change to a record, classifier or rule runs the
   regression panel before sign-off. See [docs/regression-check.md](docs/regression-check.md).

Related files:

| Path | Purpose |
|---|---|
| [db/_schema/](db/_schema/) | JSON schemas for `metadata.yaml` and `order.yml`; DuckDB schema |
| [db/Mucoromycota/classifiers/](db/Mucoromycota/classifiers/) | Idiomorph classifier: HMMs, build manifest, paralog negatives |
| [db/suppress.txt](db/suppress.txt) | Genomes that must not be searched (for example amplicon-only deposits) |
| [db/taxon_overrides.tsv](db/taxon_overrides.tsv) | Genomes whose deposited name is likely wrong (used in held-out scoring only) |
| [db/assembly_zygosity.yml](db/assembly_zygosity.yml) | Lineages where an assembly can collapse a heterozygous locus |
| [ANNOTATION_ERRORS_FIXED_REPORT.md](ANNOTATION_ERRORS_FIXED_REPORT.md) | Errors found in deposits, annotations and literature, and what was done |
| [analysis/](analysis/) | One report per study, the curator's [decisions](analysis/decisions.md), [open questions](analysis/open-questions.md) and an [index](analysis/INDEX.md) |
| [docs/publication-notable-findings/](docs/publication-notable-findings/) | Findings with citeable evidence (new loci, homothallism candidates) |
| [docs/holdout-benchmark.md](docs/holdout-benchmark.md) | How held-out recall is measured |

## Repository layout

| Path | Content |
|---|---|
| [src/MATPredict/](src/MATPredict/) | The Python package and the `matpredict` CLI |
| [src/MATPredict/detect/](src/MATPredict/detect/) | Detection pipeline (routing, search, polish, gates, classifier, reports) |
| [src/MATPredict/db/](src/MATPredict/db/) | Database curation, validation, taxonomy, drawing |
| [db/](db/) | The curated reference database |
| [scripts/](scripts/) | SLURM runners, regression check, classifier build, record tools |
| [tests/](tests/) | pytest suite |
| [testset/](testset/) | Ground-truth coordinates and the regression panel list |
| [analysis/](analysis/) | Study reports and the decision log |
| [docs/](docs/) | Design notes, handoffs, benchmark and regression documentation |

## Quickstart

### Install

MATPredict uses [pixi](https://pixi.sh) for its environment. The environment
holds Python, BLAST+, DIAMOND, miniprot, exonerate, Augustus, pyhmmer, MAFFT and
the drawing tools. Linux (x86-64) only.

```bash
git clone https://github.com/stajichlab/MATPredict.git
cd MATPredict
pixi install
pixi run matpredict --help
```

### Other ways to install

**Docker.** Each release publishes an image to the GitHub Container Registry.
The image holds the environment, the package, the curated database and an
NCBI taxonomy snapshot (see [Offline taxonomy](#offline-taxonomy)), so a run
needs no network.

```bash
docker pull ghcr.io/stajichlab/matpredict:latest     # or a release version, :<version>
docker run --rm -v "$PWD":/data \
  ghcr.io/stajichlab/matpredict:latest \
  detect --genome /data/genome.fna --taxid 4837 --out-dir /data/out
```

On an HPC system without Docker, use Apptainer/Singularity:

```bash
apptainer pull matpredict.sif docker://ghcr.io/stajichlab/matpredict:latest
apptainer run matpredict.sif detect --genome genome.fna --taxid 4837 --out-dir out
```

**Conda / mamba.** [environment.yml](environment.yml) lists the same packages
as `pixi.toml` (exported with `pixi workspace export conda-environment`). Run
from the repository root, because it installs MATPredict from the checkout:

```bash
mamba env create -f environment.yml
mamba activate matpredict
matpredict --help
```

The exact versions used for releases are in [pixi.lock](pixi.lock).

### Run one genome

This example uses *Phycomyces blakesleeanus* NRRL 1555 (NCBI taxid 4837).

```bash
# 1. Download the assembly from NCBI
curl -L -o phybl.zip \
  "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/GCF_001638985.1/download?include_annotation_type=GENOME_FASTA"
unzip -q phybl.zip

# 2. Detect MAT loci (run from the repository root, so db/ is found).
#    Set your e-mail for NCBI first (see Configuration).
export MATPREDICT_NCBI_EMAIL=you@example.org
pixi run matpredict detect \
  --genome ncbi_dataset/data/GCF_001638985.1/GCF_001638985.1_Phybl2_genomic.fna \
  --taxid 4837 \
  --out-dir phybl_out
```

The run takes about 1.5 minutes on one CPU. `phybl_out/detection_report.yaml`
reports one locus:

```
Mucoromycota:MAT  Minus  high  mat_locus  NW_017265134.1:3978283-3991589
genes: tptA, sexM, rnhA, algA
```

This genome is the source of a curated record
(`4837_nrrl1555_MAT_Minus`), so the example shows that the tool works. It is
not an independent test.

## Usage

### `matpredict detect`

```
matpredict detect --genome GENOME [--proteins PROTEINS] [--taxid TAXID]
                  --out-dir OUT_DIR [options]
```

| Option | Description |
|---|---|
| `--genome` | Genome assembly, FASTA (uncompressed) |
| `--proteins` | Optional proteome FASTA. Enables the DIAMOND fast path |
| `--taxid` | NCBI taxid. Selects the families for the lineage and the genetic code |
| `--phylum {Ascomycota,Basidiomycota,Mucoromycota}` | Search one phylum's families and skip taxid routing |
| `--exhaustive` | Search every family when no curated family covers the genome |
| `--genetic-code N` / `--genetic-code-map FILE` | Translation table, for one genome or per genome in a batch |
| `--max-polished-clusters-per-family N` | Polish cap (default 6; 0 = no cap) |
| `--min-hits N` | Distinct genes a cluster needs before polishing (default 2) |
| `--min-identity PCT` | Best-hit identity a cluster needs before polishing (default: none) |
| `--no-require-core-role` | Allow polishing of clusters without a core MAT gene hit |
| `--exclude-records IDS` | Withhold curated records (comma-separated) for leave-one-out tests |
| `--evidence-diagnostics FILE` | Write per-hit evidence as JSON lines |
| `--html` / `--no-html` | Write `report.html` (default: on). `MATPREDICT_HTML=0` turns the default off; the batch scripts set it |
| `--pdf` | Also write `report.pdf` (WeasyPrint, in the pixi environment and the image; off by default). A report or PDF error is logged and never fails the run |
| `--sample` / `--organism` | Names shown in the report and recorded in the `run` block |
| `--emit-cds-fasta` | Add CDS features with translations to the GFF3, and write `detected_loci.fasta` (full sequence of every contig with a call; can be large) |

**Outputs** (in `--out-dir`):

| File | Content |
|---|---|
| `detection_report.yaml` | Routing mode, genetic code, detected loci (family, contig, coordinates, idiomorph, confidence, `locus_class`, genes found and missing, classifier scores and margin, gene evidence), withheld loci with the reason, zygosity and two-idiomorph flags, families not detected |
| `detected_loci.gff3` | One gene feature per found or missing gene at each locus |
| `detected_loci.fasta` | Only with `--emit-cds-fasta` |
| `report.html` | Readable report (see `matpredict report genome`); not with `--no-html` |
| `report.pdf` | Only with `--pdf` |

**Routing modes** (`routing_mode` in the report): `lineage` (taxid matched a
family scope), `phylum_fallback`, `explicit_phylum` (`--phylum`),
`exhaustive`, `not_searched`.

### `matpredict report genome`

```
matpredict report genome --run OUT_DIR [--out report.html] [--pdf report.pdf] [--sample NAME]
```

`detect` writes this report itself by default; use this command to re-render a run (for example one made before
the report existed, or with `--pdf`). Writes a self-contained HTML report of one `detect` run (default `OUT_DIR/report.html`): the result in plain words,
a gene-order figure and evidence table per called locus, what was searched, withheld candidates and provenance. Open it
in any browser; "Save as PDF" prints it. `--pdf` writes the PDF directly with WeasyPrint (part of the pixi environment,
`environment.yml` and the Docker image), or, where WeasyPrint is missing, a headless Chrome or Chromium on PATH. `detect --sample` and `--organism` set the names the report
shows.

### `matpredict detect` helper commands

| Command | Purpose |
|---|---|
| `matpredict detect rollout-summary --reports-dir DIR --out FILE` | Combine the reports of a batch into one summary |
| `matpredict detect suppress-filter --list FILE` | Remove suppressed genomes from an `ASMID<TAB>...` list |
| `matpredict detect audit-scope` | Check every record's taxid against its family's `taxonomic_scope` |
| `matpredict detect benchmark` | Leave-one-out sensitivity/specificity benchmark |

### Batch runs on SLURM

| Script | Purpose |
|---|---|
| [scripts/run_clade_panel.slurm](scripts/run_clade_panel.slurm) | Run every BFD genome of one clade (`RANK=ORDER CLADE=Mucorales`). Per-genome timeout `GENOME_TIMEOUT` (default 3600 s) |
| [scripts/submit_binned_panel.py](scripts/submit_binned_panel.py) | Split a panel into genome-size bins with matching time limits |
| [scripts/run_detection_batch.slurm](scripts/run_detection_batch.slurm) | Generic batch runner |
| [scripts/run_regression_panel.sh](scripts/run_regression_panel.sh), [scripts/regression_check.py](scripts/regression_check.py) | Regression check: baseline and candidate on a fixed panel, then a diff |

Run a batch from a frozen git worktree on shared storage, not from a tree you
are editing. Set `SRC` and `MATPREDICT_DB_ROOT` to that worktree, so the code
and the database come from the same commit.

### `matpredict curate-db`

| Command | Purpose |
|---|---|
| `propose --phylum P --record-file F` | Add a candidate record to `db/candidates/` |
| `validate --phylum P --record-id ID` | Run automated validation on a candidate |
| `accept --phylum P --order-or-family O --record-id ID` | Move a validated candidate into `db/<Phylum>/<Order>/` (curator only) |
| `reject --phylum P --record-id ID --reason TEXT` | Reject a candidate and keep the reason |
| `build-gff` / `backfill-gff` | Write `locus.gff3` from the record |
| `draw-locus` / `draw-synteny` | Draw one locus, or synteny across records (pyGenomeViz, clinker) |
| `build-duckdb [--out FILE]` | Build the DuckDB query cache from the record files |
| `release` | Stamp a database release on all records |

### Configuration

MATPredict looks up NCBI taxonomy (for `--taxid` routing and the genetic code)
through NCBI E-utilities. NCBI asks every caller to send an e-mail address.
MATPredict has no default address. Set yours once, in one of two ways:

1. **Config file** (recommended). Create `~/.config/matpredict/config.toml`:

   ```toml
   [ncbi]
   email = "you@example.org"
   # api_key = "your-ncbi-api-key"   # optional: higher request rate
   ```

   MATPredict reads `$MATPREDICT_CONFIG` instead when it is set, or
   `$XDG_CONFIG_HOME/matpredict/config.toml` when `XDG_CONFIG_HOME` is set.

2. **Environment variables**, for example in `~/.bashrc` or a SLURM script:

   ```bash
   export MATPREDICT_NCBI_EMAIL=you@example.org
   export MATPREDICT_NCBI_API_KEY=your-ncbi-api-key   # optional
   ```

An environment variable wins over the config file. With no address set, the
tool still runs, sends no e-mail, and logs one warning.

| Variable | Default | Purpose |
|---|---|---|
| `MATPREDICT_NCBI_EMAIL` | none | E-mail sent to NCBI E-utilities |
| `MATPREDICT_NCBI_API_KEY` | none | NCBI API key |
| `MATPREDICT_CONFIG` | `~/.config/matpredict/config.toml` | Config file path |
| `MATPREDICT_DB_ROOT` | `./db` | Reference database to use |
| `MATPREDICT_CACHE_DIR` | `./.matpredict_cache` | Cache for NCBI responses. The e-mail and API key are not part of the cache key, so users with different settings share one cache |
| `MATPREDICT_TAXONOMY` | none (image: `/app/taxonomy/ncbi_taxonomy.tsv.zst`) | Offline NCBI taxonomy table; used before E-utilities |
| `MATPREDICT_OFFLINE` | off (image: `1`) | `1`: never call NCBI; a taxid not in the table is reported in `routing_error` |

### Offline taxonomy

`detect` needs the lineage, phylum and genetic code of a `--taxid`. A local
table answers these without NCBI E-utilities. Build it from a dated NCBI
archive and point `MATPREDICT_TAXONOMY` at it:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2026-10-01.zip
matpredict curate-db build-taxonomy --taxdump taxdmp_2026-10-01.zip \
    --snapshot 2026-10-01 --out ncbi_taxonomy.tsv.zst
export MATPREDICT_TAXONOMY=$PWD/ncbi_taxonomy.tsv.zst
export MATPREDICT_OFFLINE=1        # optional: never call NCBI
```

All taxa, about 24 MB; loading it takes about 2 s and 80 MB of memory. Each
report records the source in `taxonomy_source` (for example
`local NCBI taxonomy snapshot 2026-10-01`). The Docker image builds this table
at image build time and sets both variables.

### Tests

```bash
pixi run -e test pytest -q
```

## Validation and known limits

Measured results (details in [analysis/INDEX.md](analysis/INDEX.md)):

| Test | Result |
|---|---|
| Zygomycete ground-truth set (23 genomes) | 23/23 loci and idiomorphs |
| Held-out LCG Mucoromycotina genomes (never used in curation) | 536/621 called |
| Held-out Jena Mucoromycotina genomes | 61/64 called |
| Idiomorph typing validation | 0/540 typed wrong |
| Pezizomycotina hold-out, order radius | 7/7 |
| Basidiomycota, 3,270 BFD genomes | Agaricomycotina 83.0%, Ustilaginomycotina 95.6%, Pucciniomycotina 13.2% called |

Known limits (open items in [analysis/open-questions.md](analysis/open-questions.md)):

- Only three phyla have curated families.
- Ascomycota outside the curated orders has had no full-scale run. A genome
  far from every reference often fails the evidence bar, so an uncalled genome
  is weak evidence that the genome has no MAT locus.
- Fragmented assemblies can split a locus or lose its flanking genes.
- In Mucoromycota, *Absidia* flanking genes are often off-scaffold, and the
  gene order of the Lichtheimiaceae is not yet known.
- Basidiomycota pheromone-receptor calls that depend only on the CAAX motif are
  marked unverified. Each pheromone-receptor call also carries `receptor_array_id`,
  `receptor_array_size`, `receptor_array_members` and an `receptor_array_support` flag; receptor arrays
  occur for mating and non-mating receptors alike, so membership does not make
  a locus a mating receptor ([docs/receptor-arrays.md](docs/receptor-arrays.md)).
- A two-idiomorph result is a flag for review. It can come from homothallism,
  a mixed culture or a duplication.

## Authors

- Jason E. Stajich — University of California, Riverside
  ([jason.stajich@ucr.edu](mailto:jason.stajich@ucr.edu),
  ORCID [0000-0002-7591-0020](https://orcid.org/0000-0002-7591-0020))

## Citation

No paper describes MATPredict yet. Cite the software and the release you used.
The citation metadata is in [CITATION.cff](CITATION.cff); GitHub shows it under
"Cite this repository".

> Stajich JE. MATPredict: MAT locus identification in Fungi. Version 0.6.1.
> https://github.com/stajichlab/MATPredict

```bibtex
@software{matpredict,
  author  = {Stajich, Jason E.},
  title   = {MATPredict: MAT locus identification in Fungi},
  version = {0.6.1},
  url     = {https://github.com/stajichlab/MATPredict},
  year    = {2026}
}
```

Please also cite the publications behind the reference records you rely on.
Each record's `metadata.yaml` lists its source (PMID/DOI or accession).

## License

MATPredict is released under the MIT License. See [LICENSE](LICENSE).
