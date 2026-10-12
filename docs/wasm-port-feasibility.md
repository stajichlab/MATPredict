# Browser (WebAssembly) port: feasibility

Status: research note, revised after an independent review of the code. Nothing
here has been built or run. Items marked **verify** are assumptions to check in
a spike before relying on them.

## Question

Can `matpredict detect` run in a static HTML page, with the Python policy layer
converted or compiled to WebAssembly and the native search tools replaced or
compiled too?

## Short answer

Feasible for a genome-only v1, with one hard dependency: the BLAST+ `tblastn`
localization step (or a replacement for it). The Python orchestration can run
largely unchanged under Pyodide (CPython compiled to WebAssembly), but the
subprocess shim is more work than it first appears, and the bigger risks are
result parity and silent network fallbacks. Rust is only needed if a `tblastn`
replacement is required. Porting the Python policy layer to Rust is not
recommended.

## What the detect path depends on

### Native binaries (called through `subprocess`)

| Tool | Used for | Where | Browser plan |
|---|---|---|---|
| `makeblastdb`, `tblastn` | Genome-wide localization | `detect/search.py` (`search_localize`) | Core. Hardest piece. |
| `miniprot` | Gene-model polish on a sliced window | `detect/search.py` (`polish_with_miniprot`) | Core. Small C tool; likely the easiest WASM build. **verify** |
| `exonerate` | Gene-model polish (`--refine`) | `detect/search.py` (`polish_with_exonerate`) | Core. C, depends on glib. The code retries without `--refine` after a crash (`search.py:784-803`). |
| `diamond` | Proteome fast path (only with `--proteins`) | `detect/search.py` (`search_fast_path`) | Drop for v1. The genome-only branch (`pipeline.py:2142-2148`) does not use it. |
| `augustus` | Ab initio prediction | `pixi.toml` only | Drop. No call site in `src/`, `scripts/` or `tests/`. |
| `mafft` | Classifier build, offline only | `detect/classifier_build.py` | Drop. Not on the runtime path. |
| `datasets` | Genome download | `detect/genome_acquisition.py` | Drop. The browser takes an upload. |
| `zstd` | Decompressing the taxonomy table when no Python zstd is available | `db/local_taxonomy.py:60` | Replace: bundle `zstandard`, or ship the table uncompressed. |

`taxonkit` is **not** on the detect path. It is used by `rollout_aggregate` only
(`db/taxonomy.py:21`). Detect resolves lineage through
`default_lineage_taxids`, `default_lineage_phylum_name` and `default_genetic_code`
(`db/taxonomy.py:122-160`), which read a local table (`$MATPREDICT_TAXONOMY`,
a 25 MB zstd taxdump) or fall back to NCBI efetch (`db/taxonomy.py:49-54`).

### The subprocess shim

The wrappers take an injected `runner` callable, and the tests use that seam.
A browser build cannot just patch `subprocess.run` afterwards. The default
argument is bound when `detect/search.py` is imported
(`runner: Callable = subprocess.run`, at `search.py:272`, `:291`, `:435`, `:708`,
`:875`), so the shim must be installed before any MATPredict module is imported,
or the wrapped functions must be passed in explicitly. `run_pipeline`
(`pipeline.py:1962`) takes no `runner`, so the explicit route means injecting
functions into it.

The shim must also:

- accept `capture_output=True, text=True` and return an object with
  `.returncode`, `.stdout` and `.stderr` (`search.py:259`);
- map a WASM trap or abort to a negative return code, because the native code
  treats `returncode < 0` as a crash (`search.py:261`);
- start a fresh module instance per call, since an Emscripten `main()` is not
  safe to run twice; exonerate and miniprot run once or more per polished gene,
  so instantiation cost adds up;
- copy files between tools. Each separately compiled Emscripten module has its
  own filesystem. `makeblastdb` writes `genome_db.*` that `tblastn` then reads
  (`search.py:448-457`), and every window FASTA is written by MATPredict and read
  by a polish tool. A 60 MB genome copied on every call is a real cost. The
  alternative is to link the tools into one module.

Direct `subprocess` calls that bypass any runner, and matter here:
`db/local_taxonomy.py:60` (`zstd -dc`). Off the detect path: `report/pdf.py:97`
(`--pdf` only), `reads/kmers.py:34`, `detect/classifier_build.py:177`,
`detect/genome_acquisition.py:89`.

### Network: a silent behaviour change

The detect path reaches the network by default. Failed lookups fall back quietly,
and the fallbacks change results:

- a failed lineage lookup sends routing down the fallback chain
  (`detect/family_registry.py:651`);
- a failed genetic-code lookup falls back to table 1 (`detect/pipeline.py:2071-2080`),
  which runs a different search for CTG-clade yeasts.

So the browser must ship the taxonomy table and force offline behaviour
(**verify** the existing `MATPREDICT_OFFLINE` switch, or add one), or require
`--phylum` and an explicit genetic code. Otherwise browser results can differ
from native ones with no error. A bundled Fungi-only table is the expected
route.

### Python libraries on the detect path

- **biopython**: sequence objects, translation, genetic-code tables.
- **pyyaml**: report and metadata I/O.
- **pyhmmer**: idiomorph and paralog HMM scoring, imported lazily
  (`detect/classifier.py:138`). Only Mucoromycota ships HMMs, so Ascomycota and
  Basidiomycota runs never need it. **verify** a Pyodide-compatible build. HMMER3
  uses SIMD (SSE or NEON), so a WASM build is a real porting job if none exists.
- Standard library: pathlib, dataclasses, tempfile, gzip, re, hashlib, and others.

Not on the detect path: `jsonschema` (imported only by `db/schema.py:8`),
`duckdb`, `weasyprint`, `pygenomeviz`, `clinker` (the `viz` extra), `requests`
except through the network fallbacks above.

Two import hazards for the browser entry point:

- Import `MATPredict.detect.cli`, not `MATPredict.__main__`. `__main__` imports
  `db/cli.py`, which imports `build_duckdb` (duckdb) and `requests` at module
  scope (`db/cli.py:9-15`).
- `MATPredict/__init__.py:20` calls `importlib.metadata` on `MATPredict`. A copied
  source tree will fail at import, so the package must be installed as a wheel
  (for example with micropip).

### Concurrency

None on the detect path: no `multiprocessing`, `concurrent.futures`, or
threading was found in `src/MATPredict`. That fits a single-threaded WASM
runtime.

### Data

| Item | Size |
|---|---|
| Curated reference DB (`db/`) | 5.5 MB |
| Classifier directory (`db/Mucoromycota/classifiers/MAT`) | 372 KB |
| Package Python source (`src/MATPredict/**/*.py`) | about 875 KB |

167 `metadata.yaml` files: 166 curated records plus one in `db/candidates/`,
which searches exclude (`detect/family_registry.py:582`). 21 curated records
are in Mucoromycota. The data is small enough to ship with the page.

### Output

`report/html.py` already writes a self-contained HTML page with no network,
external CSS, fonts or scripts. It is close to a browser deliverable. Its PDF
path (WeasyPrint or Chrome) does not apply in the browser; use the page's print
stylesheet and the browser's "Save as PDF" instead.

## Options for the search tools

1. **Compile the C tools to WASM with Emscripten.** Keeps results identical to
   native. This matters because the thresholds were tuned on native scores; the
   code comments note that raw scores are comparable only within one tool.
   - miniprot: low to medium effort. **verify**
   - exonerate: medium to high, mainly because of glib.
   - BLAST+ (`tblastn`, `makeblastdb`): high. Large C++ build, no official WASM
     target. **verify** whether community WASM builds exist.
2. **Replace `tblastn` with a Rust/WASM translated search.** Plausible for about
   100 reference proteins against one genome, but hits will not match BLAST's.
   A replacement must reproduce the 11-column outfmt (`search.py:350`), the
   negative-frame and reversed-coordinate conventions for the minus strand,
   `-seg no`, the genetic-code option, and BLAST's default e-value and
   max-target cut-offs. The identity and bitscore gates and the rescue logic
   would then need recalibration against the regression panel.
2b. **Use miniprot as the genome-wide localizer.** Neither polish tool can run
   without a localizer, because polish windows come from `tblastn` hits
   (`pipeline.py:2143`, `_extract_window`). miniprot searches a whole genome
   by design and is being compiled anyway. This needs the same recalibration as
   option 2, but no new search code. Worth testing before option 2.
3. **Port the Python logic to Rust.** Not recommended. `detect/pipeline.py` alone
   is 3,318 lines of tuned policy (gates, confidence tiers, flank rules, rescue
   logic). A port would be a second implementation to keep in sync, and the
   regression panel would test the port rather than the code users run.

### Licenses

Of the bundled tools, exonerate is GPL (glib is LGPL), BLAST+ is public domain,
and miniprot is MIT. The licences of HMMER, pyhmmer and biopython are **verify**.
Check the combined bundle before shipping it.

## Recommended path

- **Phase 0, scope (v1):** genome upload only; `tblastn` (or miniprot as
  localizer, option 2b) plus miniprot and exonerate polish; taxid routing from a
  bundled Fungi-only lineage table in offline mode; no proteome path, classifier
  build, or genome download; HTML report output. Decide the classifier question
  now: either ship Ascomycota and Basidiomycota only, or drop classifier-dependent
  calls in the browser.
- **Phase 1, Pyodide spike:** install the `MATPredict` wheel in Pyodide, import
  `MATPredict.detect.cli`, and confirm biopython, pyyaml and pyhmmer load. Run the
  code path up to the first subprocess call with stubbed binaries, and time it
  on the *Phycomyces blakesleeanus* NRRL 1555 genome (taxid 4837, the README
  example).
- **Phase 2, binaries and shim:** build miniprot to WASM first, then exonerate,
  with the shim described above (import-time injection, trap mapping, per-call
  instances, file copying). Try the BLAST+ build; if it fails, test option 2b,
  then option 2.
- **Phase 3, parity:** the full regression panel (166 genomes) lives on
  `/bigdata` with the SLURM and NCBI-cache setup, not in the repo, and a browser
  cannot run it. Run a 10–20 genome subset through a Node or headless harness,
  against the native results from `scripts/regression_check.py`. Acceptance:
  identical calls, and score differences within a stated tolerance.
- **Phase 4, page:** a static HTML page that loads the Pyodide runtime, the WASM
  tools and the DB, and runs the job in a Web Worker so the UI stays responsive.
  Flow: upload genome, choose taxid or phylum, run, view report.

## Risks, in order of concern

1. **Result parity.** Gate thresholds depend on native scores. This is the main
   risk, larger than build effort.
2. **Silent network fallbacks.** Failed lineage or genetic-code lookups change
   routing and search without an error. Ship the table and force offline mode
   before anything else.
3. **Speed and memory.** Native medians from `scripts/run_regression_panel.sh:18-19`
   are about 620 s per genome for Ascomycota MAT, about 150 s for Basidiomycota
   and about 130 s for Mucoromycota. The clade-panel SLURM job recorded a maximum
   of 1,737 s, on a 1.71 Gb Pucciniales genome. tblastn already runs
   single-threaded natively, so the browser pays only the WASM slowdown, but
   expect tens of minutes per Ascomycota genome. A wasm32 module is limited to
   4 GB of memory, so the page needs a maximum genome size. Measure tblastn first.
4. **BLAST+ build.** Decides whether option 1 or option 2/2b is needed for the
   localization step.
5. **Shim cost.** Per-call instantiation, file copying between module filesystems,
   and trap handling. Measure in the spike; it could dominate runtime.
6. **Bundle size.** Pyodide plus binaries plus DB is likely tens of MB. **verify**
7. **pyhmmer in WASM.** Affects whether Mucoromycota calls are in scope.
8. **Licenses.** See above.

## Open questions for the curator

- Is a genome-only v1 acceptable, or does the browser need the proteome path
  (which would bring back diamond)?
- Is a small, documented score difference acceptable if a `tblastn` replacement
  or a miniprot-as-localizer change is needed, or must results match native
  exactly?
- Is the browser limited to Ascomycota and Basidiomycota, leaving out the
  Mucoromycota classifier?
- What is the maximum genome size the page should accept?

## Review record

The first version of this note was reviewed by an independent Opus subagent
against the code. The review found the default-argument binding of the runner,
the network fallbacks, the taxonkit and zstd misattribution, the `__init__`
metadata import, and the unrealistic full-panel parity plan. The spot checks
for those findings were done against the source before this revision.
