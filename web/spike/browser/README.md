# Browser-style detect: tool host and Python bridge (spike)

Runs `matpredict detect` under Pyodide with every external tool served by a WebAssembly module,
using the same JS host and Python bridge that a browser page would use. Node stands in for the
browser here.

## Pieces

| File | Purpose |
|---|---|
| `host/tool-host.mjs` | Runs one tool per call in a fresh in-memory module. Copies the files the call touches in from the Python filesystem, copies new or changed files back, maps traps to negative return codes. |
| `python/bridge.py` | Replaces `subprocess.run` before any MATPredict import. Routes `makeblastdb` and `tblastn` to the localizer, and `mmseqs`, `miniprot`, `exonerate` to the host. Waits on the host promise with `run_sync` (JSPI). |
| `../localizer/localizer.py` | The `tblastn`/`makeblastdb` stand-ins as an importable module (also used by the command-line stand-ins in `../localizer/bin`). |
| `run-detect.mjs` | Node harness: loads Pyodide and the wheel, mounts the database and genome, runs `detect`, copies the report out. |

## Module builds (in-memory filesystem, modular)

The tools must be built without `NODERAWFS` (so each instance has its own in-memory filesystem)
and as `MODULARIZE` modules (so the host can create one per call).

- miniprot: `build_miniprot_wasm.sh` with `EXTRA_EMCC="-sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=268435456 -sSTACK_SIZE=8388608"`.
- exonerate: relink the program with `-sMODULARIZE=1 -sEXPORT_NAME=createExonerate -sINVOKE_RUN=0 -sEXIT_RUNTIME=1 -sEXPORTED_RUNTIME_METHODS=callMain,FS -sENVIRONMENT=web,node` (no `NODERAWFS`).
- MMseqs2: relink `mmseqs` the same way, with `-sSTACK_SIZE=8MB -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB`.

These relinks are currently ad hoc, from scratch directories under `~/wasm-build`. They should go
into `build-mmseqs-wasm.sh` and the exonerate build script before this is reproducible.

## Results (Node, 2026-10-13)

- `detect` on `MAT_loci3.fasta` with all tools in WebAssembly: return code 0, about 31 s.
- Against the native-MMseqs2 localizer run: detected locus, families attempted and not-detected
  list identical; `detected_loci.gff3` identical; `suppressed_loci` has the same 153 entries as
  before (same 148-of-153 difference).
- Against the native-tblastn reference: the same differences as before (weaker recall on weak hits).

## Limits and open items

- **JSPI.** `run_sync` needs JSPI. Node 26 has it. Recent Chromium has it; other browsers need
  checking. Without JSPI the design needs a worker with `SharedArrayBuffer`, which needs
  cross-origin isolation headers.
- **Instance cost.** Every tool call creates a fresh module, and copies files in and out. The
  detect run made many calls and took 31 s on a 366 kb input. A genome-scale run will be dominated
  by this and by the 4 GB memory limit; not measured.
- **Copy scope.** The host copies every directory named by an absolute argument, so large
  directories are copied on every call. Acceptable for this spike; tighten before a real run.
- **Not browser-tested.** The Node harness shares the same code, but no browser has run it yet.
- **One test input.** `MAT_loci3.fasta` only. The pipeline needs real genomes to test.
