# Phase 1 spike: MATPredict detect under Pyodide

Prototype on branch `wasm-spike`. Not wired into the package or the release. See
[docs/wasm-port-feasibility.md](../../docs/wasm-port-feasibility.md) for the plan.

## What this checks

1. `import_check.mjs`: the MATPredict wheel imports under Pyodide 0.27.7 (CPython 3.12)
   with biopython, pyyaml, jsonschema, zstandard and requests from the Pyodide CDN.
2. `seam_probe.mjs`: `matpredict detect` runs through routing and the reference
   build, and its first external tool call is observed. `subprocess.run` is replaced
   before any MATPredict import, and every tool call is logged and refused.

## Run

```bash
python3 -m pip wheel . --no-deps -w web/spike/dist     # from the repo root
cd web/spike && npm install
node import_check.mjs dist/matpredict-0.6.3-py3-none-any.whl
node seam_probe.mjs
```

## Results (2026-10-12)

- All seven detect-path modules import: `MATPredict`, `detect.search`,
  `detect.pipeline`, `detect.cli`, `detect.classifier`, `db.taxonomy`,
  `db.local_taxonomy`.
- pyhmmer is **not** in the Pyodide 0.27.7 package set. `detect.classifier` imports
  without it because pyhmmer is loaded lazily, so the classifier path is untested.
- Wheel install: micropip refuses a local wheel path under Node. The probe unpacks
  the wheel onto `sys.path` instead. Browser code can use the same approach or serve
  the wheel over HTTP.
- The seam works. With `MATPREDICT_DB_ROOT=/db` and `--phylum Ascomycota`, the first
  external call is `makeblastdb -in <genome> -dbtype nucl -out <tmp>/genome_db`,
  which means the injected default argument was honoured.
- No network calls were needed on the `--phylum` route. The `--taxid` route was not
  tested.

## Not done

- No Emscripten build of any native tool (no `emcc` on this machine).
- Input is `testset/MAT_loci3.fasta`, a multi-record FASTA of GenBank sequences, not
  an assembly. It is enough to reach the search stage. It says nothing about speed
  or memory.
- The `--taxid` path, the classifier path, and the report path were not run.
