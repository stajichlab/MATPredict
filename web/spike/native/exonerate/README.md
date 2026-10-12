# Phase 2: exonerate 2.4.0 to WebAssembly

Builds exonerate 2.4.0 (upstream commit `8397592`, the pristine 2.4.0 import) for WebAssembly,
with glib replaced by `../glib-lite`. Spike only; not part of the package or release.

## Pieces

| File | Purpose |
|---|---|
| `../glib-lite/` | Subset of GLib that exonerate calls, plus a `pkg-config` stand-in. |
| `build-glib-lite.sh` | Builds `libglib-lite.a` for native and WebAssembly. |
| `build-exo-native-shim.sh` | Native exonerate against glib-lite (used as a check on the shim). |
| `bootstrap-wasm.sh` | Regenerates exonerate's compiled C4 models with emcc. |
| `build-exo-wasm.sh` | Full WebAssembly build. |
| `tools/ranlib` | `emranlib` wrapper; the bootstrapper calls a plain `ranlib`. |
| `../patches/exonerate-240-jobqueue-wasm.patch` | Runs jobs inline when built for WebAssembly (no threads). |

The scripts assume the working tree in `~/wasm-build` (exonerate checked out at `8397592`,
emsdk in `~/emsdk`). They are kept here so the work can be reproduced, not for direct use.

## Things that were not obvious

- The configure step calls the bare `pkg-config`, not `$PKG_CONFIG`. The stand-in must be
  first on `PATH`. Without that, the "glib-lite" build silently links Homebrew glib.
- The bootstrapper hard-codes the native `ranlib` and a codegen directory under the native
  source root, and reuses objects found there. The WebAssembly bootstrap needs `tools/ranlib`
  and a cleared codegen directory.
- Exonerate's build helpers cannot run as WebAssembly, so `BOOTSTRAPPER=/usr/bin/true` is
  passed to make; the generated header and archive come from the wasm bootstrap run.
- With `USE_PTHREADS`, the job queue starts threads that never run without pthreads, and the
  queue spins. The patch runs each job in `JobQueue_submit` instead. The non-threaded branch
  of `jobqueue.c` never compiled in the stock build, so it had not been tested.
- The output order of records is not stable across runs even for the native build, so
  comparisons are of sorted lines.
- Node needs `NODERAWFS` to read input files, and exonerate needs more than the 16 MB default
  memory.

## Results (2026-10-12)

| Check | Result |
|---|---|
| glib-lite native vs real glib native, 468 curated proteins | Sorted lines identical (16,123 lines each). |
| WebAssembly vs native 2.4.0, 15 proteins | Sorted lines identical. |
| WebAssembly vs native 2.4.0, 468 proteins | Sorted lines identical, excluding the hostname line. 30 s in WebAssembly, about 80 s native. |

Test target: `testset/MAT_loci3.fasta` (366 kb, multi-record). Not a genome-scale test.

## Open

- Only `protein2genome --refine region` was compared, the form `polish_with_exonerate`
  uses. Other models and options are untested.
- The shim covers the functions this build needs. Other glib functions would fail at link time.
- Exonerate master has 50 commits after 2.4.0; this build pins 2.4.0.
- No timing or memory measurement beyond the one run above.
- Integration with `search.py` (routing the `exonerate` calls to the WebAssembly module) is not done.
