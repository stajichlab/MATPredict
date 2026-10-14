# pyhmmer for Pyodide (in progress, not working yet)

pyhmmer 0.12.3 (the version MATPredict pins) built for Pyodide 0.27.7 (CPython 3.12, Emscripten 3.1.58).

## Build

- Sources: pyhmmer 0.12.3 sdist plus `pyhmmer-pyodide.patch` (apply with `patch -p1` in the unpacked sdist).
- Build host: Python 3.12 (Cython generates code for the build host's Python version).
- Toolchain: Pyodide build environment 0.27.7 (`xbuildenv-0.27.7.tar.bz2`), pyodide-build 0.27.3,
  emsdk 3.1.58. The layout was assembled by hand because the installer's index URL is gone.
- Flags: `CFLAGS=-I<shim>` (x86intrin.h stand-in), `LDFLAGS=-sWASM_BIGINT`.
- Result: builds `pyhmmer-0.12.3-cp312-cp312-pyodide_2024_0_wasm32.whl`.

## Status (update 2026-10-13)

- Step 1 (cause found): liblibhmmer and liblibeasel call each other through Pyodide's dynamic-link
  proxy. Easel's `esl_index_t` is int64, so many of these calls pass 64-bit values, which JavaScript
  cannot convert without BigInt. The search then fails with `type incompatibility`.
- Step 2 (done): Easel and HMMER are built as static libraries (CMake `STATIC`), so there are no
  cross-library calls. The search now runs to completion in Pyodide.
- New blocker: the MSV filter (`F1`) rejects sequences in WebAssembly that native keeps. With
  `F1=1` the WebAssembly hits match native to six decimals. Ruled out: unsigned byte SIMD intrinsics
  (`sse_epu8_test.c`, identical to x86) and exp/log math (`math_test.c`, identical). Not yet checked:
  the byte-score profile construction (`p7_oprofile`) and the MSV filter's own arithmetic.
- `filter-variants.py` runs the default and loosened filter settings in Pyodide.

## MSV filter investigation (2026-10-13)

Reference points (3 curated sequences, `sexM.hmm`):
- native arm64 pyhmmer (NEON): default filters 3 hits; with `F1=1` 3 hits.
- x86_64 pyhmmer 0.12.3 under Rosetta (SSE): default filters 3 hits.
- WebAssembly static build: default 0 hits; `F1=1` 3 hits; `F2=1`, `F3=1`, `bias_filter=False` still 0.
  So only the MSV stage (`F1`) differs, and the WebAssembly MSV p-value or score is wrong.

Ruled out (WebAssembly results identical to x86):
- unsigned and signed SIMD intrinsics used by `msvfilter.c` (`sse_epu8_test.c`, `sse_msv_test.c`)
- exp, log, rounding functions (`math_test.c` in the scratch tree, `rounding_test.c`)
- evalue parameters (`params.py`: MSV, Viterbi and Forward values identical)
- optimisation level (the `-O0` build gives the same result)

Patch tests (`run-variants.sh`), each compared with the same filter-variant test:
- revert the SSE header change: still 0 default hits (F1=1 gives 3);
- revert the stopwatch replacements: still 0 default hits;
- revert the `cpuid` stand-in: does not compile, so the original cannot be restored as a test;
- `cpuid` stand-in reporting SSE and SSE4.1 features: still 0 default hits.
So none of the patches is the cause.

Remaining suspects: the static copies of Easel and HMMER. The search runs inside one
extension, so duplication is unlikely to matter for the search itself, which makes it a weaker
candidate. The most likely remaining cause is code-level: a difference in how the emulated SSE
code runs as a whole, not any single intrinsic. Next step: a generic (non-SIMD) MSV path, or
per-stage score dumps from both builds, to locate the first differing value.

## Status (earlier)

- Imports, HMM reading, sequence reading and `Pipeline` creation work in Pyodide.
- `Pipeline.search_hmm` fails with `TypeError: type incompatibility when transforming from/to JS`.
  The call goes through Pyodide's dynamic-link stub to a function in pyhmmer's own shared library
  (`liblibeasel.so` or `liblibhmmer.so`) with a 64-bit argument. The stub mechanism and the link
  flags were checked; the specific function is not yet identified.
- The native comparison (classifier scores, `native-scores.txt`) has not been matched.

## Patches in the patch file, and why

- `esl_cpu.c`: no `cpuid` on WebAssembly; report no optional features.
- `esl_sse.c`, `esl_sse.h`: include `smmintrin.h` on Emscripten (the shim header does this for the rest).
- `platform/CMakeLists.txt`, `platform/__init__.pxd`: Emscripten uses the Linux platform file.
- `esl_stopwatch.c`: `times()` and `time()` stand-ins based on `clock()` (timing only, not results).
- `pyproject.toml`: Cython `>=3.0.10,<3.1`; no `wheel.py-api` (the stable-ABI build hides functions Cython needs).

## Next steps

1. Identify the 64-bit cross-library call. Compare the function signatures imported by
   `liblibeasel.so` / `liblibhmmer.so` with those in the search path.
2. Try linking Easel and HMMER statically into the extension modules (no separate shared libraries),
   which removes the cross-module calls.
3. Run the same score comparison (`score.py`, `pyodide-score.mjs`) once search runs.

`psutil.py` is a one-function stand-in (pyhmmer only calls `cpu_count`), and `score.py` runs
the search single-threaded (`cpus=1`), which WebAssembly needs.
