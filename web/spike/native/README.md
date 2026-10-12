# Phase 2: miniprot to WebAssembly

Builds miniprot 0.18-r281 (upstream commit `671db24`) with Emscripten and checks its
output against the native binary. Part of the browser port; see
[docs/wasm-port-feasibility.md](../../../docs/wasm-port-feasibility.md).

## Files

| File | Purpose |
|---|---|
| `build_miniprot_wasm.sh` | Builds the WASM module. Applies the patch to a scratch copy. Default `OPT=-O1`. |
| `patches/miniprot-r281-kthread-serial.patch` | Serial `kt_pipeline` under `__EMSCRIPTEN__` (no pthreads). |
| `miniprot_wasm_test.mjs` | Node harness: fresh instance per call, mounts the input files, prints stdout. |
| `sse_diff_test.c` | Differential test of the SSE intrinsics used by `nasw-sse.c`. |

## Build

Needs Emscripten (tested with 6.0.12) and a miniprot checkout at `671db24`.

```bash
source ~/emsdk/emsdk_env.sh
bash web/spike/native/build_miniprot_wasm.sh <miniprot-src> <out-dir>
node web/spike/native/miniprot_wasm_test.mjs <out-dir> <input-dir> --gff -t 1 target.fa query.faa
```

## Results

Against the native r281 build, on `testset/MAT_loci3.fasta` as target:

| Query set | Native | WASM `-O1` |
|---|---|---|
| 3 proteins (sexM, one Mucoromycota record) | 4 lines | 4 lines, identical |
| All 468 curated reference proteins | 161 lines | 161 lines, identical |

The native reference is also the x86_64 build under Rosetta (same output as arm64).

The index statistics match native exactly at every optimisation level tried
(3022 blocks, 194,206 kmer-block pairs, 137,447 distinct k-mers).

## Problems found and their status

1. **Silent no-op without pthreads. Fixed.** `kt_pipeline` always creates pthread
   workers, even with `-t 1`. Without pthreads, `pthread_create` fails and the mapping
   stage never runs, with no error. The patch replaces the pipeline with a serial loop
   under `__EMSCRIPTEN__`. It follows the same single-worker schedule as native `-t 1`.
2. **SSE4.1 emulation. Checked, not a fault.** The SSE intrinsics the kernel uses
   give identical results in WASM and on x86 (`sse_diff_test.c`, 24 of 24 lines).
   `-msse2` fails to compile under Emscripten, so the build uses `-msse4.1`, which
   selects the code path the x86 SSE4.1 build also uses.
3. **`-O2` and `-O3` fail. Open.** Symptoms: at `-O3`, heap corruption aborts inside
   `realloc` in `kseq_read`. At `-O2`, the run exits with a trap and produces only the
   header. `-O1` and `-O0` are correct on the tests above. UBSan at `-O2` reports:
   - misaligned 8-byte loads and stores in `sketch.c:57` and `misc.c:5`. The custom
     allocator returns 8-byte-aligned blocks, which native x86 tolerates;
   - a shift exponent of 32 for `size_t` at `ntseq.c:8` (the line is a macro
     expansion, so the source is not yet located). This is a 32-bit `size_t` issue,
     and it is the likeliest cause of the optimised-build failure.

   `-fno-strict-aliasing` does not fix it. Emscripten AddressSanitizer did not finish
   in a reasonable time. **Ship `-O1` for now.** Investigate the shift and the alignment
   before moving to a higher level.
4. **Harness limits.** Arguments not present in `inputDir` are ignored, so an input
   filename typo fails with `ERROR during mapping`, not with a clear message.

## Not done

- No timing or memory measurement at `-O1` beyond these tiny inputs.
- No genome-scale test. The target is a 366 kb multi-record FASTA, not an assembly.
- exonerate is not built (Phase 2 continues with it).
