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
3. **`-O2` and `-O3` failed. Fixed.** The cause was `kroundup64` (`kseq.h`, `mppriv.h`).
   Its last step is `x |= x >> 32`. For a 32-bit `size_t`, as on wasm32, that shift is
   undefined. `kseq_read` calls the macro on `seq->seq.m` while growing its buffer, so
   the optimised build computed a bad size and corrupted the heap. The patch skips the
   `>>32` step when `size_t` is 4 bytes. After the fix, `-O1`, `-O2` and `-O3` all match
   native on the 468-query set. UBSan at `-O2` no longer reports the shift.
   The AddressSanitizer route was abandoned: it did not finish in a reasonable time.
4. **Misaligned 8-byte accesses. Open, benign here.** UBSan still reports 8-byte loads and
   stores on addresses that are only 4-byte aligned, in `align.c`, `sketch.c` and `misc.c`.
   The custom allocator's blocks are not 8-byte aligned for every use. Native x86 tolerates
   this, and WASM accesses are not trapped for it, so output matches. It is still undefined
   behaviour in C. Fixing it means changing the allocator's alignment, which is not done.
5. **Harness limits.** Arguments not present in `inputDir` are ignored, so an input
   filename typo fails with `ERROR during mapping`, not with a clear message.

## Not done

- No timing or memory measurement at `-O1` beyond these tiny inputs.
- No genome-scale test. The target is a 366 kb multi-record FASTA, not an assembly.
- exonerate is not built (Phase 2 continues with it).
