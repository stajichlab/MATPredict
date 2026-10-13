# MMseqs2 for WebAssembly (localizer candidate)

Spike only. MMseqs2 is the proposed replacement for the `tblastn` localization step, because
the BLAST+ toolkit has no practical WebAssembly build. See `../../localizer/README.md` for the
recall results that motivate it.

## Status (2026-10-12)

Builds to WebAssembly with zero errors and runs under Node. Unthreaded, OpenMP off.

Verified:
- Build from a clean clone (`build-mmseqs-wasm.sh`), 210 objects, no errors.
- `mmseqs version` and `createdb` run under Node.
- Translated search output matches a native build of the same commit: 2418 of 2418 lines
  identical after sort, on 195 reference proteins against `testset/MAT_loci3.fasta`.
- Equivalence was checked with `mmseqs search`'s six underlying module calls (below).

Fixes applied (in `mmseqs2-wasm.patch`):
- Upstream CMake passes `-sRELOCATABLE=1` to the Rust link, which current Emscripten rejects.
- `find_package(OpenMP)` is skipped for Emscripten, and `-DREQUIRE_OPENMP=0`. With OpenMP on,
  the build needs atomics, and the Rust std object lacks them.
- `kmermatcher.cpp`: explicit template instantiations that rely on the defaulted `IncludeSeqLen`
  argument create duplicate definitions in clang. Spelling out all template arguments fixes it.
  The same error occurs with the host clang, so it is not Emscripten-specific.
- `Fwbw.cpp`: `static_cast<float>` on two SIMD element accesses for `std::max`.

Runtime notes:
- `mmseqs search` and `easy-search` are shell workflows run by `execvp`, which fails in WebAssembly.
  Call the modules directly instead: `createdb`, `createdb --dbtype 2`, `extractorfs`, `prefilter`,
  `align`, `offsetalignment`, `convertalis`. Parameters are in `calls.log` in the review scratch tree.
- `extractorfs` needs an 8 MB stack (Emscripten's default is 64 KB and its objects are ~42 KB).
- Node runs need `-sNODERAWFS=1`, so the module can read host files.

## Build

```bash
bash build-mmseqs-wasm.sh ~/wasm-build/mmseqs
```

The script pins upstream commit `564f40d`, applies `mmseqs2-wasm.patch`, and relinks with the runtime flags.
Output: `<work>/build/src/mmseqs-stk.js` and `.wasm`.

## Not yet checked

- The localizer stand-in (`web/spike/localizer/bin/tblastn`) still calls `mmseqs search`; it must use the module calls above.

- Runtime under Node/WebAssembly: threads (MMseqs2 defaults to multiple threads), memory
  mapping of database files, and output parity.
- Upstream does not test this target in CI as far as the repo shows (no Emscripten job in
  `azure-pipelines.yml`). Verify before relying on it.
