# MMseqs2 for WebAssembly (localizer candidate)

Spike only. MMseqs2 is the proposed replacement for the `tblastn` localization step, because
the BLAST+ toolkit has no practical WebAssembly build. See `../../localizer/README.md` for the
recall results that motivate it.

## Status (2026-10-12)

- Configure: works with `emcmake`.
- Rust block-aligner crate: fails with Homebrew's Rust (no wasm32 standard library). Works with
  rustup's stable toolchain and `-DRust_CARGO_TARGET=wasm32-unknown-emscripten`.
- Link: upstream passes `-sRELOCATABLE=1`, which the current Emscripten rejects. Removed in the patch.
- Compile: `src/alignment/Fwbw.cpp` `std::max` calls fail on SIMD element types. Fixed with a cast
  in the patch. The cast does not change values.
- Build of `mmseqs` after these fixes: see the result recorded in the feasibility note.

## Build

```bash
bash build-mmseqs-wasm.sh ~/wasm-build/mmseqs
```

The script pins upstream commit `564f40d` and applies `mmseqs2-wasm.patch`.

## Not yet checked

- Runtime under Node/WebAssembly: threads (MMseqs2 defaults to multiple threads), memory
  mapping of database files, and output parity.
- Upstream does not test this target in CI as far as the repo shows (no Emscripten job in
  `azure-pipelines.yml`). Verify before relying on it.
