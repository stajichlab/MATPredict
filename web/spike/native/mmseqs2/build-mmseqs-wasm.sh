#!/usr/bin/env bash
# Build MMseqs2 (translated search for the localizer) for WebAssembly, unthreaded. Spike only.
# Usage: build-mmseqs-wasm.sh <work-dir>
# Result: <work-dir>/build/src/mmseqs-stk.js (+ .wasm). Run it with node (NODERAWFS, so it sees host files).
set -eo pipefail
WORK="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
PINNED=564f40d   # upstream master commit used for this build
mkdir -p "$WORK" && cd "$WORK"
[ -d MMseqs2 ] || git clone -q https://github.com/soedinglab/MMseqs2.git
git -C MMseqs2 checkout -q "$PINNED"
git -C MMseqs2 apply "$HERE/mmseqs2-wasm.patch"
# Rust must come from rustup: Homebrew's rustc has no wasm32 standard library.
export PATH="$HOME/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH"
source ~/emsdk/emsdk_env.sh >/dev/null 2>&1
rm -rf build && mkdir build && cd build
# OpenMP off: upstream's serial mode (no -fopenmp, so no atomics; ips4o uses its serial sort).
emcmake cmake "$WORK/MMseqs2" -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DREQUIRE_OPENMP=0 \
  -DRust_CARGO_TARGET=wasm32-unknown-emscripten
emmake make -k -j8 mmseqs
# Relink with the runtime settings the extractorfs step needs: an 8 MB stack (its Orf/TranslateNucl
# objects are ~42 KB each) and memory growth. NODERAWFS lets Node read host files.
cd src
LINK=$(cat CMakeFiles/mmseqs.dir/link.txt)
eval "${LINK/-o mmseqs.js/-o mmseqs-stk.js} -sNODERAWFS=1 -sEXIT_RUNTIME=1 -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=8MB"
ls -la mmseqs-stk.js mmseqs-stk.wasm
