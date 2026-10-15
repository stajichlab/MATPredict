#!/usr/bin/env bash
# Build MMseqs2 (translated search for the localizer) for WebAssembly, unthreaded. Spike only.
# Usage: build-mmseqs-wasm.sh <work-dir>
# Result: <work-dir>/build/src/mmseqs-mem.js (+ .wasm), a modular in-memory-filesystem module.
set -eo pipefail
WORK="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
PINNED=564f40d   # upstream master commit used for this build
mkdir -p "$WORK" && cd "$WORK"
[ -d MMseqs2 ] || git clone -q https://github.com/soedinglab/MMseqs2.git
# Reset to the pinned commit first, so a reused work directory (already patched) builds the same tree.
git -C MMseqs2 checkout -q -f "$PINNED"
git -C MMseqs2 reset -q --hard "$PINNED"
git -C MMseqs2 apply "$HERE/mmseqs2-wasm.patch"
# Rust must come from rustup: a distro or Homebrew rustc has no wasm32 standard library.
# (ci/setup-toolchains.sh installs it and the wasm32-unknown-emscripten target.)
source ~/emsdk/emsdk_env.sh >/dev/null 2>&1
rm -rf build && mkdir build && cd build
# OpenMP off: upstream's serial mode (no -fopenmp, so no atomics; ips4o uses its serial sort).
emcmake cmake "$WORK/MMseqs2" -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DREQUIRE_OPENMP=0 \
  -DRust_CARGO_TARGET=wasm32-unknown-emscripten
emmake make -k -j8 mmseqs
# Relink as a modular, in-memory-filesystem module for the browser host: no NODERAWFS (the host
# copies files in and out), an 8 MB stack for extractorfs (its Orf/TranslateNucl objects are ~42 KB
# each), and memory growth.
cd src
LINK=$(cat CMakeFiles/mmseqs.dir/link.txt)
eval "${LINK/-o mmseqs.js/-o mmseqs-mem.js} -sEXIT_RUNTIME=1 -sALLOW_MEMORY_GROWTH=1 -sMAXIMUM_MEMORY=4GB -sSTACK_SIZE=8MB -sMODULARIZE=1 -sEXPORT_NAME=createMmseqs -sINVOKE_RUN=0 -sEXPORTED_RUNTIME_METHODS=callMain,FS -sENVIRONMENT=web,node"
ls -la mmseqs-mem.js mmseqs-mem.wasm
