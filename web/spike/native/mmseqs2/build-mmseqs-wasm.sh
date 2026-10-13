#!/usr/bin/env bash
# Build MMseqs2 (translated search for the localizer) for WebAssembly. Spike only.
# Usage: build-mmseqs-wasm.sh <work-dir>
set -eo pipefail
WORK="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
PINNED=564f40d   # upstream commit used for this build
mkdir -p "$WORK" && cd "$WORK"
[ -d MMseqs2 ] || git clone -q https://github.com/soedinglab/MMseqs2.git
git -C MMseqs2 checkout -q "$PINNED"
git -C MMseqs2 apply "$HERE/mmseqs2-wasm.patch"
# Rust must come from rustup: Homebrew's rustc has no wasm32 standard library.
export PATH="$HOME/.rustup/toolchains/stable-aarch64-apple-darwin/bin:$PATH"
source ~/emsdk/emsdk_env.sh >/dev/null 2>&1
rm -rf build && mkdir build && cd build
emcmake cmake "$WORK/MMseqs2" -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DRust_CARGO_TARGET=wasm32-unknown-emscripten
emmake make -j4 mmseqs
