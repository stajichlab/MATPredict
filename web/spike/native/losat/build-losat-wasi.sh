#!/usr/bin/env bash
# Build LOSAT (pure-Rust BLAST reimplementation, MIT) as a WASI command for the localizer. Spike only.
# Usage: build-losat-wasi.sh <work-dir>
# Result: <work-dir>/LOSAT/LOSAT/target/wasm32-wasip1/release/LOSAT.wasm (and a native binary in target/release)
#
# Notes:
#  - The crate is in the LOSAT/ subdirectory of the repository (Cargo.toml lives there).
#  - Needs the wasm32-wasip1 standard library in the rustup toolchain (ci/setup-toolchains.sh installs it).
#  - On macOS, rustup's rust-lld loads libLLVM.dylib from the toolchain lib directory (handled below).
set -eo pipefail
WORK="$1"
mkdir -p "$WORK" && cd "$WORK"
# Pinned: the commit the browser build was tested against (PR #183 merge).
[ -d LOSAT ] || git clone -q https://github.com/satoshikawato/LOSAT.git LOSAT
git -C LOSAT checkout -q 730b778eeea784e2a76df8f1c3a4e5f69ef14c5d
cd LOSAT/LOSAT
rustup target add wasm32-wasip1 > /dev/null 2>&1 || true
# macOS only: rustup's rust-lld loads libLLVM.dylib from the toolchain lib directory.
if [ "$(uname)" = Darwin ]; then export DYLD_FALLBACK_LIBRARY_PATH="$(rustc --print sysroot)/lib"; fi
cargo build --release --target wasm32-wasip1 --no-default-features
cargo build --release
ls -la target/wasm32-wasip1/release/LOSAT.wasm target/release/LOSAT
