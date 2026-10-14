#!/usr/bin/env bash
# Build LOSAT (pure-Rust BLAST reimplementation, MIT) as a WASI command for the localizer. Spike only.
# Usage: build-losat-wasi.sh <work-dir>
# Result: <work-dir>/LOSAT/LOSAT/target/wasm32-wasip1/release/LOSAT.wasm (and a native binary in target/release)
#
# Notes:
#  - The crate is in the LOSAT/ subdirectory of the repository (Cargo.toml lives there).
#  - Needs the wasm32-wasip1 standard library in the rustup toolchain: rustup target add wasm32-wasip1.
#  - rustup's rust-lld loads libLLVM.dylib from its own lib directory; the build needs that directory on
#    DYLD_FALLBACK_LIBRARY_PATH on this machine.
set -eo pipefail
WORK="$1"
mkdir -p "$WORK" && cd "$WORK"
[ -d LOSAT ] || git clone -q --depth 1 https://github.com/satoshikawato/LOSAT.git LOSAT
cd LOSAT/LOSAT
TC="$HOME/.rustup/toolchains/stable-aarch64-apple-darwin"
"/opt/homebrew/opt/rustup/bin/rustup" target add wasm32-wasip1 --toolchain stable-aarch64-apple-darwin > /dev/null 2>&1 || true
PATH="$TC/bin:$PATH" DYLD_FALLBACK_LIBRARY_PATH="$TC/lib" cargo build --release --target wasm32-wasip1 --no-default-features
PATH="$TC/bin:$PATH" cargo build --release
ls -la target/wasm32-wasip1/release/LOSAT.wasm target/release/LOSAT
