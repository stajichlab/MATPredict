#!/usr/bin/env bash
# Regenerate exonerate's compiled C4 models for WebAssembly, using the native bootstrapper
# with the emscripten compiler and archiver.
set -euo pipefail
source ~/emsdk/emsdk_env.sh >/dev/null 2>&1
cd ~/wasm-build/exo-wasm/src/program
rm -f c4_model_archive.a c4_model_archive.h
CC=emcc AR=emar C4_AR=emar CFLAGS="-O2" \
  ~/wasm-build/exo-shim-native/src/model/bootstrapper --compiled no
ls -la c4_model_archive.a c4_model_archive.h
