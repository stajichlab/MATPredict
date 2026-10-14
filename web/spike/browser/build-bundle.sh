#!/usr/bin/env bash
# Build the WebAssembly tools for the browser page into one bundle directory. Spike only.
#   miniprot   - gene-model polish (modular in-memory module)
#   mmseqs     - translated search for the localizer (modular in-memory module)
#   pyhmmer    - Pyodide wheel for the Mucoromycota classifier
# exonerate is built by web/spike/native/exonerate/build-exo-wasm.sh (it needs the native shim tree);
# see its README. Usage: build-bundle.sh <work-dir> <bundle-dir>
set -eo pipefail
WORK="$1"; BUNDLE="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
mkdir -p "$WORK" "$BUNDLE"

# miniprot: pinned r281 source with the WebAssembly patch, modular, in-memory filesystem.
[ -d "$WORK/miniprot" ] || git clone -q https://github.com/lh3/miniprot.git "$WORK/miniprot"
git -C "$WORK/miniprot" checkout -q 671db24
source "$HOME/emsdk/emsdk_env.sh" > /dev/null 2>&1
OPT=-O2 EXTRA_EMCC="-sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=268435456 -sSTACK_SIZE=8388608" \
  bash "$REPO/web/spike/native/build_miniprot_wasm.sh" "$WORK/miniprot" "$BUNDLE/miniprot"

# MMseqs2: pinned commit, patched, modular module.
bash "$REPO/web/spike/native/mmseqs2/build-mmseqs-wasm.sh" "$WORK/mmseqs" > "$WORK/mmseqs.log" 2>&1
mkdir -p "$BUNDLE/mmseqs"
cp "$WORK/mmseqs/build/src/mmseqs-mem.js" "$WORK/mmseqs/build/src/mmseqs-mem.wasm" "$BUNDLE/mmseqs/"

# pyhmmer: Pyodide wheel.
bash "$HERE/pyhmmer/build-pyhmmer-wasm.sh" "$WORK/pyhmmer" > "$WORK/pyhmmer.log" 2>&1
cp "$WORK/pyhmmer/dist/"*.whl "$BUNDLE/"

ls -la "$BUNDLE"
