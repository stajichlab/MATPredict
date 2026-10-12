#!/usr/bin/env bash
# Build miniprot 0.18-r281 (commit 671db24) to WebAssembly with Emscripten.
#
# Usage: build_miniprot_wasm.sh <miniprot-source-dir> <output-dir>
#   OPT=-O1 (default). -O2 and -O3 build but fail at run time (see README.md).
#   EXTRA_EMCC adds flags, e.g. EXTRA_EMCC="-g -fsanitize=undefined" (word-split).
# Assumes emsdk is activated: source ~/emsdk/emsdk_env.sh
set -euo pipefail

SRC_DIR="$(cd "$1" && pwd)"
OUT_DIR="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"
PATCH="$HERE/patches/miniprot-r281-kthread-serial.patch"

# Work on a scratch copy so the patch never touches the upstream checkout.
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cp -R "$SRC_DIR"/. "$WORK"/
patch -p0 -d "$WORK" -i "$PATCH" >/dev/null

mkdir -p "$OUT_DIR"

# Everything except example.c (its own main). kt_pipeline is serial under
# __EMSCRIPTEN__ (the patch), kt_for is serial because miniprot runs with -t 1.
SOURCES=(align.c bseq.c chain.c format.c hit.c index.c kalloc.c kthread.c main.c
         map.c misc.c nasw-sse.c nasw-tab.c ntseq.c options.c sketch.c sys.c)

# -msse4.1 selects the __SSE4_1__ branch in nasw-sse.c. Its SSE2-only branch calls
# _mm_insert_epi16, which Emscripten's emmintrin.h does not declare.
# shellcheck disable=SC2086
cd "$WORK"
emcc "${OPT:--O1}" ${EXTRA_EMCC:-} -std=c99 -msse4.1 -msimd128 \
  -Wno-unused-but-set-variable -Wno-unused-variable \
  "${SOURCES[@]}" \
  -o "$OUT_DIR/miniprot.mjs" \
  -sUSE_ZLIB=1 \
  -sMODULARIZE=1 -sEXPORT_NAME=createMiniprot \
  -sINVOKE_RUN=0 -sEXIT_RUNTIME=1 \
  -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=64MB -sSTACK_SIZE=1MB \
  -sENVIRONMENT=web,node -sFORCE_FILESYSTEM=1 \
  -sEXPORTED_RUNTIME_METHODS=callMain,FS \
  -lm
ls -la "$OUT_DIR"
