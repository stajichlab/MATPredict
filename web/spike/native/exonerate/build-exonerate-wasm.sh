#!/usr/bin/env bash
# Build exonerate 2.4.0 (upstream commit 8397592) for WebAssembly against glib-lite, as a modular
# in-memory-filesystem module for the browser host. Spike only.
#
# Usage: build-exonerate-wasm.sh <work-dir>
# Result: <work-dir>/out/exonerate-mem.js and out/exonerate.wasm (the loader expects that name)
#
# Steps:
#   1. glib-lite (the GLib subset exonerate needs) for native and WebAssembly, with a pkg-config stand-in.
#   2. native exonerate against glib-lite: its bootstrapper generates the C4 model sources.
#   3. WebAssembly exonerate: the native bootstrapper is run with emcc/emar to build the models, then
#      make builds the program with the job-queue patch and a modular link.
# The bootstrapper keeps objects in a codegen directory under its own source root, so each build uses its own copy.
set -eo pipefail
WORK="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
GLIB="$REPO/web/spike/native/glib-lite"
PATCH="$REPO/web/spike/native/patches/exonerate-240-jobqueue-wasm.patch"
mkdir -p "$WORK" && cd "$WORK"

# 1. glib-lite, native and WebAssembly.
mkdir -p glib-native glib-wasm
cp "$GLIB/glib.h" "$GLIB/glib-lite.c" "$GLIB/pkg-config" glib-native/
cp "$GLIB/glib.h" "$GLIB/glib-lite.c" "$GLIB/pkg-config" glib-wasm/
(cd glib-native && clang -O2 -c glib-lite.c -o glib-lite.o && ar rcs libglib-lite.a glib-lite.o)
source "$HOME/emsdk/emsdk_env.sh" > /dev/null 2>&1
(cd glib-wasm && emcc -O2 -c glib-lite.c -o glib-lite.o && emar rcs libglib-lite.a glib-lite.o)

# Source: pinned 2.4.0. The full history is needed to find the commit.
[ -d exonerate-src ] || git clone -q https://github.com/nathanweeks/exonerate.git exonerate-src
git -C exonerate-src checkout -q 8397592

# 2. native build against glib-lite. The stand-in pkg-config must be first on PATH: configure calls the bare command.
rm -rf native && cp -R exonerate-src native
cd native
PATH="$WORK/glib-native:$PATH" ./configure > ../configure-native.log 2>&1
make AUTOMAKE=: ACLOCAL=: AUTOCONF=: AUTOHEADER=: > ../make-native.log 2>&1
cd "$WORK"
test -x native/src/program/exonerate

# 3. WebAssembly build.
rm -rf wasm && cp -R exonerate-src wasm
cd wasm
make distclean > /dev/null 2>&1 || true
PATH="$WORK/glib-wasm:$PATH" emconfigure ./configure --host=x86_64-pc-linux-gnu > ../configure-wasm.log 2>&1
mkdir -p src/model
cp "$WORK/native/src/model/bootstrapper" src/model/bootstrapper
rm -rf "$WORK/native/codegen"   # the bootstrapper reuses objects in this directory; clear the native ones
mkdir -p "$WORK/tools" && printf '#!/bin/sh\nexec emranlib "$@"\n' > "$WORK/tools/ranlib" && chmod +x "$WORK/tools/ranlib"
(cd src/program && rm -f c4_model_archive.a c4_model_archive.h && \
  PATH="$WORK/tools:$PATH" CC=emcc AR=emar C4_AR=emar CFLAGS="-O2" ../model/bootstrapper --compiled no > "$WORK/bootstrap-wasm.log" 2>&1)
touch src/model/bootstrapper src/program/c4_model_archive.a src/program/c4_model_archive.h
patch -p1 -i "$PATCH" > /dev/null
# Modular in-memory module: no NODERAWFS (the host copies files), memory growth, and a stack big enough for exonerate.
sed -i.bak 's#^LDFLAGS = \(.*\)#LDFLAGS = \1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=268435456 -sSTACK_SIZE=8388608 -sEXIT_RUNTIME=1 -sMODULARIZE=1 -sEXPORT_NAME=createExonerate -sINVOKE_RUN=0 -sEXPORTED_RUNTIME_METHODS=callMain,FS -sENVIRONMENT=web,node#' src/program/Makefile
make -k BOOTSTRAPPER=/usr/bin/true AUTOMAKE=: ACLOCAL=: AUTOCONF=: AUTOHEADER=: CC=emcc AR=emar RANLIB=emranlib > ../make-wasm.log 2>&1 || true
cd "$WORK"
mkdir -p out && cp wasm/src/program/exonerate out/exonerate-mem.js && cp wasm/src/program/exonerate.wasm out/exonerate.wasm  # the JS loads exonerate.wasm by name
ls -la out
