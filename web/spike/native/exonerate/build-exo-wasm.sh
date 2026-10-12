#!/usr/bin/env bash
# Exonerate 2.4.0 for WebAssembly against glib-lite (wasm). See the README in web/spike/native.
set -eo pipefail
source ~/emsdk/emsdk_env.sh >/dev/null 2>&1
cd ~/wasm-build/exo-wasm
make distclean >/dev/null 2>&1 || true
PATH=~/wasm-build/glib-lite-wasm:$PATH emconfigure ./configure --host=x86_64-pc-linux-gnu > ../exo-wasm-config.log 2>&1
grep "^PKG_CONFIG\|^glib_cflags\|^glib_libs" Makefile | cut -c1-140
# The bootstrapper is a native host tool. It writes the compiled C4 models with emcc/emar.
mkdir -p src/model
cp ~/wasm-build/exo-shim-native/src/model/bootstrapper src/model/bootstrapper
# The bootstrapper keeps its compiled objects in a codegen directory fixed by host type under the
# native source root, and reuses them if present. Clear it so the objects are rebuilt with emcc.
rm -rf ~/wasm-build/exo-shim-native/codegen
(cd src/program && rm -f c4_model_archive.a c4_model_archive.h && \
  PATH=~/wasm-build/tools:$PATH CC=emcc AR=emar C4_AR=emar CFLAGS="-O2" ../model/bootstrapper --compiled no > ../../../bootstrap-wasm.log 2>&1)
touch src/model/bootstrapper src/program/c4_model_archive.a src/program/c4_model_archive.h
# BOOTSTRAPPER=true: the generated C4 header and archive are already current, so the rule that
# would run the bootstrapper (as wasm, which cannot run on the host) does nothing.
# Patch: no-thread job queue (repo patch). -N skips it if already applied.
patch -p1 -N -d ~/wasm-build/exo-wasm < /Users/jstajich/projects/MATPredict/web/spike/native/patches/exonerate-240-jobqueue-wasm.patch > /dev/null || true
# Link flags for Node: the real filesystem, and enough memory (exonerate OOMs at the 16 MB default).
# Appended to the program Makefile's LDFLAGS (configure regenerates that file on every run above).
sed -i.bak 's#^LDFLAGS = \(.*\)#LDFLAGS = \1 -sNODERAWFS=1 -sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=268435456 -sSTACK_SIZE=8388608#' src/program/Makefile
grep -c "NODERAWFS" src/program/Makefile
make -k BOOTSTRAPPER=/usr/bin/true AUTOMAKE=: ACLOCAL=: AUTOCONF=: AUTOHEADER=: CC=emcc AR=emar RANLIB=emranlib > ../exo-wasm-make.log 2>&1 || true
ls -la src/program/exonerate* 2>&1 || { echo "NO exonerate"; grep " error" ../exo-wasm-make.log | head -10 | cut -c1-200; }
