#!/usr/bin/env bash
# Sync glib-lite from the repo and build libglib-lite.a for native and WebAssembly.
set -eo pipefail
SRC=/Users/jstajich/projects/MATPredict/web/spike/native/glib-lite
build() {  # build <name> <cc> <ar>
  local dir=~/wasm-build/glib-lite-$1
  mkdir -p "$dir"
  cp "$SRC/glib.h" "$SRC/glib-lite.c" "$SRC/pkg-config" "$dir/"
  (cd "$dir" && rm -f glib-lite.o libglib-lite.a && $2 -O2 -c glib-lite.c -o glib-lite.o && $3 rcs libglib-lite.a glib-lite.o)
  echo "built $dir/libglib-lite.a"
}
build native clang ar
source ~/emsdk/emsdk_env.sh >/dev/null 2>&1
build wasm emcc emar
