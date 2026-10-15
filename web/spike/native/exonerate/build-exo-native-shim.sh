#!/usr/bin/env bash
# Native exonerate 2.4.0 against glib-lite. The stand-in pkg-config must come first on PATH.
# Automake/autoconf are disabled in make so that the tree is not regenerated (and configure
# is not re-run without the PATH override).
set -euo pipefail
SHIM=~/wasm-build/glib-lite-native
cd ~/wasm-build/exo-shim-native
make distclean >/dev/null 2>&1 || true
PATH="$SHIM:$PATH" ./configure > ../exo-shim-config.log 2>&1
grep "^PKG_CONFIG\|^glib_cflags\|^glib_libs" Makefile | cut -c1-140
make AUTOMAKE=: ACLOCAL=: AUTOCONF=: AUTOHEADER=: > ../exo-shim-make3.log 2>&1
echo "make=$?"
