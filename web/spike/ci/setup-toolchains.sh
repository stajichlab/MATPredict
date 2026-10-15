#!/usr/bin/env bash
# Install the toolchains the WebAssembly bundle needs, at the versions it was built and tested with.
# Linux (GitHub Actions, the Docker test) and macOS. Idempotent: existing installs are reused.
#
#   emsdk 6.0.12        -> $HOME/emsdk        (miniprot, MMseqs2, exonerate, LOSAT's native side)
#   emsdk 3.1.58        -> $HOME/emsdk-3.1.58 (pyhmmer: Pyodide 0.27.7 is built with it)
#   Rust 1.99.0 + wasm32-wasip1, wasm32-unknown-emscripten   (LOSAT, MMseqs2)
#   Python 3.12 venv with pyodide-build 0.27.3, $PYHMMER_TOOLS/pb312   (pyhmmer)
#   Pyodide 0.27.7 xbuildenv, $PYHMMER_TOOLS/xbuildenv               (pyhmmer)
#
# Usage: setup-toolchains.sh          (PYHMMER_TOOLS defaults to $HOME/pyhmmer-tools)
# Then:  export PYHMMER_TOOLS=... and run browser/build-bundle.sh <work-dir> <bundle-dir>
set -eo pipefail
TOOLS="${PYHMMER_TOOLS:-$HOME/pyhmmer-tools}"
EMSDK_NEW=6.0.12
EMSDK_OLD=3.1.58
RUST=1.99.0
PYODIDE=0.27.7
PYODIDE_BUILD=0.27.3

install_emsdk() {  # <version> <dir>
  local version="$1" dir="$2"
  [ -d "$dir" ] || git clone -q https://github.com/emscripten-core/emsdk.git "$dir"
  (cd "$dir" && git fetch -q --tags && git checkout -q "$version" && ./emsdk install "$version" > /dev/null && ./emsdk activate "$version" > /dev/null)
}
install_emsdk "$EMSDK_NEW" "$HOME/emsdk"
install_emsdk "$EMSDK_OLD" "$HOME/emsdk-$EMSDK_OLD"

if ! command -v rustup > /dev/null; then
  curl -sSf https://sh.rustup.rs | sh -s -- -y --profile minimal --default-toolchain "$RUST"
fi
# shellcheck disable=SC1091
[ -f "$HOME/.cargo/env" ] && . "$HOME/.cargo/env"
rustup toolchain install "$RUST" --profile minimal > /dev/null
rustup default "$RUST" > /dev/null
rustup target add wasm32-wasip1 wasm32-unknown-emscripten > /dev/null

PY=$(command -v python3.12 || true)
[ -n "$PY" ] || { echo "python3.12 is required (apt install python3.12 python3.12-venv)" >&2; exit 1; }
mkdir -p "$TOOLS"
[ -x "$TOOLS/pb312/bin/python3" ] || "$PY" -m venv "$TOOLS/pb312"
"$TOOLS/pb312/bin/pip" install -q "pyodide-build==$PYODIDE_BUILD" "wheel<0.46"

if [ ! -d "$TOOLS/xbuildenv" ]; then
  mkdir -p "$TOOLS/xbuildenv"
  curl -sSfL "https://github.com/pyodide/pyodide/releases/download/$PYODIDE/xbuildenv-$PYODIDE.tar.bz2" \
    | tar xj -C "$TOOLS/xbuildenv"
fi
mkdir -p "$TOOLS/.pyodide-xbuildenv-$PYODIDE_BUILD"
[ -e "$TOOLS/.pyodide-xbuildenv-$PYODIDE_BUILD/xbuildenv" ] || ln -s "$TOOLS/xbuildenv" "$TOOLS/.pyodide-xbuildenv-$PYODIDE_BUILD/xbuildenv"
touch "$TOOLS/.pyodide-xbuildenv-$PYODIDE_BUILD/.installed"

echo "toolchains ready:"
echo "  emsdk $EMSDK_NEW, $EMSDK_OLD; rust $(rustc --version | cut -d' ' -f2); python $("$TOOLS/pb312/bin/python3" --version | cut -d' ' -f2); pyodide-build $PYODIDE_BUILD"
echo "  PYHMMER_TOOLS=$TOOLS"
