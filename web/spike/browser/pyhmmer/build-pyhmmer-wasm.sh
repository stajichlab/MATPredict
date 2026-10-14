#!/usr/bin/env bash
# Build pyhmmer 0.12.3 (the version MATPredict pins) as a Pyodide 0.27.7 wheel. Spike only.
# Usage: build-pyhmmer-wasm.sh <work-dir>
# Result: <work-dir>/dist/pyhmmer-0.12.3-cp312-cp312-pyodide_2024_0_wasm32.whl
#
# Prerequisites (set up once):
#   - emsdk 3.1.58 at $HOME/emsdk-3.1.58 (Pyodide 0.27.7 is built with it)
#   - Python 3.12 (Homebrew python@3.12) for the build venv, with pyodide-build==0.27.3 and wheel<0.46
#   - the Pyodide 0.27.7 build environment: xbuildenv-0.27.7.tar.bz2 from the Pyodide release,
#     unpacked into <work-dir>/xbuildenv, with <work-dir>/.pyodide-xbuildenv-0.27.3/xbuildenv linking to it
#     and an empty .installed marker (the installer's package step is not needed for pyhmmer)
set -eo pipefail
WORK="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$WORK" && cd "$WORK"
# Toolchain from an existing setup (see the prerequisites above). PYHMMER_TOOLS points at it.
TOOLS="${PYHMMER_TOOLS:-$HOME/wasm-build/pyhmmer}"
[ -e pb312 ] || ln -s "$TOOLS/pb312" pb312
[ -e .pyodide-xbuildenv-0.27.3 ] || ln -s "$TOOLS/.pyodide-xbuildenv-0.27.3" .pyodide-xbuildenv-0.27.3
# Download with the Python 3.12 toolchain, not whatever python3 is first on PATH (the system 3.9 pip
# stalled on retries here). An already-present sdist is reused.
PY_DL="$TOOLS/pb312/bin/python3"; [ -x "$PY_DL" ] || PY_DL=python3
[ -f pyhmmer-0.12.3.tar.gz ] || "$PY_DL" -m pip download --no-binary :all: --no-deps "pyhmmer==0.12.3" -d . -q
rm -rf pyhmmer-0.12.3 && tar xzf pyhmmer-0.12.3.tar.gz
patch -p0 -i "$HERE/pyhmmer-pyodide.patch" > /dev/null
source "$HOME/emsdk-3.1.58/emsdk_env.sh" > /dev/null 2>&1
PATH="$WORK/pb312/bin:$PATH" CFLAGS="-O2 -I$HERE" LDFLAGS="-sWASM_BIGINT" \
  "$WORK/pb312/bin/pyodide" build pyhmmer-0.12.3 --outdir "$WORK/dist"
ls -la "$WORK/dist"
