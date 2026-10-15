#!/usr/bin/env bash
# Assemble the self-contained GitHub Pages site: the page, the shared host and bridge, the WebAssembly
# bundle, the curated database, a MATPredict wheel built from this checkout, and a manifest. No symlinks
# (Pages does not follow them), so every asset is copied. Needs web/spike's node_modules (npm ci).
# Usage: build-site.sh <bundle-dir> <site-dir>
set -eo pipefail
BUNDLE="$1"; SITE="$2"
HERE="$(cd "$(dirname "$0")" && pwd)"
SPIKE="$(cd "$HERE/.." && pwd)"
REPO="$(cd "$SPIKE/../.." && pwd)"
rm -rf "$SITE" && mkdir -p "$SITE/page/_assets" "$SITE/browser/host" "$SITE/browser/python" \
  "$SITE/browser/pyhmmer" "$SITE/localizer" "$SITE/node_modules/@bjorn3/browser_wasi_shim"

# The page and the shared modules it imports (paths relative to the page, as in the repo).
cp "$HERE/index.html" "$HERE/app.mjs" "$SITE/page/"
cp "$SPIKE/browser/host/tool-host.mjs" "$SPIKE/browser/host/wasi-runner.mjs" "$SITE/browser/host/"
cp "$SPIKE/browser/python/bridge.py" "$SITE/browser/python/"
cp "$SPIKE/browser/pyhmmer/psutil.py" "$SITE/browser/pyhmmer/"
cp "$SPIKE/localizer/localizer.py" "$SITE/localizer/"
cp -RL "$SPIKE/node_modules/@bjorn3/browser_wasi_shim/." "$SITE/node_modules/@bjorn3/browser_wasi_shim/"

# WebAssembly tools and the pyhmmer wheel, copied (not linked) from the bundle.
# The page requests these from _assets/bundle/ (see app.mjs and index.html).
mkdir -p "$SITE/page/_assets/bundle"
cp -RL "$BUNDLE/miniprot" "$BUNDLE/mmseqs" "$BUNDLE/exonerate" "$BUNDLE/losat" "$SITE/page/_assets/bundle/"
cp -L "$BUNDLE"/pyhmmer-*.whl "$SITE/page/_assets/bundle/"

# The curated database, copied from the checkout.
cp -RL "$REPO/db" "$SITE/page/_assets/db"

# A MATPredict wheel built from this checkout, so the site runs the code it was deployed from.
# MATPredict needs Python 3.11+ to build; PYTHON overrides the default interpreter.
PYTHON="${PYTHON:-python3.12}"
command -v "$PYTHON" > /dev/null || PYTHON=python3
WHEEL_DIR="$(mktemp -d)"
"$PYTHON" -m pip wheel --no-deps -q -w "$WHEEL_DIR" "$REPO"
cp "$WHEEL_DIR"/matpredict-*.whl "$SITE/page/_assets/matpredict.whl"

# The database manifest the page mounts, and the pyhmmer wheel's name.
"$PYTHON" - "$SITE/page/_assets" <<'PY'
import json, os, sys
assets = sys.argv[1]
db = os.path.join(assets, "db")
files = sorted(os.path.relpath(os.path.join(r, n), db) for r, _, names in os.walk(db) for n in names)
wheels = sorted(n for n in os.listdir(os.path.join(assets, "bundle")) if n.startswith("pyhmmer-") and n.endswith(".whl"))
assert wheels, "no pyhmmer wheel in the bundle"
json.dump({"db": files, "pyhmmerWheel": wheels[0]}, open(os.path.join(assets, "db-manifest.json"), "w"), indent=1)
print(f"manifest: {len(files)} database files, {wheels[0]}")
PY

# Pages root redirects to the page; .nojekyll keeps underscore directories (_assets) from being dropped.
cat > "$SITE/index.html" <<'HTML'
<!doctype html><meta charset="utf-8"><meta http-equiv="refresh" content="0; url=page/index.html">
<a href="page/index.html">MATPredict browser page</a>
HTML
touch "$SITE/.nojekyll"
rm -rf "$WHEEL_DIR"
du -sh "$SITE"
