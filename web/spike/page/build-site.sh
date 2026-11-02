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
cp "$HERE/index.html" "$HERE/app.mjs" "$HERE/genbank.mjs" "$HERE/genbank_input.mjs" "$HERE/proteins.mjs" "$HERE/fetch_retry.mjs" \
  "$HERE/gmg_guess.mjs" "$HERE/gmg_route.mjs" "$HERE/gmg_worker.mjs" "$SITE/page/"
mkdir -p "$SITE/page/examples" && cp "$HERE/examples/examples.json" "$SITE/page/examples/"
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

# GmG's WebAssembly package (the phylum guess), pinned and checked against its published SHA-256.
bash "$HERE/fetch-gmg-wasm.sh" "$SITE/page/_assets/gmg"

# The curated database, copied from the checkout.
cp -RL "$REPO/db" "$SITE/page/_assets/db"

# A MATPredict wheel built from this checkout, so the site runs the code it was deployed from.
# MATPredict needs Python 3.11+ to build; PYTHON overrides the default interpreter.
PYTHON="${PYTHON:-python3.12}"
command -v "$PYTHON" > /dev/null || PYTHON=python3
WHEEL_DIR="$(mktemp -d)"
"$PYTHON" -m pip wheel --no-deps -q -w "$WHEEL_DIR" "$REPO"
cp "$WHEEL_DIR"/matpredict-*.whl "$SITE/page/_assets/matpredict.whl"

# The database as one archive (one request per run, not one per file), and the manifest the page mounts:
# the file list, the archive's size and SHA-256, and the pyhmmer wheel's name.
"$PYTHON" "$HERE/build_db_manifest.py" "$SITE/page/_assets"
rm -rf "$SITE/page/_assets/db"   # the page reads the archive; the loose files are not served

# .nojekyll keeps underscore directories (_assets) from being dropped by Pages.
# The site root is the landing page (web/intro); it links to the tool at page/index.html.
# The landing page and results explorer live in web/intro; their data is regenerated from results/ first.
INTRO="$REPO/web/intro"
"$PYTHON" "$REPO/scripts/build_site_data.py" --out "$INTRO/data"
"$PYTHON" "$REPO/scripts/build_reads_site.py" --repo "$REPO" --out "$INTRO/data/reads"
cp "$INTRO/landing.html" "$SITE/index.html"
cp "$INTRO/results.html" "$INTRO/site.css" "$INTRO/charts.js" "$INTRO/theme.js" "$INTRO/locus-view.js" "$INTRO/about.html" "$INTRO/about.md" "$INTRO/md.js" "$SITE/"
cp -R "$INTRO/data" "$SITE/data"
cp -R "$INTRO/logo" "$SITE/logo"
cp -R "$INTRO/reads" "$SITE/reads"
touch "$SITE/.nojekyll"
rm -rf "$WHEEL_DIR"
du -sh "$SITE"
