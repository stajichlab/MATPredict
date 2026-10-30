#!/usr/bin/env bash
# Assemble the browser page's assets under page/_assets: links to the WebAssembly bundle and the curated
# database, a copy of the MATPredict wheel, and a manifest of the database files (the page cannot list
# directories over HTTP). Spike only. Usage: assemble.sh <bundle-dir> [db-dir] [matpredict-wheel]
set -eo pipefail
BUNDLE="$1"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
DB="${2:-$REPO/db}"
WHEEL="${3:-$REPO/web/spike/dist/matpredict-0.8.0-py3-none-any.whl}"
A="$HERE/_assets"
mkdir -p "$A"
ln -sfn "$(cd "$BUNDLE" && pwd)" "$A/bundle"
ln -sfn "$(cd "$DB" && pwd)" "$A/db"
cp "$WHEEL" "$A/matpredict.whl"
python3 - "$DB" "$A/db-manifest.json" "$BUNDLE" <<'PY'
import json, os, sys
db, out, bundle = sys.argv[1], sys.argv[2], sys.argv[3]
files = []
for root, _, names in os.walk(db):
    for n in names:
        files.append(os.path.relpath(os.path.join(root, n), db))
wheels = sorted(n for n in os.listdir(bundle) if n.endswith(".whl"))
json.dump({"db": sorted(files), "pyhmmerWheel": wheels[0] if wheels else None}, open(out, "w"), indent=1)
print(f"manifest: {len(files)} database files, pyhmmer wheel {wheels[0] if wheels else 'MISSING'}")
PY
ls -la "$A"
