#!/usr/bin/env bash
# Download GmG's WebAssembly package (web target) at a pinned release and put its files in <dir>, for the page's "Guess phylum"
# button. The release's SHA-256 is checked, so a changed asset is refused. Needs curl and tar.
# Usage: fetch-gmg-wasm.sh <dir>      (GMG_WASM_VERSION and GMG_WASM_SHA256 override the pin, for a newer GmG release)
set -eo pipefail
DIR="$1"
V="${GMG_WASM_VERSION:-v0.3.0}"
SHA="${GMG_WASM_SHA256:-3e9a072c1c72187bee5092699f315d6db9557995968aa17ee89e0450a01a14c8}"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
curl -fsSL --retry 3 -o "$TMP/gmg-wasm.tar.gz" "https://github.com/stajichlab/GmG/releases/download/$V/gmg-wasm-$V.tar.gz"
GOT="$( (sha256sum "$TMP/gmg-wasm.tar.gz" 2>/dev/null || shasum -a 256 "$TMP/gmg-wasm.tar.gz") | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || { echo "GmG wasm $V: sha256 $GOT does not match the pinned $SHA" >&2; exit 1; }
tar -xzf "$TMP/gmg-wasm.tar.gz" -C "$TMP"
mkdir -p "$DIR"
cp "$TMP/gmg-wasm-$V/web/gmg_wasm.js" "$TMP/gmg-wasm-$V/web/gmg_wasm_bg.wasm" "$TMP/gmg-wasm-$V/LICENSE" "$DIR/"
echo "GmG wasm $V -> $DIR"
