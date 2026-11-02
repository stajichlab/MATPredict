// Web Worker: loads the GmG index and wasm, sketches one genome and returns GmG's record. The page starts a fresh worker
// for each guess and terminates it afterwards, so the index (about 180 MB in memory) is freed before detection starts;
// the downloaded index stays in the Cache API, so the next guess does not download it again.
// In:  { bytes: Uint8Array (plain FASTA), name, indexUrl, size }
// Out: { type: "progress", loaded, total, cached } | { type: "stage", stage } | { type: "record", record, version } | { type: "error", message }
// A namespace import: `version` exists only in GmG builds after v0.2.1, and a missing named export would stop the worker starting.
import * as gmg from "./_assets/gmg/gmg_wasm.js";

const CACHE = "gmg-index-v1";
const post = (m) => self.postMessage(m);

async function fetchIndex(url, size) {
  let cache = null;
  try { cache = await caches.open(CACHE); } catch { /* no Cache API (private mode): download every time */ }
  const hit = cache && await cache.match(url);
  if (hit) {
    post({ type: "progress", loaded: 1, total: 1, cached: true });
    return new Uint8Array(await hit.arrayBuffer());
  }
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`could not fetch the GmG index (HTTP ${resp.status})`);
  // content-length counts compressed bytes when the server gzips; the reader yields decompressed ones, so then use the known size
  const total = resp.headers.get("content-encoding") ? size || 0 : Number(resp.headers.get("content-length")) || size || 0;
  const copy = resp.clone();
  const reader = resp.body.getReader();
  const chunks = [];
  let loaded = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    chunks.push(value);
    loaded += value.length;
    post({ type: "progress", loaded, total });
  }
  if (cache) cache.put(url, copy).catch(() => {}); // a failure (quota) only costs the cache
  const bytes = new Uint8Array(loaded);
  let off = 0;
  for (const c of chunks) { bytes.set(c, off); off += c.length; }
  return bytes;
}

self.onmessage = async (e) => {
  const { bytes, name, indexUrl, size } = e.data;
  try {
    await gmg.default();
    const indexBytes = await fetchIndex(indexUrl, size);
    post({ type: "stage", stage: "index" });
    const ix = new gmg.Index(indexBytes);
    post({ type: "stage", stage: "classify" });
    const record = JSON.parse(ix.classify("genome", name, bytes, "assembly", 0, 0n, 5));
    ix.free();
    post({ type: "record", record, version: gmg.version ? gmg.version() : null });
  } catch (err) {
    post({ type: "error", message: String((err && err.message) || err) });
  }
};
