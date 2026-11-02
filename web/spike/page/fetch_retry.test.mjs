// node --test web/spike/page/fetch_retry.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import { fetchBytes, FetchError, retryAfterMs } from "./fetch_retry.mjs";

const ok = (text) => ({ ok: true, status: 200, headers: new Headers(), arrayBuffer: async () => new TextEncoder().encode(text).buffer });
const bad = (status, headers = {}) => ({ ok: false, status, headers: new Headers(headers), arrayBuffer: async () => new ArrayBuffer(0) });

// A fetch that answers from a script and records how often it was called, and a sleep that only records.
function harness(script) {
  const calls = [];
  const waits = [];
  const fetchFn = async (url) => {
    calls.push(url);
    const next = script.shift();
    if (next instanceof Error) throw next;
    return next;
  };
  return { calls, waits, opts: { fetchFn, sleep: async (ms) => { waits.push(ms); }, random: () => 0 } };
}

test("a 503 is retried and the later 200 is returned", async () => {
  const h = harness([bad(503), bad(503), ok("wheel")]);
  const bytes = await fetchBytes("./_assets/matpredict.whl", h.opts);
  assert.equal(new TextDecoder().decode(bytes), "wheel");
  assert.equal(h.calls.length, 3);
  assert.deepEqual(h.waits, [500, 1500]);
});

test("the wait follows Retry-After when the server gives one, capped at 10 s", async () => {
  const h = harness([bad(429, { "retry-after": "2" }), bad(503, { "retry-after": "999" }), ok("x")]);
  await fetchBytes("a/b.bin", h.opts);
  assert.deepEqual(h.waits, [2000, 10000]);
});

test("a 404 fails at once, without a retry", async () => {
  const h = harness([bad(404), ok("never")]);
  await assert.rejects(fetchBytes("./_assets/nope.json", h.opts), (e) => e instanceof FetchError && e.status === 404 && e.attempts === 1);
  assert.equal(h.calls.length, 1);
  assert.equal(h.waits.length, 0);
});

test("a persistent 503 gives up after the retries, with a message that names the file and the status", async () => {
  const h = harness([bad(503), bad(503), bad(503), bad(503), ok("never")]);
  await assert.rejects(fetchBytes("./_assets/matpredict.whl", { ...h.opts, retries: 3 }), (e) => {
    assert.ok(e instanceof FetchError);
    assert.equal(e.message, "could not fetch matpredict.whl: HTTP 503, tried 4 times");
    return true;
  });
  assert.equal(h.calls.length, 4);
});

test("a network error is retried like a 503", async () => {
  const h = harness([new TypeError("Failed to fetch"), ok("ok")]);
  const bytes = await fetchBytes("x.txt", h.opts);
  assert.equal(bytes.length, 2);
  assert.equal(h.calls.length, 2);
});

test("onRetry is told the file, the status, the attempt and the wait", async () => {
  const h = harness([bad(502), ok("ok")]);
  const seen = [];
  await fetchBytes("dir/f.faa", { ...h.opts, onRetry: (...a) => seen.push(a) });
  assert.deepEqual(seen, [["dir/f.faa", 502, 1, 500]]);
});

test("retryAfterMs ignores a missing or non-numeric header", () => {
  assert.equal(retryAfterMs(new Headers()), null);
  assert.equal(retryAfterMs(new Headers({ "retry-after": "Wed, 21 Oct 2026 07:28:00 GMT" })), null);
  assert.equal(retryAfterMs(new Headers({ "retry-after": "0" })), 0);
});
