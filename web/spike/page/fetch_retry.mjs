// Fetch a file as bytes, retrying the failures that pass on their own. A run requests several hundred
// files in a few seconds, and GitHub Pages can answer one of them with a 429 or 503; without a retry
// that single answer ended the whole run. A 404 and the like fail at once.

export const RETRYABLE_STATUS = new Set([408, 425, 429, 500, 502, 503, 504]);
const MAX_WAIT_MS = 10000;

export class FetchError extends Error {
  constructor(url, status, attempts, detail) {
    const name = String(url).split("/").pop() || String(url);
    const what = status ? `HTTP ${status}` : `network error (${detail})`;
    const tries = attempts > 1 ? `, tried ${attempts} times` : "";
    super(`could not fetch ${name}: ${what}${tries}`);
    this.name = "FetchError";
    this.url = url;
    this.status = status;
    this.attempts = attempts;
  }
}

// Seconds in a Retry-After header, in milliseconds; null when absent or not a number of seconds.
export function retryAfterMs(headers) {
  const v = headers && headers.get ? headers.get("retry-after") : null;
  const s = v === null || v === "" ? NaN : Number(v);
  return Number.isFinite(s) && s >= 0 ? Math.min(s * 1000, MAX_WAIT_MS) : null;
}

// `retries` is the number of extra tries after the first. Waits baseMs, 3x baseMs, 9x baseMs, ... (plus up
// to 25% jitter so parallel fetches do not retry in step), or the server's Retry-After when it gives one.
// fetchFn, sleep and random are injectable for tests; onRetry(url, status, attempt, waitMs) lets the page log.
export async function fetchBytes(url, {
  retries = 3, baseMs = 500, fetchFn = (u) => fetch(u), sleep = (ms) => new Promise((r) => setTimeout(r, ms)),
  random = Math.random, onRetry = null,
} = {}) {
  for (let attempt = 1; ; attempt++) {
    let status = 0;
    let detail = "";
    let headers = null;
    try {
      const r = await fetchFn(url);
      if (r.ok) return new Uint8Array(await r.arrayBuffer());
      status = r.status;
      headers = r.headers;
    } catch (err) {
      detail = err && err.message ? err.message : String(err);
    }
    const retryable = status === 0 || RETRYABLE_STATUS.has(status);
    if (!retryable || attempt > retries) throw new FetchError(url, status, attempt, detail);
    const wait = retryAfterMs(headers) ?? Math.round(baseMs * 3 ** (attempt - 1) * (1 + 0.25 * random()));
    if (onRetry) onRetry(url, status, attempt, wait);
    await sleep(wait);
  }
}
