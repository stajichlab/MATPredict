// The "Guess" button: ask GmG (in a worker) which phylum the loaded genome is in, then select it in the phylum pulldown.
// Nothing is uploaded: the genome is sketched in this browser; only the public GmG index is downloaded (once, then cached).
import { guessFromRecord, decide } from "./gmg_route.mjs";

// Pinned, like the Pyodide version: the index MATPredict is validated with (see docs/notes/2026-11-02_gmg-phylum-routing.md).
// The GmG site serves it with CORS open. build-site.sh vendors the matching wasm package under _assets/gmg/.
export const GMG_INDEX_URL = "https://gmg.fungalgenomes.org/index/bfd-reps-2026.10.5-s20000.gmgi";
export const GMG_INDEX_BYTES = 66798702;

// Run GmG on `bytes` (plain FASTA). onProgress({loaded,total,cached,stage}). Resolves to { record, version }.
export function runGuess(bytes, name, onProgress = () => {}, workerUrl = new URL("./gmg_worker.mjs", import.meta.url)) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(workerUrl, { type: "module" });
    const done = (fn, v) => { worker.terminate(); fn(v); };   // terminate: frees the index's memory before detection
    worker.onmessage = (e) => {
      const m = e.data;
      if (m.type === "record") done(resolve, { record: m.record, version: m.version });
      else if (m.type === "error") done(reject, new Error(m.message));
      else onProgress(m);
    };
    worker.onerror = (e) => done(reject, new Error(e.message || "the GmG worker failed to start"));
    const copy = bytes.slice();   // the page keeps its own bytes for the detection run
    worker.postMessage({ bytes: copy, name, indexUrl: GMG_INDEX_URL, size: GMG_INDEX_BYTES }, [copy.buffer]);
  });
}

// Wire the controls: `ui` = { button, weak (checkbox), note (element), select (the phylum <select>) }; getGenome() -> {name, bytes}|null;
// log(msg) writes to the page log. Returns { update() } to call when the genome changes (enables/disables the button).
export function initGuess({ ui, getGenome, log, run = runGuess }) {
  const available = () => [...ui.select.options].map((o) => o.value || o.textContent);
  const update = () => { ui.button.disabled = !getGenome(); };
  ui.button.addEventListener("click", async () => {
    const genome = getGenome();
    if (!genome) return;
    ui.button.disabled = true;
    ui.note.textContent = "Starting GmG…";
    const t = performance.now();
    try {
      const { record, version } = await run(genome.bytes, genome.name, (m) => {
        if (m.stage === "classify") ui.note.textContent = "Comparing the genome with the reference…";
        else if (m.stage === "index") ui.note.textContent = "Reading the index…";
        else if (m.cached) ui.note.textContent = "Loading the GmG index from this browser's cache…";
        else if (m.total) ui.note.textContent = `Downloading the GmG index (one time): ${(m.loaded / 1e6).toFixed(0)} of ${(m.total / 1e6).toFixed(0)} MB`;
      });
      const guess = guessFromRecord(record);
      const d = decide(guess, available(), ui.weak.checked);
      if (d.phylum) ui.select.value = d.phylum;
      ui.note.textContent = d.phylum ? `${d.note} Selected below.` : `${d.note} The phylum was not changed.`;
      log(`GmG${version ? " " + version : ""}: ${d.note} (${((performance.now() - t) / 1000).toFixed(1)} s)`);
    } catch (err) {
      ui.note.textContent = `GmG could not guess: ${err.message}. Choose the phylum yourself.`;
      log(`GmG guess failed: ${err.message}`);
    } finally {
      update();
    }
  });
  update();
  return { update };
}
