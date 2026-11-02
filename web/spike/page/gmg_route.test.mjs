// node --test web/spike/page/gmg_route.test.mjs
// The same cases as tests/detect/test_gmg_guess.py (decision table and record reader), so the page and `detect --phylum auto` agree.
import test from "node:test";
import assert from "node:assert/strict";
import { guessFromRecord, decide } from "./gmg_route.mjs";
import { initGuess } from "./gmg_guess.mjs";

const AVAILABLE = ["Ascomycota", "Basidiomycota", "Mucoromycota"];
const rec = (phylum, flags = [], status = "ok", hint) => ({
  status, qc_flags: flags, guess: { calls: { phylum: { taxon: phylum, taxid: phylum ? 4890 : null } } }, ...(hint ? { routing_hint: hint } : {}),
});

test("a record without a routing_hint gets the same rule", () => {
  const cases = [
    [rec("Ascomycota"), null], [rec("Ascomycota", ["weak_phylum"], "low_confidence"), "weak_phylum"],
    [rec("Chordata", ["not_fungal"], "no_match"), "not_fungal"], [rec(null, [], "no_match"), "no_phylum_call"], [rec(null, [], "error"), "error"],
  ];
  for (const [r, reason] of cases) {
    const g = guessFromRecord(r);
    assert.equal(g.reason, reason);
    assert.equal(g.usable, reason === null);
  }
});

test("the routing_hint of a schema 9 record is used as given", () => {
  const g = guessFromRecord(rec("Ascomycota", [], "ok", { usable: true, reason: null, phylum: "Ascomycota" }));
  assert.deepEqual([g.usable, g.phylum], [true, "Ascomycota"]);
});

test("a usable call in a curated phylum is selected", () => {
  assert.equal(decide({ usable: true, phylum: "Mucoromycota", reason: null }, AVAILABLE).phylum, "Mucoromycota");
});

test("anything doubtful selects nothing and says why", () => {
  const cases = [
    [{ reason: "error" }, /could not read/], [{ reason: "no_phylum_call" }, /no phylum call/], [{ reason: "not_fungal", phylum: "Chordata" }, /not a fungus/],
    [{ reason: "weak_phylum", phylum: "Ascomycota" }, /weak call/], [{ usable: true, phylum: "Chytridiomycota", reason: null }, /no curated MAT locus family/],
    [{ reason: "something_new", phylum: "Ascomycota" }, /not usable/],
  ];
  for (const [g, re] of cases) {
    const d = decide(g, AVAILABLE);
    assert.equal(d.phylum, null);
    assert.match(d.note, re);
  }
});

test("a weak call is used only when accepted, and never for a non-fungal or uncurated phylum", () => {
  const weak = { reason: "weak_phylum", phylum: "Ascomycota" };
  assert.equal(decide(weak, AVAILABLE).phylum, null);
  const d = decide(weak, AVAILABLE, true);
  assert.deepEqual([d.phylum, d.weak], ["Ascomycota", true]);
  assert.equal(decide({ reason: "not_fungal", phylum: "Chordata" }, AVAILABLE, true).phylum, null);
  assert.equal(decide({ reason: "weak_phylum", phylum: "Chytridiomycota" }, AVAILABLE, true).phylum, null);
});

// The button: a fake DOM and a fake runner.
function fakeUi(options) {
  const listeners = {};
  const el = (extra = {}) => ({ disabled: false, textContent: "", checked: false, addEventListener: (t, f) => { listeners[t] = f; }, ...extra });
  const select = { value: "Ascomycota", options: options.map((o) => ({ value: o, textContent: o })) };
  return { button: el({ disabled: true }), weak: el(), note: el(), select, click: () => listeners.click() };
}

test("the Guess button fills the pulldown, and leaves it alone when the guess is not usable", async () => {
  let genome = null;
  const logs = [];
  const ui = fakeUi(AVAILABLE);
  let next = { record: rec("Basidiomycota"), version: "gmg 0.3.0" };
  initGuess({ ui, getGenome: () => genome, log: (m) => logs.push(m), run: async (b, n, onProgress) => { onProgress({ total: 10, loaded: 5 }); return next; } });
  assert.equal(ui.button.disabled, true, "disabled until a genome is loaded");
  genome = { name: "g.fna", bytes: new Uint8Array(4) };
  const ctl = initGuess({ ui, getGenome: () => genome, log: (m) => logs.push(m), run: async () => next });
  ctl.update();
  assert.equal(ui.button.disabled, false);
  await ui.click();
  assert.equal(ui.select.value, "Basidiomycota");
  assert.match(ui.note.textContent, /Selected below/);
  next = { record: rec("Ascomycota", ["weak_phylum"], "low_confidence"), version: "x" };
  await ui.click();
  assert.equal(ui.select.value, "Basidiomycota", "a weak guess does not change the pulldown");
  assert.match(ui.note.textContent, /not changed/);
  ui.weak.checked = true;
  await ui.click();
  assert.equal(ui.select.value, "Ascomycota");
  assert.ok(logs.some((l) => l.includes("GmG gmg 0.3.0")));
});

test("a failed guess is reported and the pulldown is untouched", async () => {
  const ui = fakeUi(AVAILABLE);
  const logs = [];
  initGuess({ ui, getGenome: () => ({ name: "g", bytes: new Uint8Array(1) }), log: (m) => logs.push(m), run: async () => { throw new Error("offline"); } });
  await ui.click();
  assert.equal(ui.select.value, "Ascomycota");
  assert.match(ui.note.textContent, /could not guess: offline/);
  assert.equal(ui.button.disabled, false);
});
