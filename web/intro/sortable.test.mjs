// node --test web/intro/sortable.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import { cellNumber, isMissing, sortRows } from "./sortable.js";

const rows = (...col) => col.map((c, i) => ({ cells: [c], id: i }));
const order = (rs, col = 0, dir = 1) => sortRows(rs, col, dir).map((r) => r.cells[col]);

test("numbers sort as numbers, not as text", () => {
  assert.deepEqual(order(rows("9", "10", "2", "1,234")), ["2", "9", "10", "1,234"]);
  assert.deepEqual(order(rows("9", "10", "2"), 0, -1), ["10", "9", "2"]);
});

test("percentages and ratios sort by their first number", () => {
  assert.deepEqual(order(rows("88%", "9%", "100%")), ["9%", "88%", "100%"]);
  assert.deepEqual(order(rows("45/46", "9/9", "100/101")), ["9/9", "45/46", "100/101"]);
});

test("run accessions sort in natural order, text case-insensitively", () => {
  assert.deepEqual(order(rows("SRR10", "SRR2", "SRR1")), ["SRR1", "SRR2", "SRR10"]);
  assert.deepEqual(order(rows("beta", "Alpha", "alpha2")), ["Alpha", "alpha2", "beta"]);
});

test("empty cells, dashes and 'not run' sort last in either direction", () => {
  const r = rows("5", "", "1", "–", "not run", "3");
  assert.deepEqual(order(r), ["1", "3", "5", "", "–", "not run"]);
  assert.deepEqual(order(r, 0, -1), ["5", "3", "1", "", "–", "not run"]);
});

test("one non-numeric cell makes the column sort as text, still in natural order", () => {
  assert.deepEqual(order(rows("10", "9", "1e")), ["1e", "9", "10"]);   // as numbers "1e" would not parse and sort wrongly
});

test("equal values keep their order (a stable sort)", () => {
  const r = [{ cells: ["a", "1"], id: 1 }, { cells: ["a", "2"], id: 2 }, { cells: ["b", "3"], id: 3 }];
  assert.deepEqual(sortRows(r, 0, 1).map((x) => x.id), [1, 2, 3]);
  assert.deepEqual(sortRows(r, 0, -1).map((x) => x.id), [3, 1, 2]);
});

test("cellNumber and isMissing", () => {
  assert.equal(cellNumber("1,234"), 1234);
  assert.equal(cellNumber("3.5"), 3.5);
  assert.ok(Number.isNaN(cellNumber("MAT1-1 4")));
  assert.ok(isMissing(" – ") && isMissing("") && isMissing("None") && !isMissing("0"));
});

test("a column missing from a row counts as missing", () => {
  const r = [{ cells: ["x", "2"] }, { cells: ["y"] }, { cells: ["z", "1"] }];
  assert.deepEqual(sortRows(r, 1, 1).map((x) => x.cells[0]), ["z", "x", "y"]);
});
