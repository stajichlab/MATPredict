// node --test web/intro/reads/grouping.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import { UNPLACED, allKeys, groupSpecies, lineageOf, rowsOf, summarize } from "./grouping.js";

const sp = (species, phylum, cls, order, confidence = "confirmed", runs = 10) => ({ species, phylum, class: cls, order, confidence, runs });
const ROWS = [
  sp("Aspergillus flavus", "Ascomycota", "Eurotiomycetes", "Eurotiales"),
  sp("Aspergillus niger", "Ascomycota", "Eurotiomycetes", "Eurotiales", "withheld", 5),
  sp("Coccidioides immitis", "Ascomycota", "Eurotiomycetes", "Onygenales", "provisional", 40),
  sp("Fusarium graminearum", "Ascomycota", "Sordariomycetes", "Hypocreales"),
  sp("Candida albicans", "Ascomycota", "Saccharomycetes", "Saccharomycetales", "provisional (homothallic)", 20),
  sp("Mystery fungus", "", "", ""),
];

test("a species is grouped under its phylum, class and order, and a blank rank gets a placeholder", () => {
  assert.deepEqual(lineageOf(ROWS[0]), ["Ascomycota", "Eurotiomycetes", "Eurotiales"]);
  assert.deepEqual(lineageOf(ROWS[5]), [UNPLACED, "(class not assigned)", "(order not assigned)"]);
  assert.deepEqual(lineageOf(sp("x", "Ascomycota", "", "Eurotiales")), ["Ascomycota", "(class not assigned)", "Eurotiales"]);
});

test("groups come out in name order with the unplaced ones last", () => {
  const tree = groupSpecies(ROWS);
  assert.deepEqual(tree.map((g) => g.name), ["Ascomycota", UNPLACED]);
  const asco = tree[0];
  assert.deepEqual(asco.children.map((c) => c.name), ["Eurotiomycetes", "Saccharomycetes", "Sordariomycetes"]);
  assert.deepEqual(asco.children[0].children.map((o) => o.name), ["Eurotiales", "Onygenales"]);
  assert.deepEqual(tree[1].children.map((c) => c.name), ["(class not assigned)"]);
});

test("a placeholder class sorts after the named classes of the same phylum", () => {
  const tree = groupSpecies([sp("a", "Ascomycota", "", "Eurotiales"), sp("b", "Ascomycota", "Zeta", "O"), sp("c", "Ascomycota", "Alpha", "O")]);
  assert.deepEqual(tree[0].children.map((c) => c.name), ["Alpha", "Zeta", "(class not assigned)"]);
});

test("the rows of an order keep the caller's order, so a sort applies inside every group", () => {
  const reversed = [...ROWS].reverse();
  const eurotiales = groupSpecies(reversed)[0].children[0].children[0];
  assert.deepEqual(eurotiales.rows.map((r) => r.species), ["Aspergillus niger", "Aspergillus flavus"]);
});

test("keys are paths, levels count down the tree, and every row is in exactly one leaf", () => {
  const tree = groupSpecies(ROWS);
  assert.deepEqual([tree[0].level, tree[0].children[0].level, tree[0].children[0].children[0].level], [0, 1, 2]);
  assert.equal(tree[0].children[0].children[0].key, "Ascomycota/Eurotiomycetes/Eurotiales");
  assert.equal(rowsOf(tree[0]).length + rowsOf(tree[1]).length, ROWS.length);
  const keys = allKeys(tree);
  assert.equal(new Set(keys).size, keys.length);
  assert.ok(keys.includes("Ascomycota") && keys.includes("Ascomycota/Sordariomycetes/Hypocreales"));
});

test("summarize counts species, runs and the confidence kinds, with both provisional kinds together", () => {
  assert.deepEqual(summarize(ROWS), { species: 6, runs: 95, confirmed: 3, provisional: 2, withheld: 1 });
  assert.deepEqual(summarize([]), { species: 0, runs: 0, confirmed: 0, provisional: 0, withheld: 0 });
});

test("filtering rows first leaves the keys of the groups that remain unchanged", () => {
  const all = groupSpecies(ROWS), some = groupSpecies(ROWS.filter((r) => r.species.startsWith("Aspergillus")));
  assert.deepEqual(allKeys(some), ["Ascomycota", "Ascomycota/Eurotiomycetes", "Ascomycota/Eurotiomycetes/Eurotiales"]);
  assert.ok(allKeys(some).every((k) => allKeys(all).includes(k)));
});
