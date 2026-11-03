// Click-to-sort for plain HTML tables: click a header to sort by that column, click it again to reverse.
// Numbers sort as numbers ("1,234", "88%", "9/9" by its first number), everything else as text with run
// numbers in natural order (SRR2 before SRR10), and an empty cell or a dash sorts last in either direction.
// The pure part (sortRows) has no DOM and is tested in sortable.test.mjs.

const MISSING = /^(|[–—-]|n\/a|none|not run)$/i;
const NUMBER = /^[-+]?(\d+\.?\d*|\.\d+)\s*(%|\/\s*\d+.*)?$/;

export const isMissing = (t) => MISSING.test(String(t).trim());
export const cellNumber = (t) => {
  const s = String(t).trim().replace(/,/g, "");
  return NUMBER.test(s) ? parseFloat(s) : NaN;
};

// rows: [{cells: [text, ...], ...}]. Returns a new array sorted by column `col`, ascending for dir = 1.
export function sortRows(rows, col, dir = 1) {
  const present = rows.filter((r) => !isMissing(r.cells[col] ?? ""));
  const numeric = present.length > 0 && present.every((r) => !Number.isNaN(cellNumber(r.cells[col])));
  const cmp = numeric
    ? (a, b) => cellNumber(a.cells[col]) - cellNumber(b.cells[col])
    : (a, b) => String(a.cells[col]).localeCompare(String(b.cells[col]), undefined, { numeric: true, sensitivity: "base" });
  const sorted = [...present].sort((a, b) => cmp(a, b) * dir);   // Array.prototype.sort is stable
  return [...sorted, ...rows.filter((r) => isMissing(r.cells[col] ?? ""))];
}

// Make `table` sortable. Returns {refresh}: call it after the table's rows are rebuilt (a filter, a regroup) to
// keep the chosen order and the header marks. opts.fixedLast(tr) marks rows that stay at the bottom (a total row).
export function makeSortable(table, { fixedLast = null } = {}) {
  const state = { col: -1, dir: 1 };
  table.classList.add("sortable");
  const headerCells = () => [...table.querySelectorAll("tr")].filter((tr) => tr.querySelector("th")).flatMap((tr) => [...tr.children].filter((c) => c.tagName === "TH"));
  const bodyRows = () => [...table.querySelectorAll("tr")].filter((tr) => !tr.querySelector("th"));
  function mark() {
    headerCells().forEach((th) => {
      th.tabIndex = 0;
      if (!th.title) th.title = "Click to sort; click again to reverse";
      if (th.cellIndex === state.col) th.setAttribute("aria-sort", state.dir === 1 ? "ascending" : "descending"); else th.removeAttribute("aria-sort");
    });
  }
  function apply() {
    if (state.col >= 0) {
      const all = bodyRows();
      const fixed = fixedLast ? all.filter(fixedLast) : [];
      const rows = all.filter((tr) => !fixed.includes(tr)).map((el) => ({ el, cells: [...el.children].map((c) => c.innerText ?? c.textContent) }));
      const parent = all.length ? all[0].parentNode : null;
      if (parent) [...sortRows(rows, state.col, state.dir).map((r) => r.el), ...fixed].forEach((tr) => parent.appendChild(tr));
    }
    mark();
  }
  const choose = (th) => { state.dir = th.cellIndex === state.col ? -state.dir : 1; state.col = th.cellIndex; apply(); };
  table.addEventListener("click", (e) => { const th = e.target.closest("th"); if (th && table.contains(th)) choose(th); });
  table.addEventListener("keydown", (e) => { const th = e.target.closest && e.target.closest("th"); if (th && (e.key === "Enter" || e.key === " ")) { e.preventDefault(); choose(th); } });
  return { refresh: apply };
}
