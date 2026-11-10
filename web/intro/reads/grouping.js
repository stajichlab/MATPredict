// Group the species of the raw-reads overview into phylum > class > order, for collapsible sections.
// Pure functions, no DOM: tested in grouping.test.mjs.

export const UNPLACED = "Unplaced";
const NO_CLASS = "(class not assigned)", NO_ORDER = "(order not assigned)";

// The three names a species is grouped under; blank ranks get a placeholder that sorts last.
export const lineageOf = (s) => [s.phylum || UNPLACED, s.class || NO_CLASS, s.order || NO_ORDER];

const placeholder = (name) => name === UNPLACED || name.startsWith("(");
const byName = (a, b) => (placeholder(a.name) - placeholder(b.name)) || a.name.localeCompare(b.name);

// Counts for a group's summary line. `confidence` values are the page's: confirmed, provisional,
// "provisional (homothallic)", unvalidated and withheld; the two provisional kinds are counted together.
export function summarize(rows) {
  const t = { species: rows.length, runs: 0, confirmed: 0, provisional: 0, unvalidated: 0, not_typable: 0, withheld: 0 };
  for (const r of rows) {
    t.runs += r.runs || 0;
    if (r.confidence === "confirmed") t.confirmed++;
    else if (r.confidence === "withheld") t.withheld++;
    else if (r.confidence === "unvalidated") t.unvalidated++;
    else if (r.confidence === "not typable by presence") t.not_typable++;
    else t.provisional++;
  }
  return t;
}

// rows -> [{key, name, level: 0, children: [{key, name, level: 1, children: [{key, name, level: 2, rows}]}]}].
// Groups are in name order with the placeholders last; the rows of a leaf keep the order they came in, so the
// caller's sort applies inside every group. `key` is a path ("Ascomycota/Eurotiomycetes/Eurotiales") that stays the
// same when a filter hides other groups, so an open or closed section is remembered.
export function groupSpecies(rows) {
  const root = new Map();
  for (const r of rows) {
    const [p, c, o] = lineageOf(r);
    if (!root.has(p)) root.set(p, new Map());
    const classes = root.get(p);
    if (!classes.has(c)) classes.set(c, new Map());
    const orders = classes.get(c);
    if (!orders.has(o)) orders.set(o, []);
    orders.get(o).push(r);
  }
  const named = (map, make) => [...map].map(([name, v]) => make(name, v)).sort(byName);
  return named(root, (p, classes) => ({
    key: p, name: p, level: 0,
    children: named(classes, (c, orders) => ({
      key: `${p}/${c}`, name: c, level: 1,
      children: named(orders, (o, rs) => ({ key: `${p}/${c}/${o}`, name: o, level: 2, rows: rs })),
    })),
  }));
}

// Every row under a group, in group order.
export const rowsOf = (g) => (g.rows ? g.rows : g.children.flatMap(rowsOf));

// Every group key of a tree, for expand all and collapse all.
export const allKeys = (groups) => groups.flatMap((g) => [g.key, ...(g.children ? allKeys(g.children) : [])]);
