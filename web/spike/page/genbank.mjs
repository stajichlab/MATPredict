// Locus records for download: a FASTA of the locus with its flanks, and an annotated GenBank file.
// Pure functions (no DOM), so they can be tested in Node. Coordinates in the GenBank file are relative
// to the extracted window; the source feature says where the window lies on the contig.

const WRAP = 60;

// Window around a locus, 1-based inclusive, clipped to the contig.
export function locusWindow(locus, flank, contigLength) {
  return { start: Math.max(1, locus.start - flank), end: Math.min(contigLength, locus.end + flank) };
}

function chunks(text, n) {
  const out = [];
  for (let i = 0; i < text.length; i += n) out.push(text.slice(i, i + n));
  return out;
}

export function locusFasta(locus, seq, flank, win = locusWindow(locus, flank, seq.length)) {
  const header = `>${locus.contig}:${win.start}-${win.end} MATPredict ${locus.id} flank=${flank}`;
  return [header, ...chunks(seq.slice(win.start - 1, win.end), WRAP)].join("\n") + "\n";
}

const sanitize = (s) => String(s ?? "").replace(/["\n\r]/g, " ").trim();

// A GenBank qualifier, wrapped at 79 columns with the 21-column continuation indent.
function qualifier(name, value) {
  const text = `/${name}="${sanitize(value)}"`;
  const lines = [];
  let rest = text;
  const width = 58;
  lines.push(rest.slice(0, width));
  rest = rest.slice(width);
  while (rest.length) {
    lines.push(rest.slice(0, width));
    rest = rest.slice(width);
  }
  return lines.map((l) => " ".repeat(21) + l).join("\n") + "\n";
}

function feature(key, location, qualifiers) {
  // GenBank: feature key at column 6, location at column 22.
  const head = `     ${key.padEnd(16)}${location}`;
  return head + "\n" + qualifiers.map(([k, v]) => qualifier(k, v)).join("");
}

function formatDate(d) {
  const months = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"];
  const p = (n) => String(n).padStart(2, "0");
  return `${p(d.getUTCDate())}-${months[d.getUTCMonth()]}-${d.getUTCFullYear()}`;
}

// The sequence block, 60 bases per line in groups of ten, lowercase as GenBank writes it.
function originBlock(seq) {
  const lines = [];
  const lower = seq.toLowerCase();
  for (let i = 0; i < lower.length; i += WRAP) {
    const groups = chunks(lower.slice(i, i + WRAP), 10).join(" ");
    lines.push(`${String(i + 1).padStart(9)} ${groups}`);
  }
  return "ORIGIN\n" + lines.join("\n") + "\n";
}

// An annotated GenBank record for one locus: the locus, its genes present in the search, and the source.
export function locusGenbank(locus, seq, flank, { win = locusWindow(locus, flank, seq.length), date = new Date() } = {}) {
  const sub = seq.slice(win.start - 1, win.end);
  const len = sub.length;
  const rel = (p) => p - win.start + 1;
  const name = locus.contig.slice(0, 16);
  const loc = (a, b) => `${a}..${b}`;
  const out = [];
  out.push(`LOCUS       ${name.padEnd(16)} ${String(len).padStart(11)} bp    DNA     linear   UNK ${formatDate(date)}`);
  out.push(`DEFINITION  ${sanitize(`MATPredict ${locus.id} (${locus.family}); ${locus.contig}:${win.start}-${win.end}, flank ${flank} bp`)}.`);
  out.push(`ACCESSION   ${locus.contig}`);
  out.push(`VERSION     ${locus.contig}`);
  out.push("KEYWORDS    MATPredict.");
  out.push("SOURCE      .");
  out.push("  ORGANISM  .");
  out.push("FEATURES             Location/Qualifiers");
  out.push(feature("source", loc(1, len), [
    ["mol_type", "genomic DNA"],
    ["note", `extracted from ${locus.contig}:${win.start}-${win.end} (1-based); coordinates in this record are relative to the extract`],
  ]).trimEnd());
  const lstart = rel(Math.max(win.start, locus.start));
  const lend = rel(Math.min(win.end, locus.end));
  const locNote = [locus.family, locus.idiomorph && `idiomorph ${locus.idiomorph}`, locus.confidence && `confidence ${locus.confidence}`]
    .filter(Boolean).join("; ");
  out.push(feature("misc_feature", loc(lstart, lend), [
    ["label", `MAT locus ${locus.id}`],
    ["note", locNote],
  ]).trimEnd());
  for (const g of locus.genes) {
    if (g.end < win.start || g.start > win.end) continue;
    const a = rel(Math.max(win.start, g.start));
    const b = rel(Math.min(win.end, g.end));
    const location = g.strand === "-" ? `complement(${loc(a, b)})` : loc(a, b);
    const note = [g.role && `role ${g.role}`, g.identity && `identity ${g.identity}%`, g.method, g.status]
      .filter(Boolean).join("; ");
    out.push(feature("gene", location, [["gene", g.name], ["note", note]]).trimEnd());
  }
  out.push(originBlock(sub) + "//");
  return out.join("\n") + "\n";
}
