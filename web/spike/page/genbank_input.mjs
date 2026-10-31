// A GenBank file as the page's genome input: its sequence becomes the genome FASTA that the search reads.
// Contig names are the VERSION accession (what Biopython's record.id gives the command line), so loci are
// reported on the same names. The annotation is not used in the browser: the protein route needs DIAMOND,
// which has no WebAssembly build, so the genome itself is searched. Pure functions, testable in Node.

export const GENBANK_SUFFIX = /\.(gb|gbk|gbff|genbank)(\.gz)?$/i;

// True for a GenBank suffix, or text whose first non-blank line starts with LOCUS.
export function isGenbank(name, text) {
  if (GENBANK_SUFFIX.test(name || "")) return true;
  const m = /^\s*(\S+)/.exec(text.slice(0, 2000));
  return !!m && m[1] === "LOCUS";
}

// -> { fasta, contigs: [{ id, length }], cds, skipped } ; throws when no record carries a sequence.
export function genbankToFasta(text) {
  const contigs = [];
  const out = [];
  let cds = 0;
  let skipped = 0;
  for (const rec of text.split(/^\/\/\s*$/m)) {
    if (!/^\s*LOCUS/m.test(rec)) continue;
    const locus = /^LOCUS\s+(\S+)/m.exec(rec);
    const version = /^VERSION\s+(\S+)/m.exec(rec);
    const accession = /^ACCESSION\s+(\S+)/m.exec(rec);
    const id = (version || accession || locus)[1];
    cds += (rec.match(/^ {5}CDS {9}/gm) || []).length;
    const origin = /^ORIGIN.*$/m.exec(rec);
    if (!origin) { skipped++; continue; }
    const seq = rec.slice(origin.index + origin[0].length).replace(/[^A-Za-z]/g, "").toUpperCase();
    if (!seq) { skipped++; continue; }
    contigs.push({ id, length: seq.length });
    out.push(`>${id}\n${seq.replace(/(.{80})/g, "$1\n").replace(/\n$/, "")}\n`);
  }
  if (!contigs.length) throw new Error("this GenBank file has no sequence (no ORIGIN section); give the genome as FASTA");
  return { fasta: out.join(""), contigs, cds, skipped };
}
