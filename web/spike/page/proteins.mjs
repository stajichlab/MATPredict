// Predicted proteins of the MAT locus genes, translated from the gene models in the run's report:
// each gene's exons spliced in genomic order, reverse-complemented on the minus strand, and translated with
// the standard genetic code (NCBI table 1) up to the first stop. Pure functions, testable in Node.

// Standard code, codons in TCAG order (the order used by NCBI's translation tables).
const TCAG = "TCAG";
const CODE1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
const CODON = new Map();
for (let i = 0; i < 64; i++) {
  const c = TCAG[Math.floor(i / 16)] + TCAG[Math.floor(i / 4) % 4] + TCAG[i % 4];
  CODON.set(c, CODE1[i]);
}

export function translate(dna) {
  let protein = "";
  for (let i = 0; i + 3 <= dna.length; i += 3) {
    const aa = CODON.get(dna.slice(i, i + 3).toUpperCase()) ?? "X";
    if (aa === "*") break;
    protein += aa;
  }
  return protein;
}

const COMP = { A: "T", C: "G", G: "C", T: "A", N: "N" };
export function reverseComplement(dna) {
  let out = "";
  for (let i = dna.length - 1; i >= 0; i--) out += COMP[dna[i]] ?? "N";
  return out;
}

// The coding sequence of one gene model: exons in genomic order, reverse-complemented on the minus strand.
export function codingSequence(gene, seq) {
  const exons = (gene.exons && gene.exons.length) ? gene.exons : [{ start: gene.start, end: gene.end }];
  const ordered = [...exons].sort((a, b) => a.start - b.start);
  const dna = ordered.map((e) => seq.slice(e.start - 1, e.end)).join("");
  return gene.strand === "-" ? reverseComplement(dna) : dna;
}

// One FASTA record per gene of a detected locus (genes not present in the search are skipped).
export function locusProteinsFasta(detectedLocus, seq, label) {
  const records = [];
  for (const g of detectedLocus.gene_evidence || []) {
    // A set-aside gene of a resolved pair is not part of the locus (as in the GFF and the figure).
    if (g.lost_to || !g.contig || g.start == null || g.end == null) continue;
    const protein = translate(codingSequence(g, seq));
    if (!protein) continue;
    const span = `${g.contig}:${g.start}-${g.end}(${g.strand || "?"})`;
    const attrs = [g.role, g.identity != null ? `identity=${g.identity}%` : null, g.status].filter(Boolean).join(" ");
    records.push(`>${label}|${g.gene} ${span} ${attrs}`);
    for (let i = 0; i < protein.length; i += 60) records.push(protein.slice(i, i + 60));
  }
  return records.length ? records.join("\n") + "\n" : "";
}
