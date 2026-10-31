// node --test web/spike/page/genbank_input.test.mjs
import test from "node:test";
import assert from "node:assert/strict";
import { genbankToFasta, isGenbank } from "./genbank_input.mjs";

const REC = (id, seq, features = "") => `LOCUS       ${id.split(".")[0]}  ${seq.length} bp    DNA     linear   CON 09-FEB-2026
ACCESSION   ${id.split(".")[0]}
VERSION     ${id}
FEATURES             Location/Qualifiers
     source          1..${seq.length}
${features}ORIGIN
${seq.toLowerCase().replace(/(.{10})/g, "$1 ").replace(/(.{66})/g, "$1\n").split("\n").map((l, i) => String(i * 60 + 1).padStart(9) + " " + l).join("\n")}
//
`;

test("contig names are the VERSION accession and the sequence is upper case", () => {
  const r = genbankToFasta(REC("CM143129.1", "ACGTACGTACGTACGTACGT") + REC("CM143130.1", "GGGGCCCC"));
  assert.deepEqual(r.contigs, [{ id: "CM143129.1", length: 20 }, { id: "CM143130.1", length: 8 }]);
  assert.match(r.fasta, /^>CM143129\.1\nACGTACGTACGTACGTACGT\n>CM143130\.1\nGGGGCCCC\n$/);
});

test("CDS features are counted, a sequence-only file gives zero", () => {
  const withCds = genbankToFasta(REC("X1.1", "ACGT".repeat(10), "     CDS             1..30\n                     /locus_tag=\"A\"\n"));
  assert.equal(withCds.cds, 1);
  assert.equal(genbankToFasta(REC("X2.1", "ACGT".repeat(10))).cds, 0);
});

test("long sequences are wrapped at 80 columns", () => {
  const r = genbankToFasta(REC("L1.1", "A".repeat(200)));
  assert.deepEqual(r.fasta.trim().split("\n").map((l) => l.length), [5, 80, 80, 40]);
});

test("a record without ORIGIN is skipped; none with sequence throws", () => {
  const noSeq = "LOCUS       Z1  0 bp    DNA\nVERSION     Z1.1\nFEATURES\n//\n";
  assert.equal(genbankToFasta(noSeq + REC("Y1.1", "ACGT")).skipped, 1);
  assert.throws(() => genbankToFasta(noSeq), /no sequence/);
});

test("GenBank is recognised by suffix, gzip suffix or a LOCUS first line", () => {
  assert.ok(isGenbank("a.gbff", ""));
  assert.ok(isGenbank("a.GBK.gz", ""));
  assert.ok(isGenbank("a.txt", "\n\nLOCUS       X  10 bp"));
  assert.ok(!isGenbank("a.fna", ">c1\nACGT\n"));
});
