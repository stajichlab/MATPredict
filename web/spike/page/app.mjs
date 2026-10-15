// MATPredict detection in the browser. Pyodide runs MATPredict's Python; the search tools run as
// WebAssembly modules through the same tool host and Python bridge as the Node harness (run-detect.mjs).
import { loadPyodide } from "https://cdn.jsdelivr.net/pyodide/v0.27.7/full/pyodide.mjs";
import { ToolHost } from "../browser/host/tool-host.mjs";
import { runWasiCommand, WORK_DIR } from "../browser/host/wasi-runner.mjs";
import createMiniprot from "./_assets/bundle/miniprot/miniprot.mjs";
import { locusFasta, locusGenbank } from "./genbank.mjs";

const ASSETS = "./_assets";
const $ = (id) => document.getElementById(id);
const logEl = $("log");
const t0 = performance.now();
const elapsed = () => ((performance.now() - t0) / 1000).toFixed(1);
function log(msg) {
  logEl.textContent += `\n[${elapsed()}s] ${msg}`;
  logEl.scrollTop = logEl.scrollHeight;
}
const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

async function fetchBytes(url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`fetch ${url}: HTTP ${r.status}`);
  return new Uint8Array(await r.arrayBuffer());
}

// Pyodide's own linear memory, plus the JS heap where the browser reports it (Chrome, with
// --enable-precise-memory-info; otherwise coarse).
function memoryReport(py) {
  const parts = [];
  try { parts.push(`Pyodide heap ${(py._module.HEAPU8.length / 2 ** 20).toFixed(0)} MB`); } catch {}
  if (performance.memory) parts.push(`JS heap ${(performance.memory.usedJSHeapSize / 2 ** 20).toFixed(0)} MB`);
  return parts.join(", ");
}

// Mount the database from the manifest. Fetches run in small batches.
async function mountDb(py, manifest) {
  const dirs = new Set();
  for (const rel of manifest.db) {
    const parts = rel.split("/");
    for (let i = 1; i < parts.length; i++) dirs.add("/db/" + parts.slice(0, i).join("/"));
  }
  py.FS.mkdirTree("/db");
  for (const d of [...dirs].sort()) py.FS.mkdirTree(d);
  for (let i = 0; i < manifest.db.length; i += 16) {
    const batch = manifest.db.slice(i, i + 16);
    const datas = await Promise.all(batch.map((rel) => fetchBytes(`${ASSETS}/db/${rel}`)));
    batch.forEach((rel, j) => py.FS.writeFile(`/db/${rel}`, datas[j]));
  }
  return manifest.db.length;
}

// Python side: install the wheel, route subprocess to the host, run `matpredict detect`, and return a
// JSON summary of the outputs. Same flow as run-detect.mjs.
function pythonRun({ genomeName, phylum, engine, pyhmmerWheel }) {
  return `
import os, sys, zipfile, argparse, traceback, json
os.environ["MATPREDICT_DB_ROOT"] = "/db"
zipfile.ZipFile("/tmp/matpredict.whl").extractall("/site_pkgs")
sys.path.insert(0, "/psutil_stub")
import micropip
await micropip.install("emfs:/${pyhmmerWheel}", deps=False)
sys.path[:0] = ["/site_pkgs", "/browser"]
import bridge
bridge.ENGINE = "${engine}"
from MATPredict.detect import cli as detect_cli
parser = argparse.ArgumentParser(prog="matpredict")
subs = parser.add_subparsers(dest="cmd")
detect_cli.register_subcommands(subs)
args = parser.parse_args(["detect", "--genome", "/data/${genomeName}", "--phylum", "${phylum}", "--out-dir", "/out"])
error = None
rc = None
try:
    rc = args.func(args)
except BaseException as exc:
    error = f"{type(exc).__name__}: {exc}\\n" + traceback.format_exc(limit=6)
import yaml
def read(path):
    return open(path).read() if os.path.exists(path) else ""
report = yaml.safe_load(read("/out/detection_report.yaml")) or {}
json.dumps({
    "rc": rc, "error": error,
    "detected": report.get("detected") or [],
    "notDetected": report.get("not_detected") or [],
    "suppressed": report.get("suppressed_unpolished"),
    "gff": read("/out/detected_loci.gff3"),
    "html": read("/out/report.html"),
})
`;
}

// FASTA text -> Map(contig name -> sequence).
function parseFasta(text) {
  const seqs = new Map();
  let name = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      name = line.slice(1).trim().split(/\s+/)[0];
      seqs.set(name, []);
    } else if (name && line.trim()) {
      seqs.get(name).push(line.trim());
    }
  }
  const out = new Map();
  for (const [k, parts] of seqs) out.set(k, parts.join("").toUpperCase());
  return out;
}

// GFF from the pipeline -> [{ id, contig, start, end, genes: [{name, role, start, end, strand, identity}] }].
function parseLoci(gff) {
  const loci = new Map();
  const genes = [];
  for (const line of gff.split("\n")) {
    if (!line || line.startsWith("#")) continue;
    const f = line.split("\t");
    if (f.length < 9) continue;
    const attrs = Object.fromEntries(f[8].split(";").map((kv) => kv.split("=")));
    if (f[2] === "MAT_locus") {
      loci.set(attrs.ID, { id: attrs.ID, contig: f[0], start: +f[3], end: +f[4], family: attrs.family,
        idiomorph: attrs.idiomorph, confidence: attrs.confidence, genes: [] });
    } else if (f[2] === "gene" && attrs.present === "true") {
      genes.push({ parent: attrs.Parent, name: attrs.Name, role: attrs.role, start: +f[3], end: +f[4], strand: f[6],
        identity: attrs.identity, method: attrs.method, status: attrs.status });
    }
  }
  for (const g of genes) if (loci.has(g.parent)) loci.get(g.parent).genes.push(g);
  return [...loci.values()];
}

// SVG track: locus shaded, genes as arrows, coordinates along the top.
function locusSvg(locus, flank, contigLength) {
  const w0 = Math.max(1, locus.start - flank);
  const w1 = Math.min(contigLength, locus.end + flank);
  const W = 900, pad = 16, H = 120;
  const x = (p) => pad + ((p - w0) / Math.max(1, w1 - w0)) * (W - 2 * pad);
  const parts = [];
  parts.push(`<rect x="${x(locus.start)}" y="30" width="${Math.max(2, x(locus.end) - x(locus.start))}" height="60" style="fill:var(--locus)"/>`);
  parts.push(`<line x1="${pad}" x2="${W - pad}" y1="60" y2="60" style="stroke:var(--muted)" stroke-width="1"/>`);
  for (const g of locus.genes) {
    const gx = x(g.start), gw = Math.max(3, x(g.end) - x(g.start));
    const dir = g.strand === "-" ? -1 : 1;
    const cy = 60, h = 14;
    parts.push(`<rect x="${gx}" y="${cy - h / 2}" width="${gw}" height="${h}" style="fill:var(--gene)" opacity="0.85"/>`);
    parts.push(`<text x="${gx + gw / 2}" y="${cy - 12}" text-anchor="middle" font-size="11" style="fill:var(--fg)">${esc(g.name)}${dir < 0 ? " ←" : " →"}</text>`);
  }
  parts.push(`<text x="${pad}" y="${H - 6}" font-size="11" style="fill:var(--muted)">${w0.toLocaleString()}</text>`);
  parts.push(`<text x="${W - pad}" y="${H - 6}" font-size="11" style="fill:var(--muted)" text-anchor="end">${w1.toLocaleString()}</text>`);
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="Locus ${esc(locus.id)} on ${esc(locus.contig)}">${parts.join("")}</svg>`;
}

// Flanking sequence, 60 per line with 1-based coordinates; the locus bases are highlighted.
// Shown inside a collapsed block; the downloads are the primary way to use it.
function sequenceHtml(locus, seq, flank) {
  const w0 = Math.max(1, locus.start - flank);
  const w1 = Math.min(seq.length, locus.end + flank);
  const lines = [];
  for (let p = w0; p <= w1; p += 60) {
    const q = Math.min(w1, p + 59);
    let row = `${String(p).padStart(9)}  `;
    let inLocus = null;
    for (let i = p; i <= q; i++) {
      const isLocus = i >= locus.start && i <= locus.end;
      if (isLocus !== inLocus) {
        if (inLocus) row += "</span>";
        if (isLocus) row += '<span class="locus">';
        inLocus = isLocus;
      }
      row += esc(seq[i - 1]);
    }
    if (inLocus) row += "</span>";
    lines.push(row);
  }
  return lines.join("\n");
}

function renderResults(res, genomeText) {
  $("results").hidden = false;
  const tbody = $("loci").querySelector("tbody");
  tbody.innerHTML = "";
  for (const d of res.detected) {
    const span = d.supported_span || d;
    const tr = document.createElement("tr");
    tr.innerHTML = [
      d.family, d.contig, `${span.start?.toLocaleString?.() ?? span.start}–${span.end?.toLocaleString?.() ?? span.end}`,
      d.idiomorph, d.confidence, (d.genes_found || []).join(", ") || "—", (d.genes_missing || []).join(", ") || "—",
    ].map((v) => `<td>${esc(v)}</td>`).join("");
    tbody.appendChild(tr);
  }
  $("none").textContent = res.detected.length ? "" : "No MAT locus detected.";

  const contigs = parseFasta(genomeText);
  const loci = parseLoci(res.gff);
  const flank = +$("flank").value;
  $("figure").innerHTML = loci.map((l) => locusSvg(l, flank, (contigs.get(l.contig) || "").length)).join("");
  // Downloads per locus: the locus with its flanks as FASTA, and as an annotated GenBank record.
  const blobLink = (text, type, name) => {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type }));
    a.download = name;
    a.textContent = name;
    return a;
  };
  const box = $("locus-downloads");
  box.innerHTML = "";
  for (const l of loci) {
    const seq = contigs.get(l.contig) || "";
    const row = document.createElement("p");
    row.className = "dl-row";
    row.append(`${l.id} (${l.contig}, ${flank.toLocaleString()} bp flank): `);
    row.append(blobLink(locusFasta(l, seq, flank), "text/plain", `${l.id}.flank${flank}.fna`), " · ");
    row.append(blobLink(locusGenbank(l, seq, flank), "text/plain", `${l.id}.flank${flank}.gbk`));
    box.appendChild(row);
  }
  $("sequence").innerHTML = loci.map((l) => `<div class="muted">${esc(l.id)} · ${esc(l.contig)}:${l.start}–${l.end} (locus in highlight)</div>${sequenceHtml(l, contigs.get(l.contig) || "", flank)}`).join("\n\n");
  $("report").srcdoc = res.html;
  const dl = (id, text, type, name) => {
    const url = URL.createObjectURL(new Blob([text], { type }));
    $(id).href = url; $(id).download = name;
  };
  dl("dl-gff", res.gff, "text/plain", "detected_loci.gff3");
  dl("dl-html", res.html, "text/html", "report.html");
}

// A genome as the pipeline needs it: plain FASTA bytes and text. A .gz upload is decompressed in the
// browser (DecompressionStream), so nothing is sent anywhere and no server is needed.
let genome = null;
async function readGenome(file) {
  const gzipped = /\.gz$/i.test(file.name);
  const bytes = gzipped
    ? new Uint8Array(await new Response(file.stream().pipeThrough(new DecompressionStream("gzip"))).arrayBuffer())
    : new Uint8Array(await file.arrayBuffer());
  return { name: file.name.replace(/\.gz$/i, ""), bytes, text: new TextDecoder().decode(bytes), gzipped };
}
$("genome").addEventListener("change", async (e) => {
  const file = e.target.files[0] || null;
  genome = file ? await readGenome(file) : null;
  $("run").disabled = !genome;
  if (genome) {
    log(`genome ${file.name} (${(file.size / 1e6).toFixed(2)} MB${genome.gzipped ? `, decompressed to ${(genome.bytes.length / 1e6).toFixed(2)} MB` : ""})`);
  }
});
$("flank").addEventListener("change", () => {
  if (window.lastResult) renderResults(window.lastResult, genome.text);
});

$("run").addEventListener("click", async () => {
  $("run").disabled = true;
  const started = performance.now();
  try {
    log("loading Pyodide 0.27.7");
    const py = await loadPyodide();
    log("loading Python packages");
    await py.loadPackage(["biopython", "pyyaml", "jsonschema", "zstandard", "requests", "micropip"]);
    const manifest = await (await fetch(`${ASSETS}/db-manifest.json`)).json();
    log(`mounting database (${manifest.db.length} files)`);
    await mountDb(py, manifest);

    py.FS.mkdirTree("/data");
    py.FS.writeFile(`/data/${genome.name}`, genome.bytes);
    py.FS.mkdirTree("/out");
    py.FS.writeFile("/tmp/matpredict.whl", await fetchBytes(`${ASSETS}/matpredict.whl`));
    py.FS.writeFile(`/${manifest.pyhmmerWheel}`, await fetchBytes(`${ASSETS}/bundle/${manifest.pyhmmerWheel}`));
    py.FS.mkdirTree("/psutil_stub");
    py.FS.writeFile("/psutil_stub/psutil.py", await fetchBytes("../browser/pyhmmer/psutil.py"));
    py.FS.mkdirTree("/browser");
    py.FS.writeFile("/browser/localizer.py", await fetchBytes("../localizer/localizer.py"));
    py.FS.writeFile("/browser/bridge.py", await fetchBytes("../browser/python/bridge.py"));

    log("compiling WebAssembly tools");
    const factories = {
      mmseqs: (opts) => createMmseqs({ ...opts, locateFile: (f) => `${ASSETS}/bundle/mmseqs/${f}` }),
      exonerate: (opts) => createExonerate({ ...opts, locateFile: (f) => `${ASSETS}/bundle/exonerate/${f}` }),
      miniprot: (opts) => createMiniprot({ ...opts, locateFile: (f) => `${ASSETS}/bundle/miniprot/${f}` }),
    };
    const wasiTools = { losat: await WebAssembly.compile(await fetchBytes(`${ASSETS}/bundle/losat/LOSAT.wasm`)) };
    const host = new ToolHost({ pyFS: py.FS, factories, wasiTools, wasiRunner: { runWasiCommand, WORK_DIR } });
    globalThis.js_run = (argvJson) => host.run(JSON.parse(argvJson)).then((r) => JSON.stringify(r));

    const engine = $("engine").value;
    const phylum = $("phylum").value;
    log(`running detect (${phylum}, localizer ${engine})`);
    const detectStart = performance.now();
    const summary = JSON.parse(await py.runPythonAsync(pythonRun({
      genomeName: genome.name, phylum, engine, pyhmmerWheel: manifest.pyhmmerWheel,
    })));
    const detectSecs = ((performance.now() - detectStart) / 1000).toFixed(1);
    if (summary.error) log(`detect failed: ${summary.error}`);
    log(`detect return code ${summary.rc}; ${summary.detected.length} locus/loci; ${detectSecs}s`);
    $("timing").textContent = `Total ${((performance.now() - started) / 1000).toFixed(1)} s (detect ${detectSecs} s). ${memoryReport(py)}.`;
    window.lastResult = summary;
    renderResults(summary, genome.text);
  } catch (err) {
    log(`error: ${err && err.stack ? err.stack : err}`);
  } finally {
    $("run").disabled = false;
  }
});
