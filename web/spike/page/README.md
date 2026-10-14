# Browser page (spike)

A standalone page that runs `matpredict detect` on an uploaded genome. Python runs in Pyodide 0.27.7
(pinned on the jsDelivr CDN); the search tools run as WebAssembly: miniprot and exonerate (Emscripten,
in-memory filesystem), MMseqs2 (Emscripten, unthreaded), LOSAT (WASI, `browser_wasi_shim`). The
Python bridge (`../browser/python/bridge.py`) routes subprocess calls to them. Nothing is uploaded.

## Build and run

    bash ../browser/build-bundle.sh <work-dir> <bundle-dir>   # WebAssembly tools + pyhmmer wheel
    bash assemble.sh <bundle-dir>                              # links, MATPredict wheel, db manifest
    python3 -m http.server 8000 --directory ..                 # serve web/spike (the page imports ../browser, ../localizer)
    open http://localhost:8000/page/index.html

`_assets/` is generated and not committed. Serve from `web/spike/` (not `page/`) so the page can import
the shared host and bridge modules.

## Output

- Results table: detected loci (family, contig, span, idiomorph, confidence, genes found and missing).
- Locus figure: genes on the contig with flanking bases (flank selectable), SVG.
- Flanking sequence: 60 bases per line with 1-based coordinates; locus bases highlighted.
- The pipeline's own `report.html` (iframe), and downloads for `detected_loci.gff3` and `report.html`.
- Timing: total and detect time. Memory: Pyodide's heap, plus the JS heap where the browser reports it.
  The WebAssembly tool modules' memory is not measured.

## Limits

- The whole genome is read into memory, and detect runs on the main thread (no worker yet), so the page
  is unresponsive during a run.
- Threaded MMseqs2 is not built; the tblastn stand-in runs single-threaded.
- The MATPredict wheel is copied from `web/spike/dist`, which must be rebuilt after source changes.
