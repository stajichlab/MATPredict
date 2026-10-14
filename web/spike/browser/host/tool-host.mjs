// Tool host: runs WebAssembly tools for the Python pipeline.
//
// Each call creates a fresh module instance (Emscripten modules are not safe to re-run), copies the
// files the call touches in from the Python filesystem (the directories named by absolute argument
// paths), runs the tool, and copies new or changed files back out. A trap (abort, out-of-bounds)
// becomes a negative return code, which search.py already treats as a crash.
//
// `pyFS` is the Pyodide filesystem object (pyodide.FS). `factories` maps a tool name to a function
// returning a module-creating factory, so the same host works in Node and in a browser.

export class ToolHost {
  constructor({ pyFS, factories, wasiTools = {}, wasiRunner = null }) {
    this.pyFS = pyFS;
    this.factories = factories; // { mmseqs: createMmseqs, exonerate: createExonerate, ... }
    this.wasiTools = wasiTools; // { losat: compiled WebAssembly.Module } for WASI commands
    this.wasiRunner = wasiRunner; // { runWasiCommand, WORK_DIR } from wasi-runner.mjs
  }

  // argv[0] is the tool name, the rest are its arguments (as subprocess would receive them).
  async run(argv) {
    const [tool, ...args] = argv;
    if (this.wasiTools && this.wasiTools[tool]) return this.runWasi(tool, args);
    const factory = this.factories[tool];
    if (!factory) return { returncode: 127, stdout: "", stderr: `no WebAssembly tool named ${tool}` };

    const roots = this.roots(args);
    const out = [], err = [];
    const mod = await factory({
      print: (s) => out.push(s),
      printErr: (s) => err.push(s),
    });
    for (const root of roots) this.copyTree(this.pyFS, mod.FS, root);

    let returncode = 0;
    try {
      mod.callMain(args);
    } catch (exc) {
      if (exc && exc.status !== undefined) returncode = exc.status;
      else returncode = -6; // SIGABRT-like: a trap or abort, not a clean exit
    }
    for (const root of roots) this.copyTree(mod.FS, this.pyFS, root, { onlyNewer: true });
    return { returncode, stdout: out.join("\n") + (out.length ? "\n" : ""), stderr: err.join("\n") };
  }

  // WASI command (for example LOSAT). Absolute path arguments are copied into the WASI work directory
  // under their base names, and the arguments are rewritten to match.
  async runWasi(tool, args) {
    const { runWasiCommand, WORK_DIR } = this.wasiRunner;
    const inputs = new Map();
    const rewritten = args.map((a) => {
      if (typeof a !== "string" || !a.startsWith("/")) return a;
      const name = a.slice(a.lastIndexOf("/") + 1);
      try {
        inputs.set(name, this.pyFS.readFile(a, { encoding: "binary" }));
      } catch {
        return a; // not a file in the Python filesystem (for example an output path); leave it alone
      }
      return `${WORK_DIR}/${name}`;
    });
    const res = await runWasiCommand(this.wasiTools[tool], [tool, ...rewritten], inputs);
    return { returncode: res.returncode, stdout: res.stdout, stderr: res.stderr };
  }

  // Directories the call can read or write: the parent of every absolute path argument.
  roots(args) {
    const dirs = new Set();
    for (const a of args) {
      if (typeof a === "string" && a.startsWith("/")) {
        const i = a.lastIndexOf("/");
        dirs.add(i <= 0 ? "/" : a.slice(0, i));
      }
    }
    // Keep only the outermost directories, so trees are walked once.
    const sorted = [...dirs].sort((x, y) => x.length - y.length);
    const keep = [];
    for (const d of sorted) {
      if (!keep.some((k) => d === k || d.startsWith(k + "/"))) keep.push(d);
    }
    return keep;
  }

  // Copy every file under `root` from `src` FS to `dst` FS. With onlyNewer, skip files whose
  // contents already match, so unchanged inputs are not rewritten.
  copyTree(src, dst, root, { onlyNewer = false } = {}) {
    let names;
    try {
      names = src.readdir(root);
    } catch {
      return; // directory absent on this side
    }
    for (const name of names) {
      if (name === "." || name === "..") continue;
      const path = root === "/" ? `/${name}` : `${root}/${name}`;
      const st = src.stat(path);
      if (src.isDir(st.mode)) {
        this.ensureDir(dst, path);
        this.copyTree(src, dst, path, { onlyNewer });
      } else {
        const data = src.readFile(path, { encoding: "binary" });
        if (onlyNewer && this.sameFile(dst, path, data)) continue;
        this.ensureDir(dst, root);
        dst.writeFile(path, data);
      }
    }
  }

  ensureDir(fs, path) {
    const parts = path.split("/").filter(Boolean);
    let cur = "";
    for (const p of parts) {
      cur += "/" + p;
      try {
        fs.stat(cur);
      } catch {
        fs.mkdir(cur);
      }
    }
  }

  sameFile(fs, path, data) {
    try {
      const old = fs.readFile(path, { encoding: "binary" });
      if (old.length !== data.length) return false;
      for (let i = 0; i < data.length; i++) if (old[i] !== data[i]) return false;
      return true;
    } catch {
      return false;
    }
  }
}
