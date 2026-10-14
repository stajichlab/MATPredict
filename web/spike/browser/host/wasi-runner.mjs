// Run a WASI command (e.g. LOSAT) with an in-memory filesystem. The same code runs in Node and in the
// browser: @bjorn3/browser_wasi_shim keeps every file in memory, so nothing reaches the host filesystem.
import { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout, WASIProcExit } from "@bjorn3/browser_wasi_shim";

// Mount point inside the WASI filesystem for the input files.
export const WORK_DIR = "/w";

// Run `argv` (argv[0] is the program name) on `wasmModule` (a compiled WebAssembly.Module).
// `inputs` maps file names (no directories) to Uint8Array contents, all placed under WORK_DIR.
// Returns { returncode, stdout, stderr, files } where files maps names to Uint8Array for files the
// command wrote into WORK_DIR.
export async function runWasiCommand(wasmModule, argv, inputs = new Map()) {
  const stdout = [];
  const stderr = [];
  const fileMap = new Map();
  for (const [name, bytes] of inputs) fileMap.set(name, new File(bytes));
  const dir = new PreopenDirectory(WORK_DIR, fileMap);
  const fds = [
    new OpenFile(new File([])),
    ConsoleStdout.lineBuffered((line) => stdout.push(line)),
    ConsoleStdout.lineBuffered((line) => stderr.push(line)),
    dir,
  ];
  const wasi = new WASI(argv, [], fds, { debug: false });
  const instance = await WebAssembly.instantiate(wasmModule, { wasi_snapshot_preview1: wasi.wasiImport });
  let returncode = 0;
  try {
    wasi.start(instance);
  } catch (exc) {
    if (exc instanceof WASIProcExit) returncode = exc.code;
    else returncode = -6; // trap or abort, reported as a crash (negative), as the Python side expects
  }
  const files = new Map();
  for (const [name, file] of dir.dir.contents) files.set(name, file.data);
  return { returncode, stdout: stdout.join("\n") + (stdout.length ? "\n" : ""), stderr: stderr.join("\n"), files };
}
