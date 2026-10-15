// Run the miniprot WASM build on the same inputs as the native r281 binary and compare.
// Usage: node miniprot_wasm_test.mjs <build-dir> <input-dir> <extra args...>
// Each call gets a fresh module instance (the runtime exits after main()).
import { existsSync, readFileSync, statSync } from "node:fs";
import { join } from "node:path";
import { pathToFileURL } from "node:url";

const [buildDir, inputDir, ...args] = process.argv.slice(2);
const { default: createMiniprot } = await import(pathToFileURL(join(buildDir, "miniprot.mjs")).href);

async function runOnce(files, argv) {
  const stdout = [];
  const stderr = [];
  const mod = await createMiniprot({
    print: (s) => stdout.push(s),
    printErr: (s) => stderr.push(s),
    locateFile: (f) => join(buildDir, f),
  });
  for (const name of Object.keys(files)) mod.FS.writeFile(`/${name}`, files[name]);
  let code;
  try {
    mod.callMain(argv);
    code = 0;
  } catch (exc) {
    // Emscripten signals exit(n) by throwing ExitStatus.
    code = exc.status ?? `trap: ${exc.message}`;
    if (process.env.TRAP_STACK) process.stderr.write(`${exc.stack}\n`);
  }
  return { code, stdout: stdout.join("\n"), stderr: stderr.join("\n") };
}

// Mount every argument that names a file in inputDir (e.g. target.fa, allq.faa).
const files = {};
for (const a of args) {
  const p = join(inputDir, a);
  if (existsSync(p) && statSync(p).isFile()) files[a] = readFileSync(p);
}
const result = await runOnce(files, args);
process.stdout.write(result.stdout + "\n");
process.stderr.write(`exit=${result.code}\n${result.stderr}\n`);
