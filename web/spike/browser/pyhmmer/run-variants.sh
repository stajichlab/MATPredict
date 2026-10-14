#!/bin/bash
# Build pyhmmer for Pyodide with one patched file reverted at a time, to find the MSV difference.
set -u
source "$HOME/emsdk-3.1.58/emsdk_env.sh" >/dev/null 2>&1
BASE="$HOME/wasm-build/pyhmmer"
for v in A B C; do
  D="$BASE/var-$v/pyhmmer-0.12.3"
  mkdir -p "$BASE/var-$v"
  cp -R "$BASE/pyhmmer-0.12.3" "$D"
  case $v in
    A) cp "$BASE/pristine/pyhmmer-0.12.3/vendor/easel/esl_cpu.c" "$D/vendor/easel/esl_cpu.c" ;;
    B) cp "$BASE/pristine/pyhmmer-0.12.3/vendor/easel/esl_sse.c" "$D/vendor/easel/esl_sse.c"
       cp "$BASE/pristine/pyhmmer-0.12.3/vendor/easel/esl_sse.h" "$D/vendor/easel/esl_sse.h" ;;
    C) cp "$BASE/pristine/pyhmmer-0.12.3/vendor/easel/esl_stopwatch.c" "$D/vendor/easel/esl_stopwatch.c" ;;
  esac
  CFLAGS="-O2 -I$BASE/shim" LDFLAGS="-sWASM_BIGINT" PATH="$BASE/pb312/bin:$PATH" \
    "$BASE/pb312/bin/pyodide" build "$D" --outdir "$BASE/var-$v/dist" > "$BASE/var-$v/build.log" 2>&1
  echo "variant $v build=$?" >> "$BASE/variants.done"
done
echo "all builds finished" >> "$BASE/variants.done"
