#!/usr/bin/env python3
"""Mean-pooled ESM-2 embeddings (last layer, special tokens excluded) per protein.
Usage: embed_esm.py IN.faa OUT_PREFIX MODEL [MODEL ...]
Writes OUT_PREFIX.<short>.npy (n x d, rows in FASTA order) and OUT_PREFIX.ids.
Proteins longer than 1022 aa are truncated (receptors are about 300-500 aa)."""
import sys
import time

import numpy as np
import torch
from transformers import AutoTokenizer, EsmModel

fa, prefix, models = sys.argv[1], sys.argv[2], sys.argv[3:]
ids, seqs = [], []
for line in open(fa):
    if line.startswith(">"):
        ids.append(line[1:].strip())
        seqs.append("")
    else:
        seqs[-1] += line.strip()
with open(prefix + ".ids", "w") as fo:
    fo.write("\n".join(ids) + "\n")
torch.set_num_threads(8)
for m in models:
    short = m.split("/")[-1].split("_")[1] + m.split("_")[2]
    t0 = time.time()
    tok = AutoTokenizer.from_pretrained(m)
    model = EsmModel.from_pretrained(m).eval()
    out = []
    with torch.no_grad():
        for s in seqs:
            enc = tok(s[:1022].replace("*", "").replace("X", "X"), return_tensors="pt")
            h = model(**enc).last_hidden_state[0, 1:-1]
            out.append(h.mean(0).numpy())
    np.save(f"{prefix}.{short}.npy", np.vstack(out))
    print(m, np.vstack(out).shape, f"{time.time() - t0:.0f}s", flush=True)
