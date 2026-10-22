"""Fetch each truth gene's accession with 3 kb flanks. python fetch.py truth.tsv OUT_FA OUT_TSV"""
import csv, sys, time, urllib.request, urllib.error
rows = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
out, meta = open(sys.argv[2], "w"), []
for r in rows:
    ex = [tuple(map(int, e.split("-"))) for e in r["exons"].split(";")]
    lo, hi = min(a for a, b in ex), max(b for a, b in ex)
    s, e = max(1, lo - 3000), hi + 3000
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={r['accession']}&rettype=fasta&retmode=text&seq_start={s}&seq_stop={e}"
    seq = None
    for attempt in range(3):
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()
            seq = "".join(txt.splitlines()[1:]).upper(); break
        except Exception as ex_:
            time.sleep(2)
    key = f"{r['record']}|{r['gene']}"
    if not seq:
        print("FAIL", key, r["accession"]); continue
    out.write(f">{key}\n{seq}\n")
    meta.append(dict(key=key, seq_start=s, seq_len=len(seq), lo=lo, hi=hi))
    time.sleep(0.4)
with open(sys.argv[3], "w", newline="") as f:
    w = csv.DictWriter(f, meta[0].keys(), delimiter="\t"); w.writeheader(); w.writerows(meta)
print(len(meta), "of", len(rows), "fetched")
