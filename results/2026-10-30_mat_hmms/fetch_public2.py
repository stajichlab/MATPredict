"""Second pass: NCBI Protein by gene symbol (MATA1, MATALPHA2, MTLA1 ...) in Saccharomycotina; merges into data/public_*.faa."""
import json, time, urllib.parse, urllib.request
E_ = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(4):
        try: return urllib.request.urlopen(u, timeout=60).read().decode()
        except Exception: time.sleep(2 + 2 * i)
    raise SystemExit("fail")
G = {"a1": ["MATA1", "MATa1", "MTLA1", "mata1", "MAT-A1", "a1"], "alpha2": ["MATALPHA2", "MATalpha2", "MTLalpha2", "MAT-ALPHA2", "alpha2", "MATa2alpha"], "a2": ["MATA2", "MATa2", "MTLA2", "mata2", "MAT-A2"]}
old = json.load(open("data/public_ids.json"))
for c, gs in G.items():
    ids = set(old[c])
    for g in gs:
        term = f'{g}[Gene Name] AND txid147537[Organism:exp]'
        j = json.loads(get(E_ + "esearch.fcgi?db=protein&retmax=300&retmode=json&term=" + urllib.parse.quote(term)))
        ids |= set(j["esearchresult"]["idlist"]); time.sleep(0.4)
    ids = sorted(ids); print(c, len(old[c]), "->", len(ids)); old[c] = ids
    with open(f"data/public_{c}.faa", "w") as fh:
        for i in range(0, len(ids), 100):
            fh.write(get(E_ + "efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + ",".join(ids[i:i + 100]))); time.sleep(0.4)
json.dump(old, open("data/public_ids.json", "w"))
