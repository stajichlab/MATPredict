"""Public a1 / alpha2 / a2 proteins of Saccharomycotina (taxid 147537) from NCBI Protein by product name (E-utilities)."""
import json, re, sys, time, urllib.parse, urllib.request
def get(u):
    for i in range(4):
        try: return urllib.request.urlopen(u, timeout=60).read().decode()
        except Exception as e: time.sleep(2 + 2 * i)
    raise SystemExit("fail " + u)
E_ = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
Q = {
 "a1": ['"mating-type protein A1"', '"mating type protein a1"', '"mating-type a1"', '"MATa1"', '"MAT a1"', '"mating type protein A1"'],
 "alpha2": ['"mating-type protein alpha2"', '"mating type protein alpha2"', '"MATalpha2"', '"MAT alpha2"', '"mating-type alpha2"', '"mating type protein ALPHA2"'],
 "a2": ['"mating-type protein A2"', '"mating type protein a2"', '"MATa2"', '"MAT a2"', '"mating-type a2"', '"mating type protein A2"'],
}
res = {}
for c, qs in Q.items():
    ids = []
    for q in qs:
        term = f'({q}[Title]) AND txid147537[Organism:exp]'
        j = json.loads(get(E_ + "esearch.fcgi?db=protein&retmax=300&retmode=json&term=" + urllib.parse.quote(term)))
        ids += j["esearchresult"]["idlist"]; time.sleep(0.4)
    res[c] = sorted(set(ids)); print(c, len(res[c]))
json.dump(res, open("data/public_ids.json", "w"))
for c, ids in res.items():
    with open(f"data/public_{c}.faa", "w") as fh:
        for i in range(0, len(ids), 100):
            fh.write(get(E_ + "efetch.fcgi?db=protein&rettype=fasta&retmode=text&id=" + ",".join(ids[i:i + 100]))); time.sleep(0.4)
    print(c, open(f"data/public_{c}.faa").read().count(">"))
