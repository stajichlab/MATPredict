"""Fetch the BioSample (all attributes) and assembly strain/isolate for assembly accessions via NCBI E-utilities.
    fetch_biosample.py OUT.tsv GCA_...(.n) ...
"""
import sys, json, time, urllib.request, urllib.parse, re
import xml.etree.ElementTree as ET
P = "tool=matpredict&email=jason.stajich@ucr.edu"
def get(u):
    for _ in range(4):
        try:
            return urllib.request.urlopen(u, timeout=60).read().decode()
        except Exception:
            time.sleep(2)
    return ""
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
out = open(sys.argv[1], "w")
out.write("assembly\tbiosample\torganism\tstrain\tisolate\tsubmitter\tattributes\n")
for acc in sys.argv[2:]:
    j = json.loads(get(f"{E}esearch.fcgi?db=assembly&term={acc}%5BAssembly%20Accession%5D&retmode=json&{P}"))
    ids = j["esearchresult"]["idlist"]; time.sleep(0.4)
    if not ids:
        out.write(f"{acc}\tNA\n"); continue
    s = json.loads(get(f"{E}esummary.fcgi?db=assembly&id={ids[0]}&retmode=json&{P}"))["result"][ids[0]]
    bs = s.get("biosampleaccn", ""); time.sleep(0.4)
    attrs = ""; sub = s.get("submitterorganization", "")
    if bs:
        x = get(f"{E}efetch.fcgi?db=biosample&id={bs}&{P}"); time.sleep(0.4)
        try:
            r = ET.fromstring(x)
            attrs = "; ".join(f"{a.get('harmonized_name') or a.get('attribute_name')}={a.text}" for a in r.iter("Attribute"))
        except Exception:
            attrs = "PARSEFAIL"
    org = s.get("speciesname", "")
    st = s.get("meta", "")
    m = re.search(r"<Infraspecie.*?>(.*?)</", s.get("meta", ""))
    out.write("\t".join([acc, bs, org, str(s.get("biosource", {}).get("infraspecieslist", "")), "", sub, attrs]).replace("\n", " ") + "\n")
    out.flush()
