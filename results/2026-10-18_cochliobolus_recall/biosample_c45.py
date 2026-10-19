"""Confirm strain / mating-type statements from NCBI assembly + BioSample for the four C4/C5 genomes (and report the attributes verbatim)."""
import json, re, urllib.request, xml.etree.ElementTree as ET
ASM = ['GCF_028858645.1', 'GCA_028891585.1', 'GCF_000354255.1', 'GCA_000338975.1']
import time
def get(u):
    time.sleep(1.5)
    return urllib.request.urlopen(u).read().decode()
E = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
out = []
for a in ASM:
    ids = json.loads(get(f'{E}esearch.fcgi?db=assembly&term={a}[Assembly+Accession]&retmode=json'))['esearchresult']['idlist']
    s = json.loads(get(f'{E}esummary.fcgi?db=assembly&id={ids[0]}&retmode=json'))['result'][ids[0]]
    bs = s.get('biosampleaccn')
    x = ET.fromstring(get(f'{E}efetch.fcgi?db=biosample&id={bs}&retmode=xml'))
    attrs = {at.get('harmonized_name') or at.get('attribute_name'): (at.text or '') for at in x.iter('Attribute')}
    title = ''.join(x.find('.//Title').itertext()) if x.find('.//Title') is not None else ''
    row = dict(assembly=a, name=s.get('assemblyname'), organism=s.get('speciesname'), strain=s.get('biosource', {}).get('infraspecieslist'), biosample=bs,
               bioproject=s.get('gb_bioprojects', [{}])[0].get('bioprojectaccn') if s.get('gb_bioprojects') else '', level=s.get('assemblystatus'),
               contig_n50=json.loads(s['meta'].replace("'", '"')) if False else '', title=title, attrs=attrs)
    out.append(row)
    print(json.dumps(row)[:900])
json.dump(out, open('results/2026-10-18_cochliobolus_recall/c45_biosample.json', 'w'), indent=1)
