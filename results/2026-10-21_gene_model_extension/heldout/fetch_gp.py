import urllib.request, urllib.parse, time, re, json
ids=open('ids.txt').read().split()
out=[]
for i in range(0,len(ids),150):
    b=ids[i:i+150]
    data=urllib.parse.urlencode(dict(db='protein',id=','.join(b),rettype='gp',retmode='text')).encode()
    for a in range(3):
        try:
            txt=urllib.request.urlopen(urllib.request.Request("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",data=data),timeout=120).read().decode(); break
        except Exception as e: txt=''; time.sleep(3)
    out.append(txt); time.sleep(0.5)
    if i%900==0: print(i,len(txt),flush=True)
open('gp_all.txt','w').write('\n'.join(out))
