"""For genomes whose reports (minus run) differ between two variants: which top-level keys / paths differ.
usage: report_diff.py RUNS VA VB GENOME..."""
import sys, yaml, collections
def walk(a, b, path, out):
    if type(a) != type(b): out.append(path); return
    if isinstance(a, dict):
        for k in set(a) | set(b): walk(a.get(k), b.get(k), f'{path}.{k}', out)
    elif isinstance(a, list):
        if len(a) != len(b): out.append(path + '[len]'); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, f'{path}[]', out)
    elif a != b: out.append(path)
cnt = collections.Counter()
for g in sys.argv[4:]:
    A = yaml.load(open(f'{sys.argv[1]}/{sys.argv[2]}/{g}/detection_report.yaml'), Loader=yaml.CSafeLoader)
    B = yaml.load(open(f'{sys.argv[1]}/{sys.argv[3]}/{g}/detection_report.yaml'), Loader=yaml.CSafeLoader)
    A.pop('run'); B.pop('run'); out = []
    walk(A, B, '', out)
    cnt.update(set(out))
for k, v in cnt.most_common(): print(v, k)
