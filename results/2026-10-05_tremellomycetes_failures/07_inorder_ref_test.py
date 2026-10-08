#!/usr/bin/env python3
"""Direct test of an in-order reference: HD1-class and STE3 protein fragments taken from called
Trichosporonales genomes (stop-to-stop six-frame segments carrying a Pfam KN or STE3 hit, so
exon-level fragments, not gene models) are aligned with miniprot to every genome of a target list.

usage: 07_inorder_ref_test.py MAKE REF_GENOMES_TXT OUT_FAA     (pixi python on HPCC; scan_dir has genes.tsv)
       07_inorder_ref_test.py RUN ASMID REF_FAA OUTDIR
RUN writes OUTDIR/<ASMID>.inorder.tsv: contig, strand, start, end, kind, query, identity, qcov.
"""
import os, sys, subprocess
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scan_genome as sg
import scan_ste3_hd as sh


def make(refs, out, scan_dir):
    fo = open(out, 'w')
    for a in open(refs).read().split():
        seqs = sg.read_fasta(os.path.join(sg.LIB, a + '.fa.gz'))
        want = []
        for line in open(os.path.join(scan_dir, a + '.genes.tsv')):
            f = line.rstrip('\n').split('\t')
            if f[5] in ('hmm:PF05920', 'hmm:PF02076'):
                want.append((f[0], f[1], int(f[3]), int(f[4]), f[5]))
        for kind, ctg, s, e, src in want:
            # stop-to-stop segments overlapping the hit on the right strand
            for c, strand, fr, aa0, seg, L in sh.six_frame_segments({ctg: seqs[ctg]}):
                nt0 = fr + 3 * aa0
                nt1 = nt0 + 3 * len(seg)
                a_, b_ = (nt0 + 1, nt1) if strand == '+' else (L - nt1 + 1, L - nt0)
                if a_ <= e and s <= b_ and 40 <= len(seg) <= 1200:
                    k = 'KN' if src.endswith('PF05920') else 'STE3'
                    fo.write('>%s|%s|%s:%d-%d\n%s\n' % (a, k, ctg, a_, b_, seg))
    fo.close()


def run(asm, ref, outdir):
    gz = os.path.join(sg.LIB, asm + '.fa.gz')
    paf = os.path.join(outdir, asm + '.mp.paf')
    with open(paf, 'w') as fo:
        sg.run([sg.BIN + '/miniprot', '-t', '2', '-I', '--outn=30', gz, ref], stdout=fo, stderr=subprocess.DEVNULL)
    with open(os.path.join(outdir, asm + '.inorder.tsv'), 'w') as fo:
        for line in open(paf):
            f = line.rstrip('\n').split('\t')
            if len(f) < 12:
                continue
            qlen, qs, qe = int(f[1]), int(f[2]), int(f[3])
            ident = int(f[9]) / max(1, int(f[10]))
            kind = f[0].split('|')[1]
            fo.write('\t'.join(map(str, [f[5], f[4], int(f[7]) + 1, int(f[8]), kind, f[0], '%.3f' % ident, '%.2f' % ((qe - qs) / qlen)])) + '\n')
    os.remove(paf)


if __name__ == '__main__':
    if sys.argv[1] == 'MAKE':
        make(sys.argv[2], sys.argv[3], sys.argv[4])
    else:
        run(sys.argv[2], sys.argv[3], sys.argv[4])
