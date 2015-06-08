#!/usr/bin/env python
import argparse
import tabix
import re

def main_vcf2hm450k(args):

    probes = []
    for line in open(args.p):
        fields = line.strip().split('\t')
        pname = fields[0]
        if pname.startswith('rs'):
            continue
        ptype = fields[1]
        chrm = 'chr'+fields[2]
        beg = int(fields[3])
        strand = fields[4]

        probes.append((pname, ptype, chrm, beg, strand))

    probes.sort()
    fh = tabix.open(args.f)
    for pname, ptype, chrm, beg, strand in probes:
        fnC,fnT,fnMM,fnG,fnA,fnNN = [0]*6
        rnC,rnT,rnMM,rnG,rnA,rnNN = [0]*6
        fBRC6 = 'NA'
        rBRC6 = 'NA'
        for ret in fh.query(chrm, beg, beg+1):
            if re.search(r'Context=C[GH]', ret[7]):
                if int(ret[1]) == beg:
                    fBRC6 = ret[9].split(':')[2]
                    fnC, fnT, fnMM, fnG, fnA, fnNN = map(int, fBRC6.split(','))
                if int(ret[1]) == beg+1 and ret[3] == 'G':
                    rBRC6 = ret[9].split(':')[2]
                    rnC, rnT, rnMM, rnG, rnA, rnNN = map(int, rBRC6.split(','))

        nC = fnC + rnC
        nT = fnT + rnT
        cov = nC + nT
        if cov >= 5:            # coverage threshold
            beta = '%1.3f' % (float(nC) / float(cov))
        else:
            beta = 'NA'
        print '%s\t%s\t%d\t%s\t%s' % (pname, beta, cov, fBRC6, rBRC6)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='convert WGBS vcf to hm450k')
    
    parser.add_argument('-f', help='vcf file, must be tabix indexed')
    parser.add_argument('-p', default='/data/largeS2/pl-bs/data/450k_probe_design')
    parser.set_defaults(func=main_vcf2hm450k)

    args = parser.parse_args()
    args.func(args)
