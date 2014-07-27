#!/RIS/HPC_apps/AMD/python/python-2.7.2/bin/python

import argparse
import os, sys
import gzip
from subprocess import Popen, PIPE, check_call
import re
from glob import glob

class VRecord:

    def __init__(self, chrm, pos):
        self.chrm = chrm
        self.pos = int(pos)

    def format_vcf(self, samples):
        return "%s\t%d\t%s\t%s\n" % (self.chrm, self.pos, "\t".join(self.text2to8), "\t".join([self.data[_] for _ in samples]))

class VCF:

    def __init__(self, fn):
        self.fn = fn
        self.header = ''
        self.chrmline = ''

    def open(self):
        if self.fn.endswith(".gz"):
            self.fh = gzip.open(self.fn, "r")
        else:
            self.fh = open(self.fn, "r")

    def close(self):
        self.fh.close()

    def read_header(self):

        while (True):
            line = self.fh.readline()
            if not line:
                break
            
            if line.startswith("#CHROM"):
                self.samplelist = line.strip().split()[9:]
                self.chrmline = '\t'.join(line.split()[:9])
                break
            else:
                self.header += line

    def format_header(self):
        return self.header+self.chrmline+'\t'+'\t'.join(self.samplelist)

    def format1(self, r):

        s = r.chrm
        s += '\t%d\t' % r.pos
        s += '\t'.join(r.text2to8)
        s += '\t'
        s += '\t'.join([r.data[_] for _ in self.samplelist])

        return s

    def read1(self):

        while (True):
            line = self.fh.readline()
            if not line:
                break
            if line[0] == '#':
                continue
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.id = pair[2]
            r.data = dict(zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]]))
            yield r

    def fetch_region(self, chrm, beg, end):

        for line in Popen(["tabix", self.fn, "%s:%d-%d" % (chrm, beg, end)], stdout=PIPE).stdout:
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.data = dict(zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]]))
            yield r

def main_compare(args):

    vcf1 = VCF(sys.argv[1])
    vcf2 = VCF(sys.argv[2])

    vcf1.open()
    vcf2.open()

    vcf1.read_header()
    vcf2.read_header()

    intersample = set(vcf2.samplelist) & set(vcf1.samplelist)
    print len(intersample), "overlapping samples"
    
    cs2g2 = {}
    for r in vcf2.read1():
        for s in intersample:
            cs2g2[(r.id, s)] = r.data[s]

    cs2g1 = {}
    for r in vcf1.read1():
        for s in intersample:
            cs2g1[(r.id, s)] = r.data[s]

    print len(cs2g1), "call-sample pairs in vcf1"
    print len(cs2g2), "call-sample pairs in vcf2"

    overlap = set(cs2g1.keys()) & set(cs2g2.keys())
    print len(overlap), "overlapping call-sample pairs between the two vcfs"

    vars = ["./.", "0/0", "0/1", "1/1"]
    print 'vcf1\\vcf2', '\t'.join(vars)
    for var1 in vars:
        sys.stdout.write(var1)
        for var2 in vars:
            sys.stdout.write('\t%d' % len([_ for _ in overlap if cs2g1[_] == var1 and cs2g2[_] == var2]))
        sys.stdout.write('\n')


def main_filter(args):

    vcf = VCF(args.v)
    vcf.open()
    vcf.read_header()

    cids = set()
    with open(args.cid) as fh:
        for line in fh:
            pair = line.split('\t')
            cids.add(pair[args.cidcol-1])

    print vcf.format_header()
    for r in vcf.read1():
        if r.id in cids:
            print vcf.format1(r)

    return
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='vcf tool')
    subparsers = parser.add_subparsers()

    psr_compare = subparsers.add_parser("compare", help=""" compare vcfs """)
    psr_compare.add_argument('-v1', help='VCF file 1')
    psr_compare.add_argument('-v2', help='VCF file 2')
    psr_compare.set_defaults(func=main_compare)

    psr_filter = subparsers.add_parser("filter", help=""" filter vcf """)
    psr_filter.add_argument('-v', help='VCF file')
    psr_filter.add_argument('--cid', help='call id list')
    psr_filter.add_argument('--cidcol', type=int, default=1, help='call id column index (1-based) [1]')
    psr_filter.set_defaults(func=main_filter)

    args = parser.parse_args()
    args.func(args)

