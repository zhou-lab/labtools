#!/usr/bin/env python

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
                self.n = len(self.samplelist)
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
        s += '\t'.join([':'.join([r.data[_][_f] for _f in r.dataformat]) for _ in self.samplelist])

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
            r.ref = pair[3]
            r.alt = pair[4]
            r.info = dict([_.split('=',1) for _ in pair[7].split(';')])
            r.text2to8 = pair[2:9]
            r.id = pair[2]
            r.dataformat = pair[8].split(":")
            r.data = dict(zip(self.samplelist, [dict(zip(r.dataformat, _.strip('\n').split(":"))) for _ in pair[9:]]))
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

def bed_format_cg(vcf, p):
    print '%s\t%d\t%d\t%s' % (
        p.chrm, p.pos-1, p.end, '\t'.join(['%1.2f' % (float(p.retention[_]) /p.cov[_]) for _ in xrange(vcf.n)]))

def main_vcf2cg(args):

    vcf = VCF(args.i)
    vcf.open()
    vcf.read_header()
    p=None
    n = len(vcf.samplelist)
    for r in vcf.read1():
        merged=False
        if p:
            if p.pos+1==r.pos and p.ref=='C' and r.ref=='G': # merge CG
                if 'CX' in r.info and r.info['CX'] == 'CG':
                    p.end += 1
                    for i, sample in enumerate(vcf.samplelist):
                        p.retention[i] += float(r.data[sample]['BT']) * int(r.data[sample]['CV'])
                        p.cov[i] += int(r.data[sample]['CV'])
                        merged = True

            if any(_ > args.k for _ in p.cov):
                bed_format_cg(vcf, p)

        if not merged and 'CX' in r.info and r.info['CX'] == 'CG':
            r.end = r.pos
            r.retention = [0]*vcf.n
            r.cov = [0]*vcf.n
            for i, sample in enumerate(vcf.samplelist):
                r.retention[i] = float(r.data[sample]['BT']) * int(r.data[sample]['CV'])
                r.cov[i] = int(r.data[sample]['CV'])
            p = r
        else:
            p = None

    if p and any(_ > args.k for _ in p.cov):
        bed_format_cg(vcf, p)

    return

def main_vcf2hcg(args):

    vcf = VCF(args.i)
    vcf.open()
    vcf.read_header()
    p=None
    n = len(vcf.samplelist)
    for r in vcf.read1():
        merged=False
        if p:
            if p.pos+1==r.pos and p.ref=='C' and r.ref=='G': # merge CG
                if 'N5' in r.info and r.info['N5'][2:4] == 'CG' and r.info['N5'][1] != 'G':
                    p.end += 1
                    for i, sample in enumerate(vcf.samplelist):
                        p.retention[i] += float(r.data[sample]['BT']) * int(r.data[sample]['CV'])
                        p.cov[i] += int(r.data[sample]['CV'])
                        merged = True

            if any(_ > args.k for _ in p.cov):
                bed_format_cg(vcf, p)

        if not merged and 'N5' in r.info and r.info['N5'][2:4] == 'CG' and r.info['N5'][1] != 'G':
            r.end = r.pos
            r.retention = [0]*vcf.n
            r.cov = [0]*vcf.n
            for i, sample in enumerate(vcf.samplelist):
                r.retention[i] = float(r.data[sample]['BT']) * int(r.data[sample]['CV'])
                r.cov[i] = int(r.data[sample]['CV'])
            p = r
        else:
            p = None

    if any(_ > args.k for _ in p.cov):
        bed_format_cg(vcf, p)

    return

def main_vcf2gch(args):

    vcf = VCF(args.i)
    vcf.open()
    vcf.read_header()
    n = len(vcf.samplelist)
    for r in vcf.read1():
        if 'N5' in r.info and r.info['N5'][1:3] == 'GC' and r.info['N5'][3] != 'G':
            r.end = r.pos
            r.retention = [0]*vcf.n
            r.cov = [0]*vcf.n
            for i, sample in enumerate(vcf.samplelist):
                r.retention[i] = float(r.data[sample]['BT']) * int(r.data[sample]['CV'])
                r.cov[i] = int(r.data[sample]['CV'])

            if any(_ > args.k for _ in r.cov):
                bed_format_cg(vcf, r)

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='vcf tool')
    subparsers = parser.add_subparsers()

    parser_vcf2cg = subparsers.add_parser('vcf2cg', help='biscuit vcf to CpG file')
    parser_vcf2cg.add_argument('-i', required=True, help='input vcf')
    parser_vcf2cg.add_argument('-o', help='output prefix', default=sys.stdout)
    parser_vcf2cg.add_argument('-k', type=int, default=3, help='minimum coverage [3]')
    parser_vcf2cg.set_defaults(func=main_vcf2cg)

    parser_vcf2hcg = subparsers.add_parser('vcf2hcg', help='biscuit vcf to HCG bed file')
    parser_vcf2hcg.add_argument('-i', required=True, help='input vcf')
    parser_vcf2hcg.add_argument('-o', help='output prefix', default=sys.stdout)
    parser_vcf2hcg.add_argument('-k', type=int, default=3, help='minimum coverage [3]')
    parser_vcf2hcg.set_defaults(func=main_vcf2hcg)

    parser_vcf2gch = subparsers.add_parser('vcf2gch', help='biscuit vcf to GCH bed file')
    parser_vcf2gch.add_argument('-i', required=True, help='input vcf')
    parser_vcf2gch.add_argument('-o', help='output prefix', default=sys.stdout)
    parser_vcf2gch.add_argument('-k', type=int, default=3, help='minimum coverage [3]')
    parser_vcf2gch.set_defaults(func=main_vcf2gch)
    
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

