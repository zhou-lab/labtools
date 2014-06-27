#!/usr/bin/env python

import pysam
import argparse
import sys, __builtin__

class Read():

    def __init__(self, y = None, x = None, xob = None):
        """ x is pysam.AlignedRead,
        y is read file first line """
        self.keys = []
        if y:
            self.qname = y[0]
            self.read1or2 = int(y[1])
        elif x:
            self.qname = x.qname
            self.read1or2 = 1 if x.is_read1 else 2

            self.chrm = xob.getrname(x.tid)
            self.keys.append('chrm')

            self.pos = x.pos+1
            self.keys.append('pos')

            if x.cigar:
                self.cigar = x.cigar
                self.keys.append('cigar')
                self.cgstr = x.cigarstring
                self.keys.append('cgstr')

            self.mapq = x.mapq
            self.keys.append('mapq')

            self.seq = x.seq
            self.keys.append('seq')

            self.qual = x.qual
            self.keys.append('qual')

            self.unmapped = 1 if x.is_unmapped else 0
            self.keys.append('unmapped')

            self.duplicate = 1 if x.is_duplicate else 0
            self.keys.append('duplicate')

            if hasattr(x, "tags"):
                for opt, val in x.tags:
                    if opt == "NM":
                        self.edit = x.opt('NM')
                        self.keys.append('edit')
                

    def format(self):
        s = ">READ\t{0}\t{1}".format(self.qname, self.read1or2)
        for key in self.keys:
            val = getattr(self, key)
            vt = type(val).__name__
            if vt == 'float64':
                vt = 'float'
            s += "\n{0}\t{1}\t{2}".format(key, vt, val)

        s+='\n'

        return s

    def set(self, key, kt, val):

        if key not in self.keys:
            self.keys.append(key)

        if kt == "list":
            val = eval(val)
        else:
            if kt == 'float64':
                kt = 'float'
            val = getattr(__builtin__, kt)(val)

        setattr(self, key, val)

        return

def parse_readfile(fh):

    read = None
    for line in fh:
        pair = line.strip().split("\t")
        if not pair:
            continue
        if pair[0] == '>READ':
            if read:
                yield read
            read = Read(pair[1:])
        elif read and pair[0]:
            try:
                read.set(pair[0], pair[1], pair[2])
            except:
                raise Exception("%s\n%s\n" % (str(pair), line))
                
    if read:
        yield read

def main_tabulate(args):
    
    cols = args.p.split(',')
    for call in parse_readfile(args.read_file):
        print '\t'.join([str(getattr(call, k)) if hasattr(call, k) else 'NA' for k in cols])

def main_filter(args):

    for read in parse_readfile(args.read_file):
        if eval(args.s):
            print read.format()

def main_bam_region(args):

    bam = pysam.Samfile(args.bam_fn)
    for x in bam.fetch(region=args.reg):
        print Read(x=x, xob=bam).format()

def main_bam_reads(args):
    
    qnames = set()
    for read in parse_readfile(args.reads):
        qnames.add((read.qname, read.read1or2))

    bam = pysam.Samfile(args.bam_fn)
    for x in bam.fetch():
        read = Read(x=x, xob=bam)
        if (read.qname, read.read1or2) in qnames:
            print read.format()

def main_compare(args):

    r2aln1 = {}
    for r in parse_readfile(args.r1):
        r2aln1[(r.qname, r.read1or2)] = r

    r2aln2 = {}
    for r in parse_readfile(args.r2):
        r2aln2[(r.qname, r.read1or2)] = r

    reads1 = set(r2aln1.keys())
    reads2 = set(r2aln2.keys())

    print "%d reads in common" % len(reads1 & reads2)
    print "%d reads only in bam 1" % len(reads1 - reads2)
    print "%d reads only in bam 2" % len(reads2 - reads1)

    if args.d == 1 and args.o:
        for rid in reads1 - reads2:
            args.o.write('%s\n' % r2aln1[rid].format())
    elif args.d == 2:
        for rid in reads2 - reads1:
            args.o.write('%s\n' % r2aln2[rid].format())

    return

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description="bam utilities")

   subparsers = parser.add_subparsers()

   # bam
   parser_bam_region = subparsers.add_parser("bam_region", help=""" parse bam file for read file """)
   parser_bam_region.add_argument('-reg', help="region, e.g., chr12:34567-45678")
   parser_bam_region.add_argument("bam_fn", help="bam file name")
   parser_bam_region.set_defaults(func=main_bam_region)

   parser_bam_reads = subparsers.add_parser("bam_reads", help=""" parse bam file with given read names """)
   parser_bam_reads.add_argument('-reads', type = argparse.FileType('r'), default='-', help='read file name')
   parser_bam_reads.add_argument('bam_fn', help="bam file name")
   parser_bam_reads.set_defaults(func=main_bam_reads)

   # compare alignments
   parser_compare = subparsers.add_parser("compare", help="""compare two bam files""")
   parser_compare.add_argument('-r1', type = argparse.FileType('r'), default='-', help="read file 1")
   parser_compare.add_argument('-r2', type = argparse.FileType('r'), default='-', help="read file 2")
   parser_compare.add_argument('-d', type=int, default=0, help="which side of difference to print, 1 for 1-2, 2 for 2-1");
   parser_compare.add_argument('-o', help='output difference file', type=argparse.FileType('a'), default=None)
   parser_compare.set_defaults(func=main_compare)

   # filter read file
   parser_filter = subparsers.add_parser("filter", help="""filter reads""")
   parser_filter.add_argument('-s', help="filter statement. E.g., 'not r.is_duplicate' or 'not r.is_proper_pair' or 'is_unmapped'")
   parser_filter.add_argument('read_file', type = argparse.FileType('r'), default='-', help='read file name')
   parser_filter.set_defaults(func=main_filter)

   # tabulate read file
   parser_tabulate = subparsers.add_parser("tabulate", help="""tabulate read file""")
   parser_tabulate.add_argument('read_file', type = argparse.FileType('r'), default='-', help='read file name')
   parser_tabulate.add_argument('-p', default='', help='columns to display in the table')
   parser_tabulate.set_defaults(func=main_tabulate)

   args = parser.parse_args()
   args.func(args)
