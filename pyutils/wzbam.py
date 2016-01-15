#!/usr/bin/env python

import sys, re
import pysam
import argparse
import faidx
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

            self.mpos = x.mpos + 1
            self.keys.append('mpos')

            if x.cigar:
                self.cigar = x.cigar
                self.keys.append('cigar')
                self.cgstr = x.cigarstring
                self.keys.append('cgstr')

            self.mapq = x.mapq
            self.keys.append('mapq')

            self.reverse = x.is_reverse
            self.keys.append('is_reverse')

            self.seq = x.seq
            self.keys.append('seq')

            self.qual = x.qual
            self.keys.append('qual')

            self.unmapped = True if x.is_unmapped else False
            self.keys.append('unmapped')

            self.duplicate = True if x.is_duplicate else False
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


    def getbase(self, tpos):

        if self.pos > tpos:
            return (4, None, None)
            
        rpos = self.pos
        qpos = 0
        for op, clen in self.cigar:
            if op == 0:         # match
                if rpos + clen > tpos:
                    qpos += tpos - rpos
                    return ((0, self.seq[qpos], self.qual[qpos]))
                else:
                    qpos += clen
                    rpos += clen
            elif op == 1:       # insertion
                if rpos == tpos:
                    return ((1, None, self.qual[qpos]))
                else:
                    qpos += clen
            elif op == 2:       # deletion, quality is the base before
                if rpos + clen > tpos:
                    return ((2, None, self.qual[qpos]))
                else:
                    rpos += clen
            elif op == 4:
                qpos += clen
            else:
                raise Exception("unknown cigar: %d" % op)

        return (4, None, None)

    def calend(self):
        end = self.pos
        for op, clen in self.cigar:
            if op == 0:         # match
                end += clen
            elif op == 1:       # insertion
                pass
            elif op == 2:       # deletion, quality is the base before
                end += clen
            elif op == 4 or op == 5:
                pass
            else:
                raise Exception("unknown cigar: %d" % op)
        return end

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

    for r in parse_readfile(args.read_file):
        if eval(args.s):
            print r.format()

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

def main_pileupone(args):

    """ pileup and investigate base qualities at one position """
    m = re.match(r'(\S+):(\d+)', args.t)
    chrm = m.group(1)
    pos = int(m.group(2))
    bam = pysam.Samfile(args.bam_fn)
    for x in bam.fetch(reference=chrm, start=pos, end=pos+1):
        read = Read(x=x, xob=bam)
        op, base, qual = read.getbase(pos)
        if op == 0:
            read.tbase = base
            read.tqual = ord(qual) - 33
            read.keys.append('tbase')
            read.keys.append('tqual')
            print read.format()
            # print "%s\t%s\t%d" % (read.qname, base, ord(qual)-33)

    return

def main_bed6(args):

    bam = pysam.Samfile(args.bam)
    out = open(args.o,"w") if args.o else sys.stdout
    for x in bam.fetch():
        read = Read(x=x, xob=bam)
        if read.unmapped:
            continue
        if read.duplicate:
            continue
        if read.mapq < args.q:
            continue
        out.write('%s\t%d\t%d\t%s\t1\t%s\n' % (read.chrm, read.pos, read.calend(), read.qname, '-' if read.reverse else '+'))

    return
    
def main_editing(args):

    """ pileup editing """

    import subprocess
    proc = subprocess.Popen(['samtools', 'mpileup', '-f', args.ref,
                             '--min-MQ', '30', '--min-BQ', '30', args.bam], stdout=subprocess.PIPE)
    while True:
        line = proc.stdout.readline()
        if not line:
            break
        fields = line.strip('\n').split('\t')
        chrm = fields[0]
        loc = int(fields[1])
        ref = fields[2].upper()
        # cov = int(fields[3])
        bases = fields[4]
        strands = fields[5]

        if ref != 'A' and ref != 'T':
            continue

        i = 0
        n = len(bases)
        bases2 = []
        while True:
            if i >= n:
                break
            b = bases[i]
            if b == '-':        # deletion
                num = ''
                while i+1<n and bases[i+1].isdigit():
                    i += 1
                    num += bases[i]
                num = int(num)
                i += num
            elif b == '+':      # insertion
                num = ''
                while i+1<n and bases[i+1].isdigit():
                    i += 1
                    num += bases[i]
                num = int(num)
                i += num
            elif b == '^':
                i += 1
            elif b != '$':
                bases2.append(b)
            i += 1

        if len(bases2) != len(strands):
            print chrm, loc
            print bases
            print bases2, len(bases2)
            print strands, len(strands)
            sys.exit(1)

        # this works on the stranded protocol
        # for unstranded protocol, one needs
        # to consider T on the positive strand
        nA = 0
        nG = 0
        for b,s in zip(bases2, strands):

            if args.stranded:
                if ref == 'A':
                    if b == '.':
                        nA += 1
                    if b == 'G':
                        nG += 1
                if ref == 'T':
                    if b == ',':
                        nA += 1
                    if b == 'c':
                        nG += 1
            else:
                if ref == 'A':
                    if b == '.' or b == ',':
                        nA += 1
                    if b == 'g' or b == 'G':
                        nG += 1
                if ref == 'T':
                    if b == '.' or b == ',':
                        nA += 1
                    if b == 'c' or b == 'C':
                        nG += 1

        if nA + nG > 0:
            print '%s\t%d\t%d\t%s\t%d\t%d\t%s' % (chrm, loc-1, loc, ref, nA, nG, bases)
    
if __name__ == '__main__':
   parser = argparse.ArgumentParser(description="bam utilities")

   subparsers = parser.add_subparsers()

   # bam_region
   parser_bam_region = subparsers.add_parser("bam_region", help=""" parse bam file for read file """)
   parser_bam_region.add_argument('-reg', help="region, e.g., chr12:34567-45678")
   parser_bam_region.add_argument("bam_fn", help="bam file name")
   parser_bam_region.set_defaults(func=main_bam_region)

   # pileup one position in a bam file
   parser_pileupone = subparsers.add_parser("pileupone", help=""" pileup one position in a bam file """)
   parser_pileupone.add_argument('-t', help="target position. e.g. chr1:1235342")
   parser_pileupone.add_argument('bam_fn', help="bam file name")
   parser_pileupone.set_defaults(func=main_pileupone)

   # bam_reads
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

   # to as the input for BCP, GEM etc.
   parser_bed6 = subparsers.add_parser('bed6', help='convert bam to bed6 format')
   parser_bed6.add_argument('-bam', required=True, help='input bam')
   parser_bed6.add_argument('-o', help='output', default=None)
   parser_bed6.add_argument('-q', type=int, default=20, help='minimum mapq')
   parser_bed6.set_defaults(func=main_bed6)

   # pileup RNA-editting
   parser_editing = subparsers.add_parser('editing', help='pileup editing from bam file')
   parser_editing.add_argument('-bam', required=True, help='input bam')
   parser_editing.add_argument('-ref', required=True, help='fai-indexed reference')
   parser_editing.add_argument('-stranded', action='store_true', help="consider only A>G on positive strand and T>C on negative strand (otherwise consider both A>G and T>C on both strand)")
   parser_editing.set_defaults(func=main_editing)

   args = parser.parse_args()
   args.func(args)
