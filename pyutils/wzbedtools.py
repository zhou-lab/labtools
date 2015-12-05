#!/usr/bin/env python
""" this is meant to complement bedtools functionality """

import argparse
from collections import deque
import sys

class Bed1(object):

    def __init__(self, line=None):
        if line is not None:
            self.fields = line.strip('\n').split('\t')
            self.chrm = self.fields[0]
            self.beg = int(self.fields[1])
            self.end = int(self.fields[2])
            
    def __repr__(self):
        return '\t'.join([self.chrm, str(self.beg), str(self.end)] + self.fields)

    def dup(self):

        d = Bed1()
        d.chrm = self.chrm
        d.beg  = self.beg
        d.end  = self.end
        d.fields = self.fields[:]

        return d
    
def subtract(a, b):
    """ this will return [] if subtract to nothing """
    beg_max = max(a.beg, b.beg)
    end_min = min(a.end, b.end)
    if a.end < b.beg or b.end < a.beg: # a, b non-overlap
        return [a]
    elif b.beg <= a.beg and b.end >= a.end: # b cover a
        return []
    else:
        ass = []
        if b.beg > a.beg:
            a2 = Bed1()
            a2.chrm = a.chrm
            a2.beg  = a.beg
            a2.end  = b.beg
            a2.fields = a.fields[:]
            ass.append(a2)

        if b.end < a.end:
            a3 = Bed1()
            a3.chrm = a.chrm
            a3.beg  = b.end
            a3.end  = a.end
            a3.fields = a.fields[:]
            ass.append(a3)

        return ass

class Bed(object):

    def __init__(self, fh):
        self.fh = fh
        line = next(self.fh, None)
        if line is None:
            self.nxt = None
        else:
            self.nxt = Bed1(line)

    def read1(self):
        if self.nxt is None:
            return
        line = next(self.fh, None)
        if line is None:
            self.nxt = None
        else:
            self.nxt = Bed1(line)


def update_b1s(a, bbed, b1s):

    _b1s = []
    for b1 in b1s:
        if b1.chrm < a.chrm or (b1.chrm == a.chrm and b1.end < a.beg):
            continue
        _b1s.append(b1)
    b1s = _b1s

    while True:
        b1 = bbed.nxt
        if bbed.nxt is None:
            break
        if bbed.nxt.chrm > a.chrm or (bbed.nxt.chrm == a.chrm and bbed.nxt.beg > a.end):
            break
        if bbed.nxt.chrm < a.chrm or (bbed.nxt.chrm == a.chrm and bbed.nxt.end < a.beg):
            bbed.read1()
            continue
        b1s.append(bbed.nxt)
        bbed.read1()

    return b1s

def main_subtract(args):

    """ subtract A by B but appending original record of A
    both A and B are sorted
    """

    a = Bed(args.a)
    b = Bed(args.b)
    
    b1s = []
    while True:
        if a.nxt is None:
            break
        a1 = a.nxt.dup()
        b1s = update_b1s(a1, b, b1s) # element to be subtracted
        a1s = [a1]
        for b1 in b1s:          # subtract each element
            _a1s = []
            for _a1 in a1s:
                _a1s.extend(subtract(_a1, b1))
            a1s = _a1s
            # print len(_a1s)

        if args.c:
            cnts=0
            for a1 in a1s:
                cnts+=a1.end-a1.beg
            print "%s\t%d" % ('\t'.join(a.nxt.fields), cnts)
        else:
            for a1 in a1s:
                if args.v:
                    # report all the elements subtracted
                    print '%s\t%s' % (a1, ','.join(['%d-%d' % (b1.beg,b1.end) for b1 in b1s]))
                else:
                    print a1
        a.read1()

def main_merge(args):

    """ merge A and B """

    a = Bed(args.a)
    b = Bed(args.b)

    b1s = []
    while True:
        if a.nxt is None:
            break
        a1 = a.nxt
        b1s = update_b1s(a1, b, b1s)
        a1s = [a1]
        if len(b1s) == 0:
            continue
        assert(len(b1s)==1)
        print a1, b1
        b1 = b1s[0]
        assert(a1.chrm == b1.chrm and a1.beg == b1.beg and a1.end == b1.end)
        b1_beta = float(b1.fields[4])
        b1_covg = int(b1.fields[5])
        a1_beta = float(a1.fields[4])
        a1_covg = int(a1.fields[5])
        covg = a1_covg + b1_covg
        beta = (a1_beta * a1_covg + b1_beta * b1_covg) / float(covg)

        print '\t'.join(map(str, [a1.chrm, beta, covg, a1_beta, a1_covg, b1_beta, b1.covg]))

def main_sample(args):

    """ sample lines from a file (usually a bed file) """
    lines = []
    for line in args.i:
        lines.append(line)

    import random
    for line in random.sample(lines, args.n):
        print line

def main_space(args):
    """ calculate space distance between (sorted) bed file records """
    p = None
    pc = None
    a = Bed(args.i)

    while True:
        if a.nxt is None:
            break

        a1 = a.nxt
        if a1.chrm != pc:
            pc = a1.chrm
        else:
            args.o.write('%d\t%s\n' % ((a1.beg - p),str(a1)))
        p = a1.end
        a.read1()

def main_einclude(args):
    """ include from top of the list but keeping included distant, no sorting requirement """

    included = []
    a = Bed(args.i)

    while True:

        if a.nxt is None:
            break

        isnew = True
        for b in included:
            if b.chrm == a.nxt.chrm and b.end + args.d > a.nxt.beg and a.nxt.end + args.d > b.beg:
                isnew = False
                break
        if isnew:
            included.append(a.nxt)
            if args.n > 0 and len(included) >= args.n:
                break
        a.read1()

    for b in included:
        args.o.write(str(b)+'\n')

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Manipulate bed files')
    subparsers = parser.add_subparsers()

    parser_subtract = subparsers.add_parser('subtract', help='subtract A by B, append original record of A')
    parser_subtract.add_argument('-a', required=True, type=argparse.FileType('r'), help='file to subtract from')
    parser_subtract.add_argument('-b', required=True, type=argparse.FileType('r'), help='file to subtract')
    parser_subtract.add_argument('-c', action='store_true', help='report base count only')
    parser_subtract.add_argument('-v', action='store_true', help='verbose')
    parser_subtract.set_defaults(func=main_subtract)

    parser_merge = subparsers.add_parser('merge', help='merge two bed files with beta-coverage information')
    parser_merge.add_argument('-a', required=True, type=argparse.FileType('r'), help='file A')
    parser_merge.add_argument('-b', required=True, type=argparse.FileType('r'), help='file B')
    parser_merge.set_defaults(func=main_merge)

    parser_sample = subparsers.add_parser('sample', help='sample n lines from bed file')
    parser_sample.add_argument('-i', type=argparse.FileType('r'), default='-', help='bed file to sample from')
    parser_sample.add_argument('-n', type=int, required=True, help='number of records to sample')
    parser_sample.set_defaults(func=main_sample)

    
    parser_space = subparsers.add_parser('space', help='calculate space between (sorted) bed records')
    parser_space.add_argument('-i', type=argparse.FileType('r'), default='-', help='input bed file')
    parser_space.add_argument('-o', help='output', default=sys.stdout)
    parser_space.set_defaults(func=main_space)
    
    parser_einclude = subparsers.add_parser('einclude', help='sequential even include')
    parser_einclude.add_argument('-i', type=argparse.FileType('r'), default='-', help='input table')
    parser_einclude.add_argument('-n', type=int, default=-1, help='number of inclusion (-1 for infinity)')
    parser_einclude.add_argument('-d', type=int, default=100, help='spacing (100bp)')
    parser_einclude.add_argument('-o', help='output', default=sys.stdout)
    parser_einclude.set_defaults(func=main_einclude)
    
    args = parser.parse_args()

    try:
        args.func(args)
    except IOError as e:
        sys.exit()
