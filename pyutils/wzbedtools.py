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

def subtract(a, b):
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

    a = Bed(args.a)
    b = Bed(args.b)
    
    b1s = []
    while True:
        if a.nxt is None:
            break
        a1 = a.nxt
        b1s = update_b1s(a1, b, b1s)
        a1s = [a1]
        for b1 in b1s:
            _a1s = []
            for _a1 in a1s:
                _a1s.extend(subtract(_a1, b1))
            a1s = _a1s

        for a1 in a1s:
            if args.v:
                print '%s\t%s' % (a1, ','.join(['%d-%d' % (b1.beg,b1.end) for b1 in b1s]))
            else:
                print a1
        a.read1()

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Manipulate bed files')
    subparsers = parser.add_subparsers()

    parser_subtract = subparsers.add_parser('subtract', help='subtract A by B, append original record of A')
    parser_subtract.add_argument('-t', help="data table", type = argparse.FileType('r'),
                                 default='-')
    parser_subtract.add_argument('-a', required=True, type=argparse.FileType('r'), help='file to subtract from')
    parser_subtract.add_argument('-b', required=True, type=argparse.FileType('r'), help='file to subtract')
    parser_subtract.add_argument('-v', action='store_true', help='verbose')
    parser_subtract.set_defaults(func=main_subtract)
    
    args = parser.parse_args()
    args.func(args)

