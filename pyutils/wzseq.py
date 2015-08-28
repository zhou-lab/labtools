#!/usr/bin/env python
import sys
import argparse
import numpy as np
import scipy.stats as stats
import math
import collections as colls
from wzcore import *

class Indices:

    def __init__(self):
        self.spans = []

    def extend(self, start, end):
        self.spans.append((start, end))

    def extract(self, lst):
        result = []
        for start, end in self.spans:
            if not end:
                end = len(lst)
            result.extend([lst[_] for _ in xrange(start, end)])

        return result

def parse_indices(indstr):
    indices = Indices()
    if not indstr: return indices
    rgs = indstr.split(',')
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = None
            indices.extend(int(pair[0])-1 if pair[0] else 0,
                           int(pair[1]) if pair[1] else None)
        else:
            indices.extend(int(rg)-1, int(rg))

    return indices

def main_printc(args):

    import faidx
    genome = faidx.RefGenome(args.i)
    out = open(args.o, 'w') if args.o is not None else sys.stdout
    for c in genome.faidx:
        if args.v:
            err_print(c)
        gseq = genome.fetch_chrmseq(c)
        for i in xrange(len(gseq)-2):

            # for CG, print the position of CG (2-bases)
            if gseq[i] == 'C' and gseq[i+1] == 'G':
                tprint([c, i, i+2, 'CG', '+', 'CG'],out)

            # for CHG, print the position of CHG (3-bases)
            if gseq[i] == 'C' and gseq[i+1] != 'G' and gseq[i+2] == 'G':
                if gseq[i+1] != 'N':
                    tprint([c, i, i+3, 'CHG', '+', gseq[i:i+3]],out)

            # for CHH, print the position of C
            if gseq[i] == 'C' and gseq[i+1] != 'G' and gseq[i+2] != 'G':
                if gseq[i+1] != 'N' and gseq[i+2] != 'G':
                    tprint([c, i, i+1, 'CHH', '+', gseq[i:i+3]],out)
            if gseq[i] != 'G' and gseq[i+1] != 'C' and gseq[i+2] == 'G':
                if gseq[i] != 'N' and gseq[i+1] != 'N':
                    tprint([c, i+2, i+3, 'CHH', '-', reverse_complement(gseq[i:i+3])],out)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='calculate sequence')
    subparsers = parser.add_subparsers()

    psr_printc = subparsers.add_parser("printc", help=""" print C """)
    psr_printc.add_argument('-i', required=True, help="sequence format in fasta")
    psr_printc.add_argument('-o', default=None, help='output file (default stdout)')
    psr_printc.add_argument('-v', action='store_true', help='print status')
    psr_printc.set_defaults(func=main_printc)

    args = parser.parse_args()
    try:
        args.func(args)
    except IOError as e:
        sys.exit()
