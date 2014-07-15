#!/usr/bin/env python

import sys
import argparse
import numpy as np
import scipy.stats as stats
import math

def main_Utest(args):
    
    if args.skipheader:
        args.table.readline()

    if args.p:
        prn_indices = [int(_)-1 for _ in args.p.split(',')]

    for i, line in enumerate(args.table):
        pair = line.strip().split(args.delim)
        c1 = map(float,pair[args.c1-1].split(args.eldelim))
        c2 = map(float,pair[args.c2-1].split(args.eldelim))

        z, p = stats.ranksums(c1, c2)
        if args.onetail:
            if z>0:
                p = 1.0
            else:
                p /= 2.0

        op = []
        if args.p:
            op += [pair[_] for _ in prn_indices]
            
        op.append("%.3g" % p)

        if args.outz:
            op.append("%.3f" % z)

        print "\t".join(op)

    return

def main_basic(args):
    
    if args.skipheader:
        args.table.readline()

    rgs = args.c.split(',')
    indices = []
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = 99
            indices.extend(range(int(pair[0])-1, int(pair[1])))
        else:
            indices.append(int(rg)-1)

    vals = []
    for i, line in enumerate(args.table):
        fields = line.strip().split(args.delim)
        vals.append([float(fields[_]) for _ in indices])

    data = zip(*vals)

    print "sum:\t%s" % '\t'.join([str(sum(d)) for d in data])
    print "sum: %.2f" % sum(vals)
    print "median: %.2f" % np.median(vals)
    print "max: %.2f" % max(vals)
    print "min: %.2f" % min(vals)

def main_overlap(args):
    """ find overlap of two list """

    set1 = set()
    for line in args.t1:
        fields = line.strip()
        set1.add(fields[args.c1])

    set2 = set()
    for line in args.t2:
        fields = line.strip()
        set2.add(fields[args.c2])
    
    print '\n'.join(set1 & set2)

def add_std_options(psr):

    psr.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    psr.add_argument('--delim', default="\t", 
                     help="table delimiter [\\t]")
    psr.add_argument('--skipheader', action='store_true', help='skip header')
    psr.add_argument('-p', default=None, help="columns to be printed in the output, 1-based. E.g., -p 1,3,4 [None]")

    # psr.add_argument('-o', dest='outfn', default="wz.out",
    #                  help='output file name [stdout]')

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='calculator of different statistics')
    subparsers = parser.add_subparsers()
    
    # basic statistics
    psr_basic = subparsers.add_parser("basic", help=""" Basic statistics (mean, median, std, max, min) """)
    psr_basic.add_argument('-c', default=None, help="columns to study, 1-based. E.g., -c 1,3-4 [None]")
    psr_basic.add_argument('-r', action="store_true", help="row mode, apply to each row in the input table")
    psr_basic.add_argument('--eldelim', default=",", help="delimiter between elements, only effective under row mode [,]")
    add_std_options(psr_basic)
    psr_basic.set_defaults(func=main_basic)

    # Utest
    psr_Utest = subparsers.add_parser("Utest", help=""" Mann-Whitney U (rank-sum) test """)
    psr_Utest.add_argument('-c1', type=int, required=True, help="first column to compare")
    psr_Utest.add_argument('-c2', type=int, required=True, help="second column to compare")
    psr_Utest.add_argument('--eldelim', default=",", help="delimiter between elements, only effective under row mode [,]")
    psr_Utest.add_argument('--onetail', action="store_true", help="perform one-tailed analysis, c1<c2.")
    psr_Utest.add_argument('--outz', action="store_true", help="output z-score.")
    add_std_options(psr_Utest)
    psr_Utest.set_defaults(func=main_Utest)


    args = parser.parse_args()
    args.func(args)
