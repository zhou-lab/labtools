#!/usr/bin/env python

import sys
import argparse
import numpy as np
import scipy.stats as stats
import math
import collections as colls
import wzcore

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
    
    indices = parse_indices(args.c)

    if args.skipheader:
        header = args.table.readline()
        headerfields = header.strip().split(args.delim)
        print('stats\t%s' % '\t'.join(indices.extract(headerfields)))

    vals = []
    for i, line in enumerate(args.table):
        fields = line.strip().split(args.delim)
        if i == 0:
            vals = [[] for _ in indices.extract(fields)]
        for j, val in enumerate(indices.extract(fields)):
            if val != 'NA':
                vals[j].append(float(val))


    print('valid_count\t%s' % '\t'.join(map(str, map(len, vals))))
    print('sum\t%s' % '\t'.join(map(str, map(sum, vals))))
    print('median\t%s' % '\t'.join(map(str, map(np.median, vals))))
    print('mean\t%s' % '\t'.join(map(str, map(np.mean, vals))))
    print('max\t%s' % '\t'.join(map(str, map(max, vals))))
    print('min\t%s' % '\t'.join(map(str, map(min, vals))))
    
    # data = zip(*vals)
    # print "sum:\t%s" % '\t'.join([str(sum(d)) for d in data])
    # print "sum: %.2f" % sum(vals) '\t'.join([str(sum(d)) for d in data])
    # print "median:\t%s" % '\t'.join([str(np.median(d)) for d in data])
    # print "max:\t%s" % '\t'.join([str(max(d)) for d in data])
    # print "min:\t%s" % '\t'.join([str(min(d)) for d in data])
    # print "median: %.2f" % np.median(vals)
    # print "max: %.2f" % max(vals)
    # print "min: %.2f" % min(vals)

def main_rowstd(args):

    """ do complex statistics on each row """

    indices = parse_indices(args.c)
    for line in args.table:
        fields = line.strip('\n').split(args.delim)
        vals = [float(_) for _ in indices.extract(fields) if _ != args.x]

        print '%s\t%d\t%s' % ('\t'.join(fields), len(vals), '%1.3f' % np.std(vals) if len(vals)>1 else '-1')

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

def main_rowentropy(args):

    """ compute entropy of each row """

    indices = parse_indices(args.c)
    output = parse_indices(args.o)

    for line in args.table:
        fields = line.strip().split(args.delim)
        print '%s\t%1.2f' % (
            '\t'.join(output.extract(fields)),
            stats.entropy(colls.Counter(list(indices.extract(fields))).values()))

def main_zscore(args):

    """ compute Z-score for each row """
    vals = []
    fieldss = []
    for line in args.table:
        fields = line.strip('\n').split(args.delim)
        fieldss.append(fields)
        if args.x and eval(args.x):
            vals.append(None)
            v = None
        else:
            val = float(fields[args.c-1])
            vals.append(val)
            v = val

    d = [v for v in vals if v]
    vmean = np.mean(d)
    vstd = np.std(d)
    for i, val in enumerate(vals):
        fields = fieldss[i]
        if val is None:
            print '%s\tNA' % ('\t'.join(fields), )
        else:
            print '%s\t%1.3f' % ('\t'.join(fields), (val - vmean) / vstd)

def add_std_options(psr):

    psr.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    psr.add_argument('--delim', default="\t", help="table delimiter [\\t]")
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
    # psr_basic.add_argument('-r', action="store_true", help="row mode, apply to each row in the input table")
    # psr_basic.add_argument('--eldelim', default=",", help="delimiter between elements, only effective under row mode [,]")
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

    p = subparsers.add_parser('rowentropy', help='compute entropy on each row')
    p.add_argument('-c', default=None, help='columns to compute')
    p.add_argument('-o', default='-', help='columns to output')
    add_std_options(p)
    p.set_defaults(func=main_rowentropy)

    p = subparsers.add_parser('rowstd', help='compute statistics on each row')
    p.add_argument('-c', required=True, help='columns to compute')
    p.add_argument('-x', default=None, help='character to exclude, e.g., NA')
    add_std_options(p)
    p.set_defaults(func=main_rowstd)

    p = subparsers.add_parser('zscore', help='calculate the Z-score for each row in the distribution of all rows')
    p.add_argument('-c', type=int, help='column to compute')
    p.add_argument('-x', default=None, help='exclusion criterion')
    add_std_options(p)
    p.set_defaults(func=main_zscore)

    args = parser.parse_args()
    try:
        args.func(args)
    except IOError as e:
        sys.exit()
