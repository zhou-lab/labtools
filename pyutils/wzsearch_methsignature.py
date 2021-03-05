#!/usr/bin/env python

from __future__ import division # make sure division yields floating
import sys, re
import argparse
import numpy as np


def main(args):

    samples = []
    group = []
    for line in args.g:
        fields = line.strip().split('\t')
        samples.append(fields[0])
        group.append(fields[1])

    if args.i != '':
        groups = args.i.split(',')
    else:
        groups = [g for g in sorted(set(group)) if g not in args.x.split(',')]

    if args.u is not None:
        is_test = {}
        for g in groups:
            if g in args.u.split(','):
                is_test[g] = True
            else:
                is_test[g] = False

    print("chrm", "beg", "end", "test", "\t".join(groups), sep="\t")
    for ii, line in enumerate(args.m):
        if (ii == 0):
            continue
        group2meths = {}
        fields = line.strip().split('\t')
        for i, m in enumerate(fields):
            if (i<3) or (m == "NA"):
                continue
            if (group[i-3] in group2meths):
                group2meths[group[i-3]].append(float(m))
            else:
                group2meths[group[i-3]] = [float(m)]

        print("%s\t" % "\t".join(fields[:3]), end="")
        # test
        num_test_true_u = 0
        num_test_true = 0
        num_test_false_m = 0
        num_test_false = 0
        for j,g in enumerate(groups):
            if is_test[g]:
                if g in group2meths:
                    if np.mean(group2meths[g]) < args.fu:
                        num_test_true_u += 1
                    num_test_true += 1
            elif g in group2meths:
                if np.mean(group2meths[g]) > args.fm:
                    num_test_false_m += 1
                num_test_false += 1

        if (num_test_true > 0 and num_test_false > 0 and
            num_test_true_u / num_test_true >= args.pu and
            num_test_false_m / num_test_false >= args.pm):
            print("1\t", end="")
        else:
            print("0\t", end="")

        # output data
        for j,g in enumerate(groups):
            if (j!=0):
                print("\t", end="")
            if (g in group2meths):
                print("%1.2f" % np.mean(group2meths[g]),end="")
            else:
                print("NA",end="")

                
        print()

if __name__ == '__main__':

    parser = argparse.ArgumentParser('search methylation signature')
    parser.add_argument('-g', help='grouping file', type = argparse.FileType('r'), default='-')
    parser.add_argument('-m', help='methylation matrix from tbk', type = argparse.FileType('r'), default='-')
    parser.add_argument('-u', help='signature unmethylation group', default = '')
    parser.add_argument('-x', help='group to ignore', default = '')
    parser.add_argument('-i', help='group to include', default = '')
    parser.add_argument('--fu', help='def of unmethylated', type=float, default = 0.2)
    parser.add_argument('--fm', help='def of methylated', type=float, default = 0.8)
    parser.add_argument('--pu', help='minimum percentage of unmethylated', type=float, default = 1)
    parser.add_argument('--pm', help='minimum percentage of methylated', type=float, default = 1)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
