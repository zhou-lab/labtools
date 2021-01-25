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

    groups = sorted(set(group))

    print("chrm", "beg", "end", "\t".join(groups), sep="\t")
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
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
