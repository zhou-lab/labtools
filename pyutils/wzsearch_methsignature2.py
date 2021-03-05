#!/usr/bin/env python

from __future__ import division # make sure division yields floating
import sys, re
import argparse
import numpy as np


def main(args):

    is_include = {}
    is_test = {}
    for ii, line in enumerate(args.m):
        fields = line.strip().split('\t')
        if (ii == 0):
            groups = fields[6:]
            for g in groups:
                is_include[g] = True
                if g in args.u.split(','):
                    is_test[g] = True
                else:
                    is_test[g] = False

            if args.i != '':
                for g in groups:
                    if g in args.i.split(','):
                        is_include[g] = True
                    else:
                        is_include[g] = False

            if args.x != '':
                for g in args.x.split(','):
                    is_include[g] = False
            print(line.strip(), end='\n')

            continue

        min_v = np.min([float(fields[k]) for k in range(6, len(fields)) if fields[k] != "NA" and is_include[groups[k-6]]])
        n_notna = len([k for k in range(6, len(fields)) if fields[k] != "NA" and is_include[groups[k-6]]])
        n_included = len([k for k in range(6, len(fields)) if is_include[groups[k-6]]])
        if n_notna / n_included < args.pn:
            continue

        cnt_in_pass = 0.0
        cnt_in_all = 0.0
        cnt_out_pass = 0.0
        cnt_out_all = 0.0
        for jj in range(6, len(fields)):
            g = groups[jj-6]
            if (not is_include[g]) or fields[jj] == "NA":
                continue
            if is_test[g]:
                cnt_in_all += 1
                if float(fields[jj])-min_v < args.fu:
                    cnt_in_pass += 1
            else:
                cnt_out_all += 1
                # print(min_v, fields[jj], args.fm)
                if float(fields[jj])-min_v > args.fm:
                    cnt_out_pass += 1
        
        # print(cnt_in_pass, cnt_in_all, cnt_out_pass, cnt_out_all)
        if (cnt_in_all == 0 or cnt_out_all == 0):
            continue
        
        if (cnt_in_pass / cnt_in_all >= args.pu) and (cnt_out_pass / cnt_out_all >= args.pm):
            print(line.strip(), end='\n')

if __name__ == '__main__':

    parser = argparse.ArgumentParser('search methylation signature')
    parser.add_argument('-m', help='methylation matrix', type = argparse.FileType('r'), default='-')
    parser.add_argument('-u', help='signature unmethylation group', default = '')
    parser.add_argument('-x', help='group to ignore', default = '')
    parser.add_argument('-i', help='group to include', default = '')
    parser.add_argument('--fu', help='def of unmethylated', type=float, default = 0.2)
    parser.add_argument('--fm', help='def of methylated', type=float, default = 0.7)
    parser.add_argument('--pu', help='minimum percentage of unmethylated', type=float, default = 1)
    parser.add_argument('--pm', help='minimum percentage of methylated', type=float, default = 1)
    parser.add_argument('--pn', help='minimum fraction of non-na', type=float, default=0.8)
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
