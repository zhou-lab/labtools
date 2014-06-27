#!/usr/bin/env python

# /RIS/HPC_apps/AMD/python/python-2.7.2/bin/python

import sys
import argparse

def main(args):

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

    for i, line in enumerate(args.table):
        pair = line.strip().split(args.delim)
        print '\t'.join([pair[_] for _ in indices if _ < len(pair)])

    return

if __name__ == '__main__':
    
    
    psr = argparse.ArgumentParser(description='Usage: wzreorder.py -c 1,2,4-6,3 table')
    
    psr.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    psr.add_argument('--delim', default="\t", 
                     help="table delimiter [\\t]")
    psr.add_argument('-c', default=None, help="columns to be printed in the output, 1-based. E.g., -c 1,3-4 [None]")

    args = psr.parse_args()
    main(args)
