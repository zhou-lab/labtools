#!/usr/bin/env python
import argparse
import tabix

def main_vcf2hm450k(args):

    probes = []
    for line in open(args.p):
        fields = line.split('\t')
        chrm = fields[0]
        beg = int(fields[1])
        end = int(fields[2])
        pname = fields[3]
        if pname.startswith('rs'):
            continue
        probes.append((pname, chrm, beg))

    probes.sort()
    fh = tabix.open(args.f)
    for p in probes:
        print p

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='convert WGBS vcf to hm450k')
    
    parser.add_argument('-f', help='vcf file, must be tabix indexed')
    parser.add_argument('-p', default='/data/largeS2/pl-bs/data/450k_probes')
    parser.set_defaults(func=main_vcf2hm450k)

    args = parser.parse_args()
    args.func(args)
