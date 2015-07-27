#!/usr/bin/env python
import argparse
import faidx
import sys
from collections import deque

def main_scan(args):
    genome = faidx.RefGenome(args.g)
    chrms = sorted(genome.faidx.keys())
    for chrm in chrms:
        seq = genome.fetch_chrmseq(chrm).upper()
        sys.stderr.write('%s\n' % chrm)
        window300 = deque(seq[:300])
        n = {}
        for b in 'ATCGN':
            n[b] = window300.count(b)
        nCG = window300.count('CG')
        nGC = window300.count('GC')
        for i in xrange(300,len(seq)):
            n_all = n['A']+n['T']+n['G']+n['C']
            if n_all == 300:
                GCcont = (n['C']+n['G']) / float(n_all)
                if GCcont > 0.5:
                    CGodds = nCG / float(GCcont*GCcont/4.0*n_all)
                    if CGodds > 0.65:
                        print '%s\t%d\t%d' % (chrm, i-300, i)

            b1 = window300.popleft()
            b2 = seq[i]
            n[b1] -= 1
            n[b2] += 1
            if b1 == 'C' and window300[0] == 'G':
                nCG -= 1
            if window300[-1] == 'C' and b2 == 'G':
                nCG += 1
            if b1 == 'G' and window300[0] == 'C':
                nGC -= 1
            if window300[-1] == 'G' and b2 == 'C':
                nGC += 1
            window300.append(b2)

def main_max(args):
    genome = faidx.RefGenome(args.g)
    for line in args.i:
        chrm, beg, end = line.strip().split('\t')
        seq = genome.fetch_sequence(chrm, int(beg)+1, int(end)).upper()
        maxCGodds = -99999
        beg = int(beg)
        end = int(end)
        for i in xrange(len(seq)-299):
            window = seq[i:i+300]
            nCG = window.count('CG')
            nC = window.count('C')
            nG = window.count('G')
            nA = window.count('A')
            nT = window.count('T')
            n_all = nC + nG + nA + nT
            if n_all == 0:
                continue
            GCcont = (nC+nG) / float(n_all)
            if GCcont <= 0.5:
                continue
            CGodds = float(nCG) / (GCcont*GCcont/4.0*n_all)
            if CGodds > maxCGodds:
                maxCGodds = CGodds
                maxbeg = beg + i
                maxend = maxbeg + 300
                maxnCG = nCG
                maxGCcont = GCcont
                maxnC = nC
                maxnG = nG
                maxnA = nA
                maxnT = nT

        if maxCGodds < 0:
            print chrm, beg, end, seq, len(seq), seq.count('C'), seq.count('G'), seq.count('T'), seq.count('A'), seq.count('CG'), seq.count('GC')
            break
            
        print '%s\t%d\t%d\tCG:%d\tGCcontent:%1.3f\t%d\t%d\t%d\t%d' % (
            chrm, maxbeg, maxend, maxnCG, maxGCcont, maxnA, maxnC, maxnG, maxnT)

def main_methlevelaverage(args):

    import tabix
    import numpy as np
    p = tabix.open(args.p)
    for line in open(args.c):
        loc = line.strip('\n').split('\t')
        chrm = loc[0]
        beg = int(loc[1]) + 1
        end = int(loc[2])
        betas = []
        for r in p.query(chrm, beg, end):
            if r[5].endswith('CG') and r[7] != '.':
                retn = int(r[7])
                conv = int(r[8])
                if retn + conv == 0:
                    continue
                beta = retn / float(retn+conv)
                betas.append(beta)
        if len(betas) == 0:
            mbetas = 'NA'
        else:
            mbetas = '%1.2f' % np.mean(betas)
        print '%s\t%d\t%s' % ('\t'.join(loc), len(betas), mbetas)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CGI')
    subparsers = parser.add_subparsers()

    parser_scan = subparsers.add_parser('scan', help='scan the genome, should couple with | bedtools merge -i - -d 1000')
    parser_scan.add_argument('-g', help='genome', default='/home/wzhou/genomes/mm10/mm10.fa')
    parser_scan.set_defaults(func=main_scan)

    parser_max = subparsers.add_parser('max', help='find window that maximizes CpG/GpC')
    parser_max.add_argument('-g', help='genome', default='/home/wzhou/genomes/mm10/mm10.fa')
    parser_max.add_argument('-i', type = argparse.FileType('r'), default='-', help='input table')
    parser_max.set_defaults(func=main_max)

    parser_met = subparsers.add_parser('methlevelaverage', help='find average methylation level')
    parser_met.add_argument('-p', required=True, help='pileup file')
    parser_met.add_argument('-c', required=True, help='bed file')
    parser_met.set_defaults(func=main_methlevelaverage)
    
    args = parser.parse_args()
    args.func(args)


    
