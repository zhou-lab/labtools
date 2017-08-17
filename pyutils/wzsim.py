#!/usr/bin/env python
""" simulate reads """

import argparse
import faidx

def complement(base):

    return {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N'
    }[base]

def reverse_complement(seq):

    return ''.join([complement(base) for base in reversed(seq)])

def depos(rpos, chrmoffsets):

    for chrm, off, l in chrmoffsets:
        if off + l > rpos:
            break
    return chrm, rpos - off

def main_mutation(args):

    genome = faidx.RefGenome(args.genome)
    chromlens = []
    fmut = open(args.prefix + '.mutations.tsv', 'w')
    for chrm in genome.keys():         # process each chromosome
        chrmseq = list(genonme.fetch_chrmseq(chrm))
        original = range(1,len(chrmseq)+1)

        nmut = args.m * len(chrmseq)
        original

        # introduce mutations
        save(original, open(args.prefix + '.' + chrm + '.original.pkl', 'w'))

        # methylation
        save(methyl, open(args.prefix + '.' + chrm + '.methyl.pkl', 'w'))

        for mut in mutations:
            fmut.write('\t'.join(map(str, mut)))

        # chromseq
        save(chrmseq, open(args.prefix + '.' + chrm + '.chrmseq.pkl', 'w'))

        chromlens.append(chrm, len(chrmseq))

    close(fmut)
    save(chromlens, open(args.prefix + '.meta_chroms.pkl', 'w'))

def main_read(args):

    isize_mean, isize_sd = map(float, ','.split(args.I))
    chrmoffsets, genome, original, methyl = load(args.genome)

    genome random.rand()
    chrm, pos = depos(int(random.rand()*genomsize), chrmoffsets)
    # sample isize normal distribution
    seq = faidx.fetch_sequence(chrm, pos, pos + isize)

    # load mutated genome and methylation rate

    # bisulfite conversion and mutation happen here

    if random.rand() < 0.5:
        seq = reverse_complement(seq)
        strand = '-'
    else:
        strand = '+'
    
    qname1 = '%s:%d:%s:%d' % (chrm, original[pos], strand, isize)
    seq[1:args.l], reverse_complement(seq)[1:args.l]
    
    readname = 


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Read Simulator')
    subparsers = parser.add_subparsers()

    parser_mutation = subparsers.add_parser('mutation', help='simulate mutation and methylation')
    parser_mutation.add_argument('genome', help = '.fa genome to start with')
    parser_mutation.add_argument('prefix', help = 'output prefix')
    parser_mutation.add_argument('-m', default = 1e-3, help = 'mutation rate')
    parser_mutation.add_argument('-i', default = 1e-5, help = 'indel rate')
    parser_mutation.add_argument('-c', default = 0.6, help = 'CpG methylation rate')
    parser_mutation.add_argument('-d', default = 1e-10, help = 'CpH methylation rate')
    parser_mutation.set_defaults(func = main_mutation)

    parser_read = subparsers.add_parser('read', help='simulate short read')
    parser_read.add_argument('genome', help='mutated fasta file, should have .methyl,.original and .mutation')
    parser_read.add_argument('-I', default="200,50", help='insert size distribution (mean,sd)')
    parser_read.add_argument('-l', default=100, help='read length')
    parser_read.add_argument('-p', default=None, help='Output prefix, read output to _1/2.fq')
    parser_read.set_defaults(func = main_read)


