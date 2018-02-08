#!/usr/bin/env python
""" simulate reads """

import os, re
import argparse
import faidx
import random
import numpy.random
from cPickle import dump, load

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

def main_mutation(args):

    """ generate genome with mutation from a mutation list """
    
    genome = faidx.RefGenome(args.genome)
    fmut = open(args.prefix + '/mutations.tsv', 'w')
    for chrm in genome.faidx.keys():         # process each chromosome
        print(chrm)
        chrmseq = list(genome.fetch_chrmseq(chrm))
        original = range(len(chrmseq))
        deletion = [False]*len(chrmseq)
        insertion = ['']*len(chrmseq)

        # sample SNP
        nmut = int(args.m * len(chrmseq))
        print("sampling %d SNPs" % nmut)
        i = 0
        for coord in random.sample(original, 2*nmut):
            if i >= nmut:
                break
            source = chrmseq[coord]
            if source == 'N':
                continue
            target = random.choice([_ for _ in ['A','C','G','T'] if _ != chrmseq[coord]])
            fmut.write('%s\t%d\t%d\t%s\t%s\tSNP\n' % (chrm, coord, coord+1, source, target))
            chrmseq[coord] = target
            i += 1
        
        # sample insertion
        nIns = int(args.i * len(chrmseq) / 2)
        print("sampling %d insertion" % nIns)
        lambdaIns = 1
        lenIns = numpy.random.poisson(lambdaIns, nIns)+1
        i = 0
        for coord in random.sample(original, 2*nIns):
            if i >= nIns:
                break
            if insertion[coord] != '': # already have insertion
                continue
            source = chrmseq[coord]
            if source == 'N':
                continue
            insertion[coord] = ''.join(numpy.random.choice(['A','T','C','G'],lenIns[i]))
            fmut.write('%s\t%d\t%d\t%s\t%s\tINS\n' % (chrm, coord, coord+1, source, source+insertion[coord]))
            i += 1

        # sample deletion
        nDel = int(args.i * len(chrmseq) / 2)
        print("sampling %d deletions" % nDel)
        lambdaDel = 1
        lenDel = numpy.random.poisson(lambdaDel, nDel)+1
        i = 0
        for coord in random.sample(original, 2*nDel):
            if i >= nDel:
                break

            source = chrmseq[coord]
            if source == 'N':
                continue

            toDelete = True
            for j in xrange(lenDel[i]):
                if insertion[coord+j] != '' or deletion[coord+j]: # avoid mixing insertion with deletion
                    toDelete = False
            if not toDelete:
                continue

            delSeq = ''.join([chrmseq[coord+j+1] for j in xrange(lenDel[i])])
            for j in xrange(lenDel[i]):
                deletion[coord+j+1] = True
            fmut.write('%s\t%d\t%d\t%s\t%s\tDEL\n' % (chrm, coord, coord+lenDel[i]+1, source+delSeq, source))
            i += 1

        chrmseq1 = ''.join([(b+insertion[i]) if insertion[i] != '' else b for i,b in enumerate(chrmseq) if not deletion[i]])
        coords = []
        for i, b in enumerate(chrmseq):
            if not deletion[i]:
                coords.append(i) 
                if insertion[i] != '':
                    coords.extend([-1]*len(insertion[i]))

        # the original coordinates
        dump(coords, open(args.prefix + '/' + chrm + '.coords.pkl', 'w'))

        # chromseq
        # chrmseq0 = ''.join(chrmseq)
        dump(chrmseq1, open(args.prefix + '/' + chrm + '.chrmseq.pkl', 'w'))

    fmut.close()

def main_reads(args):

    """ sample bisulfite reads from a genome with mutation """
    isize_mean, isize_sd = map(float, args.I.split(','))

    chrm2len = {}
    genomelen = 0
    for pklfile in os.listdir(args.genomeDir):
        m = re.match(r'([^.]*).chrmseq.pkl', pklfile)
        if m:
            chrmseq1 = load(open(args.genomeDir+'/'+pklfile))
            chrm2len[m.group(1)] = len(chrmseq1)
            genomelen += len(chrmseq1)
            print ("Priming chromosomes: %s (%d)." % (m.group(1), len(chrmseq1)))

    fout1 = open(args.prefix+"_R1.fastq", 'w')
    fout2 = open(args.prefix+"_R2.fastq", "w")
    gi = 1
    chrms = sorted(chrm2len.keys())
    for chrm in chrms:
        l1 = chrm2len[chrm]
        n1 = int(float(l1) / genomelen * args.n) # number of reads to sample
        print ("Sampling %d reads from chromosome %s" % (n1, chrm))

        chrmseq1 = load(open(args.genomeDir+'/'+chrm+'.chrmseq.pkl'))
        coord1 = load(open(args.genomeDir+'/'+chrm+'.coords.pkl'))
        poses = numpy.random.choice(range(1,(l1-args.l)), n1*2)
        isizes = numpy.random.randn(n1*2)*isize_sd + isize_mean
        bsstrands = ['CT' if numpy.random.rand() > 0.5 else 'GA' for i in xrange(len(poses))]
        ii = 0
        for i, pos in enumerate(poses):

            if chrmseq1[pos] == 'N':
                continue
            
            isize = int(isizes[i])
            if isize < args.l:
                isize = args.l

            # guard agaist edge
            if pos + isize >= l1:
                pos = l1 - isize - 1

            bsseq = []
            if bsstrands[i] == 'CT':
                for b in xrange(pos, pos+isize):
                    if chrmseq1[b] == 'C':
                        if chrmseq1[b+1] == 'G':
                            if numpy.random.rand() < args.c:
                                bsseq.append('C')
                            else:
                                bsseq.append('T')
                        else:
                            if numpy.random.rand() < args.d:
                                bsseq.append('C')
                            else:
                                bsseq.append('T')
                    else:
                        bsseq.append(chrmseq1[b])
            else:
                for b in xrange(pos, pos+isize):
                    if chrmseq1[b] == 'G':
                        if chrmseq1[b-1] == 'C':
                            if numpy.random.rand() < args.c:
                                bsseq.append('G')
                            else:
                                bsseq.append('A')
                        else:
                            if numpy.random.rand() < args.d:
                                bsseq.append('G')
                            else:
                                bsseq.append('A')
                    else:
                        bsseq.append(chrmseq1[b])

            if ii > n1:
                break
            ii += 1
            if bsstrands[i] == 'CT':
                fout1.write('>%d_%s_%d_%d_%s_R1\n' % (gi, chrm, coord1[pos]+1, coord1[pos+isize]+1, bsstrands[i]))
                fout1.write(''.join(bsseq[:args.l])+'\n')
                fout1.write('+\n')
                fout1.write('9'*args.l+'\n')
                fout2.write('>%d_%s_%d_%d_%s_R2\n' % (gi, chrm, coord1[pos]+1, coord1[pos+isize]+1, bsstrands[i]))
                fout2.write(reverse_complement(''.join(bsseq[-args.l:]))+'\n')
                fout2.write('+\n')
                fout2.write('9'*args.l+'\n')
            elif bsstrands[i] == 'GA':
                fout2.write('>%d_%s_%d_%d_%s_R2\n' % (gi, chrm, coord1[pos]+1, coord1[pos+isize]+1, bsstrands[i]))
                fout2.write(''.join(bsseq[:args.l])+'\n')
                fout2.write('+\n')
                fout2.write('9'*args.l+'\n')
                fout1.write('>%d_%s_%d_%d_%s_R1\n' % (gi, chrm, coord1[pos]+1, coord1[pos+isize]+1, bsstrands[i]))
                fout1.write(reverse_complement(''.join(bsseq[-args.l:]))+'\n')
                fout1.write('+\n')
                fout1.write('9'*args.l+'\n')
            gi += 1
                
            
    # chrmoffsets, genome, original, methyl = load(args.genome)

    # # genome random.rand()
    # chrm, pos = depos(int(random.rand()*genomsize), chrmoffsets)
    # # sample isize normal distribution
    # seq = faidx.fetch_sequence(chrm, pos, pos + isize)

    # # load mutated genome and methylation rate

    # # bisulfite conversion and mutation happen here

    # if random.rand() < 0.5:
    #     seq = reverse_complement(seq)
    #     strand = '-'
    # else:
    #     strand = '+'
    
    # qname1 = '%s:%d:%s:%d' % (chrm, original[pos], strand, isize)
    # seq[1:args.l], reverse_complement(seq)[1:args.l]
    
    # # readname = 


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Read Simulator')
    subparsers = parser.add_subparsers()

    parser_mutation = subparsers.add_parser('mutation', help='simulate mutation and methylation')
    parser_mutation.add_argument('genome', help = '.fa genome to start with')
    parser_mutation.add_argument('prefix', help = 'output prefix')
    parser_mutation.add_argument('-m', default = 1e-2, help = 'mutation rate')
    parser_mutation.add_argument('-i', default = 1e-3, help = 'indel rate')
    parser_mutation.set_defaults(func = main_mutation)

    parser_reads = subparsers.add_parser('reads', help='simulate short read')
    parser_reads.add_argument('genomeDir', help='folder containing mutated sequence files, i.e. *.coords.pkl,*.chrmseq.pkl')
    parser_reads.add_argument('prefix', default='output', help='Output prefix, read output to _1/2.fq')
    parser_reads.add_argument('-I', default="200,50", help='insert size distribution (mean,sd)')
    parser_reads.add_argument('-l', default=100, help='read length')
    parser_reads.add_argument('-n', type=int, default=10000, help='number of reads to sample')
    parser_reads.add_argument('-c', default = 0.6, help = 'CpG retention rate')
    parser_reads.add_argument('-d', default = 0.01, help = 'CpH retention rate')
    parser_reads.set_defaults(func = main_reads)

    args = parser.parse_args()
    args.func(args)

