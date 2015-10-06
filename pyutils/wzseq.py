#!/usr/bin/env python
import sys
import argparse
import numpy as np
import scipy.stats as stats
import math
import collections as colls
import re
from wzcore import *


# Heng Li's readfq
# usage
# if __name__ == "__main__":
#     import sys
#     n, slen, qlen = 0, 0, 0
#     for name, seq, qual in readfq(sys.stdin):
#         n += 1
#         slen += len(seq)
#         qlen += qual and len(qual) or 0
#     print n, '\t', slen, '\t', qlen
def readfq(fp): # this is a generator function
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

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

def main_printc(args):

    import faidx
    genome = faidx.RefGenome(args.i)
    out = open(args.o, 'w') if args.o is not None else sys.stdout
    for c in genome.faidx:
        if args.v:
            err_print(c)
        gseq = genome.fetch_chrmseq(c)
        for i in xrange(len(gseq)-2):

            # for CG (symmetric), print the position of CG (2-bases)
            if gseq[i] == 'C' and gseq[i+1] == 'G':
                tprint([c, i, i+2, 'CG', '+', 'CG'],out)

            # CHGs are assymetric
            # for example, CCG is a CHG but its reverse complement CGG is not a CHG
            if gseq[i] == 'C' and gseq[i+1] != 'G' and gseq[i+2] == 'G':
                if gseq[i+1] != 'N':
                    tprint([c, i, i+1, 'CHG', '+', gseq[i:i+3]],out)

            if gseq[i] == 'C' and gseq[i+1] != 'C' and gseq[i+2] == 'G':
                if gseq[i+1] != 'N':
                    tprint([c, i+2, i+3, 'CHG', '-', reverse_complement(gseq[i:i+3])], out)

            # for CHH, print the position of C
            if gseq[i] == 'C' and gseq[i+1] != 'G' and gseq[i+2] != 'G':
                if gseq[i+1] != 'N' and gseq[i+2] != 'N':
                    tprint([c, i, i+1, 'CHH', '+', gseq[i:i+3]],out)
            if gseq[i] != 'C' and gseq[i+1] != 'C' and gseq[i+2] == 'G':
                if gseq[i] != 'N' and gseq[i+1] != 'N':
                    tprint([c, i+2, i+3, 'CHH', '-', reverse_complement(gseq[i:i+3])], out)

def main_comp(args):

    import faidx, re
    genome = faidx.RefGenome(args.i)
    m = re.match(r'([^:]*):(\d+)-(\d+)', args.g)
    chrm = m.group(1)
    beg = int(m.group(2))
    end = int(m.group(3))
    seq = genome.fetch_sequence(chrm, beg, end, uppercase=True)
    print 'n:%d' % (end-beg)
    print 'C:%d' % seq.count('C')
    print 'G:%d' % seq.count('G')
    print 'CG:%d' % seq.count('CG')

def main_orphan(args):

    import faidx
    genome = faidx.RefGenome(args.i)
    out = open(args.o, 'w') if args.o is not None else sys.stdout
    for c in genome.faidx:
        if args.v:
            err_print(c)
        gseq = genome.fetch_chrmseq(c)
        prev = None
        prev_is_good_left = True
        for i in xrange(len(gseq)-2):

            if gseq[i] == 'C' and gseq[i+1] == 'G':
                if prev and prev_is_good_left and i-prev >= args.l:
                    tprint([c, prev, prev+1, '+'], out)
                if prev:
                    prev_is_good_left = i-prev >= args.l
                prev = i
        if prev_is_good_left and prev is not None:
            tprint([c, prev, prev+1, '+'], out)

def main_internal(args):

    """ get all the internal CpG with flanking sequence """
    from Bio import SeqIO
    for seq_record in SeqIO.parse(args.i, "fasta"):
        
        seq = seq_record.seq.upper()
        for i in xrange(args.l, len(seq)-args.l-2):
            if seq[i:i+2] == 'CG':
                print '\t'.join(map(str, [".",i,i+2, '|'.join(seq_record.description.split()), seq[i-args.l:i]+'['+seq[i:i+2]+']'+seq[i+2:i+args.l+2]]))

def main_filter(args):

    """ filter fasta by keyword """
    from Bio import SeqIO
    for seq_record in SeqIO.parse(args.i, "fasta"):

        if re.search(args.k, seq_record.description):
            print seq_record.format("fasta")

def main_cleanhead(args):

    """ replace white space in head of fasta by "_" """

    from Bio import SeqIO
    for seq_record in SeqIO.parse(args.i, "fasta"):

        seq_record.description = '|'.join(seq_record.description.split())
        print seq_record.format("fasta")

def wrap(seq, l = 80):
    s2 = ''
    sl = len(seq)/l
    for i in xrange(sl+1):
        s2 += seq[i*l:(i+1)*l]
        if i != sl:
            s2 += '\n'

    return s2
        
        
def main_consensus(args):

    """ generate consensus sequence from a multiple sequence alignment """

    from Bio import AlignIO
    from collections import Counter
    align = AlignIO.read(args.i, "clustal")

    support = []
    nseqs = len(align)
    consensus = []
    consupport = []
    for i in xrange(align.get_alignment_length()):
        a = align[:,i]
        support.append(nseqs - a.count('-'))
        most_common, num_most_common = Counter([_ for _ in a if _ != '-']).most_common(1)[0]
        consensus.append(most_common)
        consupport.append(num_most_common/float(nseqs))

    seq = []
    k = 1
    j = 0
    for i in xrange(len(consensus)):
        if consupport[i]>args.s:
            if not seq:
                j = i
            seq.append(consensus[i])
        else:
            if seq and len(seq) > args.l:
                print '>%sconsensus%d_%d' % (args.p, k,j)
                print wrap(''.join(seq))
                k += 1
            seq = []
            
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='calculate sequence')
    subparsers = parser.add_subparsers()

    psr_printc = subparsers.add_parser("printc", help=""" print all C in fasta """)
    psr_printc.add_argument('-i', required=True, help="sequence format in fasta")
    psr_printc.add_argument('-o', default=None, help='output file (default stdout)')
    psr_printc.add_argument('-v', action='store_true', help='print status')
    psr_printc.set_defaults(func=main_printc)

    psr_orphan = subparsers.add_parser('orphan', help=""" orphan CpG """)
    psr_orphan.add_argument('-i', required=True, help='genomic sequence format in fasta')
    psr_orphan.add_argument('-o', default=None, help='output file (default stdout)')
    psr_orphan.add_argument('-v', action='store_true', help='print status')
    psr_orphan.add_argument('-l', type=int, default=35, help='minimum flanking distance for defining orphan (default:35bp)')
    psr_orphan.set_defaults(func=main_orphan)

    psr_comp = subparsers.add_parser('comp', help=""" sequence composition """)
    psr_comp.add_argument('-i', required=True, help='genomic sequence format in fasta')
    psr_comp.add_argument('-g', help='region, e.g., chr1:3342-4312341')
    psr_comp.add_argument('-l', help='region in bed files')
    psr_comp.set_defaults(func=main_comp)

    psr_internal = subparsers.add_parser('internal', help=""" internal CpG """)
    psr_internal.add_argument('-l', default=50, type=int, help='flanking length (50bp)')
    psr_internal.add_argument('-i', type=argparse.FileType('r'), default='-', help='sequence in fasta')
    psr_internal.set_defaults(func=main_internal)

    psr_filter = subparsers.add_parser('filter', help=""" filter fasta sequence based on keyword in description""")
    psr_filter.add_argument('-i', type=argparse.FileType('r'), default='-', help='fasta sequence')
    psr_filter.add_argument('-k', help='keyword to filter description')
    psr_filter.set_defaults(func=main_filter)

    psr_cleanhead = subparsers.add_parser('cleanhead', help=""" clean headers of fasta""")
    psr_cleanhead.add_argument('-i', type=argparse.FileType('r'), default='-', help='fasta sequence')
    psr_cleanhead.set_defaults(func=main_cleanhead)

    psr_consensus = subparsers.add_parser('consensus', help="""build consensus sequence from multiple sequence alignment""")
    psr_consensus.add_argument('-i', type=argparse.FileType('r'), default='-', help='multiple alignment in clustal format')
    psr_consensus.add_argument('-p', default='', help='prefix')
    psr_consensus.add_argument('-s', type=float, default=0.8, help='min support [0.8]')
    psr_consensus.add_argument('-l', type=int, default=20, help='min length [20]')
    psr_consensus.set_defaults(func=main_consensus)

    args = parser.parse_args()
    try:
        args.func(args)
    except IOError as e:
        sys.exit()
