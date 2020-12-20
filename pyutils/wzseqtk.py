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

def main_getfasta(args):

    import faidx
    genome = faidx.RefGenome(args.f)
    out = open(args.o, 'w') if args.o is not None else sys.stdout
    for line in args.i:
        fields = line.strip().split('\t')
        chrm = fields[0]
        beg = int(fields[1])
        end = int(fields[2])

        out.write('%s\t%s\n' % (line.strip(), genome.fetch_sequence(chrm, beg+1, end)))

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
    print('n:%d' % (end-beg))
    print('C:%d' % seq.count('C'))
    print('G:%d' % seq.count('G'))
    print('CG:%d' % seq.count('CG'))

def main_runningcomp(args):

    """ export local composition at every base """

    def compute_comp(seq):
        return (seq.count('C'), seq.count('G'), len(seq)-seq.count('N'), seq.count('CG'))
    
    import faidx
    genome = faidx.RefGenome(args.i)
    out = open(args.o, 'w') if args.o is not None else sys.stdout
    for chrm in genome.faidx:
        if args.v:
            err_print(chrm)
        gseq = genome.fetch_chrmseq(chrm).upper()
        if len(gseq) <= args.k*2:
            continue
        c,g,n,cg = compute_comp(gseq[:args.k*2])
        for i in xrange(len(gseq)-args.k*2):
            b1 = gseq[i]
            b2 = gseq[i+args.k*2]
            if b1 == 'N':
                pass
            elif b1 == 'C':
                n -= 1
                c -= 1
                if gseq[i+1] == 'G':
                    cg -= 1
            elif b1 == 'G':
                n -= 1
                g -= 1
            elif b1 == 'A' or b1 == 'T':
                n -= 1
            else:
                raise Exception("Unknown base: %s" % b1)
            
            if b2 == 'N':
                pass
            elif b2 == 'C':
                n += 1
                c += 1
            elif b2 == 'G':
                n += 1
                g += 1
                if gseq[i+args.k*2-1] == 'C':
                    cg += 1
            elif b2 == 'A' or b2 == 'T':
                n += 1
            else:
                raise Exception('Unknown base: %s' % b2)

            if i%args.s == 0:
                out.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (chrm, i+args.k-1, i+args.k, c, g, cg, n))

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
                    tprint([c, prev, prev+2, '+'], out)
                if prev:
                    prev_is_good_left = i-prev >= args.l
                prev = i
        if prev_is_good_left and prev is not None:
            tprint([c, prev, prev+2, '+'], out)

def main_internal(args):

    """ get all the internal CpG with flanking sequence """
    from Bio import SeqIO
    for seq_record in SeqIO.parse(args.i, "fasta"):
        
        seq = seq_record.seq.upper()
        # for i in xrange(args.l, len(seq)-args.l-2):
        for i in xrange(len(seq)-2):
            if args.strict:
                if (len(seq[max(0,i-args.l):i]) != args.l or len(seq[i+2:i+args.l+2]) != args.l):
                    continue
            elif len(seq[max(0,i-args.l):i]) != args.l and len(seq[i+2:i+args.l+2]) != args.l:
                continue
            if seq[i:i+2] == 'CG':
                if args.nobrackets:
                    print('\t'.join(map(
                        str, [".",i,i+2, '|'.join(seq_record.description.split()),
                              seq[max(0,i-args.l):i+args.l+2]])))
                else:
                    print('\t'.join(map(
                        str, [".",i,i+2, '|'.join(seq_record.description.split()),
                              seq[max(0,i-args.l):i]+'['+seq[i:i+2]+']'+seq[i+2:i+args.l+2]])))

def main_filter(args):

    """ filter fasta by keyword """
    from Bio import SeqIO
    for seq_record in SeqIO.parse(args.i, "fasta"):

        if re.search(args.k, seq_record.description):
            print(seq_record.format("fasta"))

def main_cleanhead(args):

    """ replace white space in head of fasta by "_" """

    from Bio import SeqIO
    for seq_record in SeqIO.parse(args.i, "fasta"):

        seq_record.description = '|'.join(seq_record.description.split())
        print(seq_record.format("fasta"))

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
                # [the index of concensus sequence]_[start of concensus sequence column]
                print('>%sconsensus%d_%d' % (args.p, k,j))
                print(wrap(''.join(seq)))
                k += 1
            seq = []


def main_ensembl2name(args):

    gtffh=opengz(args.g)
    id2name = {}
    for line in gtffh:

        fields = line.strip('\n').split('\t')
        if len(fields) < 2:
            continue
        ontol = fields[2]
        if args.gene and ontol == "gene":
            m = re.match(r'.*gene_id "([^"]*)"', line)
            n = re.match(r'.*gene_name "([^"]*)"', line)
            tt = re.match(r'.*gene_biotype "([^"]*)"', line)
            id2name[m.group(1)] = (n.group(1), tt.group(1) if tt else '.')
            m = None
            n = None
            tt = None
        if args.transcript and ontol == "transcript":
            m = re.match(r'.*transcript_id "([^"]*)"', line)
            n = re.match(r'.*gene_name "([^"]*)"', line)
            tt = re.match(r'.*transcript_biotype "([^"]*)"', line)
            id2name[m.group(1)] = (n.group(1), tt.group(1) if tt else '.')
            m = None
            n = None
            tt = None

    sys.stderr.write("Found %d id-name associations\n" % len(id2name))

    for ii, line in enumerate(args.i):

        fields = line.strip('\n').split('\t')
        if args.header and ii==0:
            args.o.write('\t'.join(fields)+'\tGeneName\tGeneType\n')
            continue

        if len(fields) <= args.c:
            sys.stderr.write(line+"not recognized\n")
            sys.exit()
            continue
        k = fields[args.c-1]

        gg = []
        tt = []
        for kk in k.split('+'):
            if kk in id2name:
                gg.append(id2name[kk][0])
                tt.append(id2name[kk][1])
            else:
                gg.append('.')
                tt.append('.')
        
        args.o.write("%s\t%s\t%s\n" % (
            line.strip('\n'), '+'.join(gg), '+'.join(tt)))


def main_cnt2rpkm(args):

    totalreads = []
    for i, line in enumerate(args.i):
        if i==0:
            continue
        elif i==1:
            fields = line.strip().split('\t')
            bams = fields[6:]
            args.o.write('ID\tchrm\tbeg\tend\tlength')
            for bam in bams:
                
                with open(bam+'.flagstat') as flagstat:
                    for _line in flagstat:
                        m = re.match(r'(\d+) \+ \d+ mapped \(', _line)
                        if m:
                            totalreads.append(int(m.group(1)))
                args.o.write('\t'+bam)
            args.o.write('\n')
        else:
            fields = line.strip().split('\t')
            tlen = int(fields[5])
            beg_coords = map(int, fields[2].split(';'))
            end_coords = map(int, fields[3].split(';'))
            args.o.write('%s\t%s\t%d\t%d\t%d\t%s\n' % (
                fields[0], fields[1].split(';')[0],  min(beg_coords), max(end_coords), tlen,
                '\t'.join(['%1.3f' % (float(cnt)/tlen*1000/totalreads[j]*1000000) for j,cnt in enumerate(fields[6:])])))

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='calculate sequence')
    subparsers = parser.add_subparsers()

    psr_printc = subparsers.add_parser("printc", help=""" print all C in fasta """)
    psr_printc.add_argument('-i', required=True, help="sequence format in fasta")
    psr_printc.add_argument('-o', default=None, help='output file (default stdout)')
    psr_printc.add_argument('-v', action='store_true', help='print status')
    psr_printc.set_defaults(func=main_printc)

    parser_getfasta = subparsers.add_parser('getfasta', help=' get sequence from indexed fasta ')
    parser_getfasta.add_argument('-f', required=True, help='reference in fasta format')
    parser_getfasta.add_argument('-i', type=argparse.FileType('r'), default='-', help='input bed')
    parser_getfasta.add_argument('-o', help='output', default=None)
    parser_getfasta.set_defaults(func=main_getfasta)

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
    psr_internal.add_argument('--nobrackets', action='store_true', help='no output of brackets around CpG')
    psr_internal.add_argument('--strict', action='store_true', help='require length on both flanking regions')
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

    parser_runningcomp = subparsers.add_parser('runningcomp', help='compute running composition')
    parser_runningcomp.add_argument('-i', required=True, help='genomic sequence format in fasta')
    parser_runningcomp.add_argument('-k', type=int, default=1000, help='flanking length')
    parser_runningcomp.add_argument('-s', type=int, default=100, help='step size')
    parser_runningcomp.add_argument('-v', action='store_true', help='print status')
    parser_runningcomp.add_argument('-o', help='output', default=None)
    parser_runningcomp.set_defaults(func=main_runningcomp)

    parser_ensembl2name = subparsers.add_parser('ensembl2name', help='convert ENSEMBL id to gene name')
    parser_ensembl2name.add_argument('-i', type=argparse.FileType('r'), default='-', help='input table')
    parser_ensembl2name.add_argument('-g', type=str, default=None, help='GTF file')
    parser_ensembl2name.add_argument('-G', '--gene', action="store_true", help='map gene id')
    parser_ensembl2name.add_argument('-T', '--transcript', action="store_true", help='map transcript id')
    parser_ensembl2name.add_argument('-H', '--header', action='store_true', help='treat first row as header')
    parser_ensembl2name.add_argument('-c', type=int, default=1, help='column id for ENSEMBL id in the input table (1-based)')
    parser_ensembl2name.add_argument('-o', type=argparse.FileType('w'), help='output', default=sys.stdout)
    parser_ensembl2name.set_defaults(func=main_ensembl2name)

    parser_cnt2rpkm = subparsers.add_parser('cnt2rpkm', help='count table from feature count to rpkm')
    parser_cnt2rpkm.add_argument('-i', type=argparse.FileType('r'), default='-', help='count table')
    parser_cnt2rpkm.add_argument('-o', type=argparse.FileType('w'), help='output', default=sys.stdout)
    parser_cnt2rpkm.set_defaults(func=main_cnt2rpkm)

    args = parser.parse_args()
    try:
        args.func(args)
    except IOError as e:
        sys.exit()
