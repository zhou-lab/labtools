#!/usr/bin/env python
import pysam
import argparse
import re
import faidx

b2c = {
    'A': '\033[32m',
    'G': '\033[33m',
    'C': '\033[36m',
    'T': '\033[31m'
}

GAPC = '\033[34m'
ENDC = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'

def q2c(qual):
    if qual < 10:
        return '\033[34m'
    elif qual < 20:
        return '\033[33m'
    elif qual < 30:
        return '\033[32m'
    else:
        return '\033[37m'
    
# class bcolors:
#     HEADER = '\033[95m'
#     OKGREEN = 
#     WARNING = 
#     FAIL = 
#     ENDC = 

flank_len = 10

def rprint(seq):
    s = ''
    for base in seq:
        s += b2c[base.upper()]+base+ENDC
    return s

def qprint(read, beg, end, refseq, refbeg, mod=''):
    s = ''
    for i, base in enumerate(read.seq[beg:end]):
        c = b2c[base.upper()]
        if read.is_reverse:
            if base == refseq[refbeg+i]:
                base = c+mod+','+ENDC
            else:
                base = c+mod+base.lower()+ENDC
        else:
            if base == refseq[refbeg+i]:
                base = c+mod+'.'+ENDC
            else:
                base = c+mod+base.upper()+ENDC
        s += base

    return s

def qualprint(read, beg, end):
    s = ''
    for i, base in enumerate(read.seq[beg:end]):
        c = q2c(ord(read.qual[beg+i])-33)
        if read.is_reverse:
            base = c+base.lower()+ENDC
        else:
            base = c+base.upper()+ENDC
        s += base

    return s

def main(args):

    samfile=pysam.Samfile(args.bam)
    ref = faidx.RefGenome(args.ref)

    read1seen = read2seen = False
    for read in samfile.fetch(args.reg):
	if read.qname != args.qname:
            continue

        chrm = samfile.getrname(read.tid)
        refbeg = read.pos - 100
        refend = read.pos + len(read.seq) + 100
        refseq = ref.fetch_sequence(chrm, refbeg, refend).upper()
        rpos = read.pos - refbeg + 1
        qpos = 0

        op,oplen = read.cigar[0]
        if op == 4:
            pr_beg = rpos - flank_len - oplen
        else:
            pr_beg = rpos - flank_len
        pp = pr = pq = pqual = ''
        pp += " "*(rpos-pr_beg)
        pr += rprint(refseq[pr_beg:rpos])
        pq += ' '*flank_len
        pqual += ' '*flank_len

	for i, (op, clen) in enumerate(read.cigar):
            if op == 0:
                if pp.isspace():
                    pp += "|{}:{}".format(chrm, rpos+refbeg-1)
                pr += rprint(refseq[rpos:rpos+clen])
                pq += qprint(read, qpos, qpos+clen, refseq, rpos)
                pqual += qualprint(read, qpos, qpos+clen)
                rpos += clen
                qpos += clen
            elif op == 1:
                pr += GAPC+'*'*clen+ENDC
                pq += qprint(read, qpos, qpos+clen, refseq, rpos)
                pqual += qualprint(read, qpos, qpos+clen)
                qpos += clen
            elif op == 2:
                pr += rprint(refseq[rpos:rpos+clen])
                pq += GAPC+'*'*clen+ENDC
                pqual += '*'*clen
                rpos += clen
            elif op == 4:
                pq += qprint(read, qpos, qpos+clen, refseq, rpos, UNDERLINE)
                pqual += qualprint(read, qpos, qpos+clen)
                qpos += clen
            else:
                raise Exception("unknown cigar: %d" % op)
        pr += rprint(refseq[rpos:rpos+flank_len])

        print "\n"+"="*5+read.qname+"="*5
        print pp
        print pr
        print pq
        print pqual
        if args.pread:
            print read

        if read.is_read1: read1seen = True
        if read.is_read2: read2seen = True
        if args.pair == '12' and read1seen and read2seen:
            break
        elif args.pair == '1' and read1seen:
            break
        elif args.pair == '2' and read2seen:
            break
        elif args.pair == '0':
            break


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-qname')
    parser.add_argument('-reg')
    parser.add_argument('-ref')
    parser.add_argument('-bam')
    parser.add_argument('-pqual', action='store_true')
    parser.add_argument('-pread', action='store_true')
    parser.add_argument('-pair', default='12')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)

