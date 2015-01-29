#!/usr/bin/env python

import argparse
import pysam
import faidx


def count_read_retention(x, refseq, min_qual=0):

   nG2A = nC2T = nG2G = nC2C = 0
   rpos = x.pos
   qpos = 0
   for op, clen in x.cigar:
      if op == 0:            # match
         for p in xrange(clen):
            rc = refseq[rpos+p]
            qc = x.seq[qpos+p]
            # print ord(x.qual[qpos+p])
            if ord(x.qual[qpos+p]) <= min_qual:
               continue
            if rc == 'G':
               if qc == 'G':
                  nG2G += 1
               elif qc == 'A':
                  nG2A += 1
            if rc == 'C':
               if qc == 'C':
                  nC2C += 1
               elif qc == 'T':
                  nC2T += 1
         rpos += clen
         qpos += clen
      elif op == 1:
         qpos += clen
      elif op == 2:
         rpos += clen
      elif op == 4:
         qpos += clen
      else:
         raise Exception("unknown cigar: %d" % op)

   return (nG2A, nC2T, nG2G, nC2C)

def main_icc(args):

   """ incomplete conversion """

   ref = faidx.RefGenome(args.ref)
   refseq = ref.fetch_chrmseq('chrM')

   bam = pysam.Samfile(args.bam)
   n_inc = 0
   n_com = 0
   for x in bam.fetch(reference='chrM'):

      # only consider primary alignment
      if x.is_secondary:
         continue

      strand_tag = dict(x.tags)['ZS']

      (nG2A, nC2T, nG2G, nC2C) = count_read_retention(x, refseq, args.m)
      if strand_tag == '++' or strand_tag == '+-':
         n_inc += nC2C
         n_com += nC2T
      else:                     # -+ or --
         n_inc += nG2G
         n_com += nG2A

      if args.v:
         if strand_tag in ['++', '+-'] and nC2C > 0:
            print '\t'.join(map(str, [nC2C, (nC2C+nC2T), float(nC2C)/(nC2C+nC2T), strand_tag, nC2C, nC2T, nG2G, nG2A, x.qname, x.pos, x.tid, x.flag]))
         elif strand_tag in ['-+', '--'] and nG2G > 0:
            print '\t'.join(map(str, [nG2G, (nG2G+nG2A), float(nG2G)/(nG2G+nG2A), strand_tag, nC2C, nC2T, nG2G, nG2A, x.qname, x.pos, x.tid, x.flag]))

   if not args.v:
      print '\t'.join(map(str, [n_inc, n_com, float(n_inc) / (n_inc+n_com)]))

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description="bisulfite sequence utilities")

   subparsers = parser.add_subparsers()
   parser_icc = subparsers.add_parser("icc", help=""" estimate incomplete conversion estimation from mitchocondrial genome """)
   parser_icc.add_argument("-bam", help="bam file name")
   parser_icc.add_argument("-ref", help="reference fasta")
   parser_icc.add_argument("-v", action="store_true",
                           help="verbose mode, output all the incomplete converted reads, otherwise, output the incomplete conversion ratio")
   parser_icc.add_argument("-m", type=int, default=40, help="minimum base quality (default: 40 before 33 deduction)")
   parser_icc.set_defaults(func=main_icc)

   args = parser.parse_args()
   args.func(args)
