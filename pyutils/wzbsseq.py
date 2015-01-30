#!/usr/bin/env python

import argparse
import pysam
import faidx
import multiprocessing
import traceback
import sys

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


def occurrences(string, sub):
    count = start = 0
    while True:
        start = string.find(sub, start) + 1
        if start > 0:
            count+=1
        else:
            return count

class Consumer(multiprocessing.Process):

    def __init__(self, window_queue, result_queue, winlen, vcf_path, reference_path):
        multiprocessing.Process.__init__(self)
        self.window_queue = window_queue
        self.result_queue = result_queue
        self.reference = faidx.RefGenome(reference_path)
        self.vcf = pysam.Tabixfile(vcf_path)
        self.winlen = winlen

    def run(self):

        while True:
            window = self.window_queue.get()
            if window is None:
                self.window_queue.task_done()
                break
            result = self.process_window(window)
            if result:  # this is important!
                self.result_queue.put(result)
                # print 'still putting'
            self.window_queue.task_done()

        return

    def process_window(self, (chrm, winbeg)):

        try:
            return self._process_window(chrm, winbeg)
        except:
            traceback.print_exc()
            raise Exception()

    def _process_window(self, chrm, winbeg):

        seqbeg = winbeg-1
        seqend = winbeg+self.winlen+1
        global mutex1, mutex2
        seq = self.reference.fetch_sequence(chrm, seqbeg, seqend).upper()
        cN = seq.count('N')
        cG = seq.count('G')
        cC = seq.count('C')
        cA = seq.count('A')
        cT = seq.count('T')
        if cN > self.winlen * 0.1:
            return None
        cCG = occurrences(seq, 'CG')
        cCC = occurrences(seq, 'CC')
        cCT = occurrences(seq, 'CT')
        cCA = occurrences(seq, 'CA')
        cAG = occurrences(seq, 'AG')
        cTG = occurrences(seq, 'TG')
        cGG = occurrences(seq, 'GG')
        rCG1 = rCG2 = rCA = rCC = rCT = rAG = rTG = rGG = 0
        for line in self.vcf.fetch(chrm, winbeg, winbeg+self.winlen):
            fields = line.split('\t')
            if len(fields) < 7:
                print line
                raise Exception(line)
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            if alt == '.':
                if ref == 'C':
                    rs1, rs2 = seq[pos-seqbeg:pos-seqbeg+2]
                    # assert(ref == rs1)
                    if rs2 == 'G':
                        rCG1 += 1
                    elif rs2 == 'C':
                        rCC += 1
                    elif rs2 == 'A':
                        rCA += 1
                    elif rs2 == 'T':
                        rCT += 1
                if ref == 'G':
                    if pos == winbeg:
                        continue
                    rs0, rs1 = seq[pos-seqbeg-1:pos-seqbeg+1]
                    # assert(ref == rs1)
                    if rs0 == 'C':
                        rCG2 += 1
                    elif rs0 == 'A':
                        rAG += 1
                    elif rs0 == 'T':
                        rTG += 1
                    elif rs0 == 'G':
                        rGG += 1

        # print [cCG, cCC, cCT, cCA, cAG, cTG, cGG]
        # print [rCG, rCC, rCT, rCA, rAG, rTG, rGG]
        # print seq

        return [chrm, winbeg, winbeg+self.winlen,
                cG, cC, cA, cT, float(cG+cC) / float(cG+cC+cA+cT),
                cCG, cCC, cCT, cCA, cAG, cTG, cGG,
                rCG1, rCG2, rCC, rCT, rCA, rAG, rTG, rGG]


class Writer(multiprocessing.Process):

    def __init__(self, result_queue):
        multiprocessing.Process.__init__(self)
        self.result_queue = result_queue

    def run(self):
        sys.stdout.write('\t'.join([
            'chrm', 'winbeg', 'winend',
            'cG', 'cC', 'cA', 'cT', 'GCcontent',
            'cCG',
            'cCC', 'cCT', 'cCA',
            'cAG', 'cTG', 'cGG',
            'rCG1', 'rCG2',
            'rCC', 'rCT', 'rCA',
            'rAG', 'rTG', 'rGG'])+'\n')
        while True:
            result = self.result_queue.get()
            if result is None:
                break
            sys.stdout.write('\t'.join(map(str, result))+'\n')

def main_CpH(args):

    # processes = multiprocessing.cpu_count() - 2
    num_processes = 25

    if args.chrm:
        chrms = [args.chrm]
    else:
        chrms = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chrM']
        

    windows = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    consumers = []
    for i in xrange(num_processes):
        c = Consumer(windows, results, args.winlen, args.vcf, args.ref)
        c.start()
        consumers.append(c)

    writer = Writer(results)
    writer.start()

    # create windows
    r = faidx.RefGenome(args.ref)
    for chrm in chrms:
        chrmlen = r.faidx[chrm][0]
        step = args.winlen
        # start from the second base, since we are looking at binucleotides
        for winbeg in xrange(2, chrmlen-args.winlen, step):
            windows.put((chrm, winbeg))

    # put poison pills
    for i in xrange(num_processes):
        windows.put(None)

    # print 'before joining'
    windows.join()
    # print 'windows joined'

    # consumers need to be joined so that
    # the items put on the queue got flushed
    for c in consumers:
        c.join()

    results.put(None)
    writer.join()
    # print 'writer joined'

    return


if __name__ == '__main__':
   parser = argparse.ArgumentParser(description="bisulfite sequence utilities")

   subparsers = parser.add_subparsers()
   parser_icc = subparsers.add_parser("icc", help="estimate incomplete conversion estimation from mitchocondrial genome")
   parser_icc.add_argument("-bam", help="bam file name")
   parser_icc.add_argument("-ref", help="reference fasta")
   parser_icc.add_argument("-v", action="store_true",
                           help="verbose mode, output all the incomplete converted reads, otherwise, output the incomplete conversion ratio")
   parser_icc.add_argument("-m", type=int, default=40, help="minimum base quality (default: 40 before 33 deduction)")
   parser_icc.set_defaults(func=main_icc)

   parser_cph = subparsers.add_parser('CpH', help="compute CpH spectrum")
   parser_cph.add_argument("-vcf", help="vcf file name")
   parser_cph.add_argument("-ref", help="reference file name")
   parser_cph.add_argument("-winlen", type=int, default=1000, help='window size')
   parser_cph.add_argument("-chrm", default=None, help='chromosome if not specified, all chromosome')
   parser_cph.set_defaults(func=main_CpH)

   args = parser.parse_args()
   args.func(args)
