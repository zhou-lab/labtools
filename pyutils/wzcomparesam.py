#!/usr/bin/env python
import sys
import argparse
# compare two sam file output by aligner without sorting

class Alignment:

  def __str__(self):

    return self.line

class SAMFile:

  def __init__(self, fname):

    self.fh = open(fname, 'r')
    self.a = None

  def read_aln1(self):

    line = self.fh.readline()

    if line == '':
      return False
    fields = line.strip('\r\n').split('\t')
    self.a = Alignment()
    self.a.qname = fields[0]
    self.a.flag = fields[1]
    if int(self.a.flag) & 0x40:
      self.a.qname += '_1'
    elif int(self.a.flag) & 0x80:
      self.a.qname += '_2'
    self.a.chrm = fields[2]
    self.a.pos = fields[3]
    self.a.mapq = int(fields[4])
    self.a.cigar = fields[5]

    self.a.line = line

    if len(fields) <= 6:
      self.a.mchrm = '.'
      self.a.mpos = -1
      self.a.tlen = -1
      self.a.seq = '.'
      self.a.qual = '.' 
      self.a.aux = '.' 
    else:
      self.a.mchrm = fields[6]
      self.a.mpos = fields[7]
      self.a.tlen = fields[8]
      self.a.seq = fields[9]
      self.a.qual = fields[10]
      self.a.aux = ' '.join(fields[11:])

    return True

  def read_read1_alns(self):
    # read all alignment corresponding to a read name

    alns = []
    if self.a is not None:
      alns.append(self.a)
      qname, readInPair = (self.a.qname, int(self.a.flag) & 0x40)
      self.a = None
    else:
      qname, readInPair = ('', True)

    while self.read_aln1():
      if qname != '' and (qname != self.a.qname or readInPair != int(self.a.flag) & 0x40):
        break
      alns.append(self.a)
      qname = self.a.qname
      readInPair = int(self.a.flag) & 0x40

    return alns

  def read_readall_alns(self):

    qname2alns = {}
    while True:
      as1 = self.read_read1_alns()
      if len(as1) == 0:
        break
      qname2alns[(as1[0].qname, int(as1[0].flag)&0x40)] = as1
      
    return qname2alns

def compare_sam_reads(as1, as2, n):

    if ''.join([_.line for _ in as1]) == ''.join([_.line for _ in as2]):
      n['all_match'] += 1
      return
  
    # single hits
    if len(as1) == len(as2) == 1:
      a1 = as1[0]
      a2 = as2[0]
  
      majorsame = True
      if a1.cigar != a2.cigar:
        n['diff_cigar'] += 1
        majorsame = False
        if args.t == 'cigar':
          print('file1: %sfile2: %s\n' % (a1,a2))

      if a1.chrm != a2.chrm or a1.pos != a2.pos:
        n['diff_pos'] += 1
        majorsame = False
        if args.t == 'pos' and ((args.m is None) or a1.mapq >= args.m or a2.mapq >= args.m):
          print('file1: %sfile2: %s\n' % (a1,a2))
  
  
      if a1.mapq != a2.mapq:
        n['diff_mapq'] += 1
        majorsame = False
        if args.t == 'mapq':
          print('file1: %sfile2: %s\n' % (a1,a2))
  
  
      if a1.seq != a2.seq:
        n['diff_seq'] += 1
        majorsame = False
        if args.t == 'seq':
          print('file1: %sfile2: %s\n' % (a1,a2))
  
        
      if a1.qual != a2.qual:
        n['diff_qual'] += 1
        majorsame = False
        if args.t == 'qual':
          print('file1: %sfile2: %s\n' % (a1,a2))
  
      if majorsame:
        n['diff_other'] += 1
        if args.t == 'other':
          print('file1: %sfile2: %s\n' % (a1,a2))
  
    # multi-hits
    if len(as1) != len(as2):
      n['diff_nhits'] += 1
  
    if len(as1) > 1:
      n['nmulti1'] += 1
  
    if len(as2) > 1:
      n['nmulti2'] += 1
  

def main(args):

  sam1 = SAMFile(args.samfn1)
  sam2 = SAMFile(args.samfn2)
  n = {
    'all': 0,                   # number of reads processed
    'all_match':0,              # number of all matched
    '1_only': 0,                # alignment only in sam file 1
    '2_only': 0,                # alignment only in sam file 2
    'diff_cigar':0,
    'diff_pos':0,
    'diff_mapq':0,
    'diff_seq':0,
    'diff_qual':0,
    'diff_other':0,
    'diff_nhits':0,             # number of hits
    'nmulti1':0,           # number of multi-hit alignment in sam file 1
    'nmulti2':0            # number of multi-hit alignment in sam file 2
  }


  if args.m:                    # in memory processing

    qname2alns1 = sam1.read_readall_alns()
    qname2alns2 = sam2.read_readall_alns()
    q1 = set(qname2alns1.keys())
    q2 = set(qname2alns2.keys())
    n['1_only'] = len(q1-q2)
    n['2_only'] = len(q2-q1)
    n['all'] = len(q1 & q2)
    for k in q1&q2:
      compare_sam_reads(qname2alns1[k], qname2alns2[k], n)

  else:                         # on-line processing, require names to be sorted

    while True:
      as1 = sam1.read_read1_alns()
      as2 = sam2.read_read1_alns()

      if len(as1) == 0 and len(as2) != 0:
        n['2_only'] += 1
        continue

      if len(as2) == 0 and len(as1) != 0:
        n['1_only'] += 1
        continue

      if len(as1) == 0 and len(as2) == 0:
        break
      
      if as1[0].qname != as2[0].qname:
        print(as1[0].qname, " in 1 is different from ", as2[0].qname, "in 2")
        print(sam1.a.qname, ' is next in 1')
        print(sam2.a.qname, ' is next in 2')
        break

      compare_sam_reads(as1, as2, n)
      n['all'] += 1
    
  # report
  if n['all'] == 0:
    return
  
  print('All done. ')
  print('%d processed.' % n['all'])
  print('%d reads only aligned in file 1' % n['1_only'])
  print('%d reads only aligned in file 2' % n['2_only'])
  print('%d (%1.2f %%) all matched.' % (n['all_match'], float(n['all_match']) / n['all'] * 100))
  for k in ['pos','cigar','mapq','seq','qual','other']:
    print('%d (%1.2f %%) single hits differs by %s' % (n['diff_'+k], float(n['diff_'+k]) / n['all'] *100, k))
  print('')
  print('%d diff hits from 1 is multi' % n['nmulti1'])
  print('%d diff hits from 2 is multi' % n['nmulti2'])
  print('%d diff hits differs by number of nhits' % n['diff_nhits'])

if __name__ == "__main__":
  try:
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('samfn1', help='sam file 1 to compare')
    parser.add_argument('samfn2', help='sam file 2 to compare')
    parser.add_argument('-t', default=None, help='target to print')
    parser.add_argument('-q', default=None, help='minimum mapq to print')
    parser.add_argument('-m', action='store_true', help='read everything into memory, this avoids sorting')
    args = parser.parse_args()
    main(args)
  except IOError:
    pass
  
  
