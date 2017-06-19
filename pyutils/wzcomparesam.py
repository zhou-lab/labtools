#!/usr/bin/env python
import sys
# compare two sam file output by aligner without sorting
samfn1 = sys.argv[1]
samfn2 = sys.argv[2]
target = None if len(sys.argv) <= 3 else sys.argv[3] # target to print

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
    self.a.mapq = fields[4]
    self.a.cigar = fields[5]
    self.a.mchrm = fields[6]
    self.a.mpos = fields[7]
    self.a.tlen = fields[8]
    self.a.seq = fields[9]
    self.a.qual = fields[10]
    self.a.aux = ' '.join(fields[11:])
    self.a.line = line

    return True

  def read_read1_alns(self):

    alns = []
    if self.a is not None:
      alns.append(self.a)
      qname = self.a.qname
      self.a = None
    else:
      qname = ''

    while self.read_aln1():
      if qname != '' and qname != self.a.qname:
        break
      alns.append(self.a)
      qname = self.a.qname

    return alns

def main():

  print samfn1, samfn2
  sam1 = SAMFile(samfn1)
  sam2 = SAMFile(samfn2)
  n = 0
  n_allmatch = 0
  n_diff = {'cigar':0, 'pos':0, 'mapq':0, 'seq':0, 'qual':0, 'other':0, 'nhits':0} 

  n_diff_nmulti1 = 0
  n_diff_nmulti2 = 0

  while True:
    as1 = sam1.read_read1_alns()
    as2 = sam2.read_read1_alns()
    if len(as1) == 0 and len(as2) == 0:
      print 'All done. '
      print n, 'processed.'
      print n_allmatch, '(%1.2f %%) all matched.' % (float(n_allmatch) / n * 100)
      for k in ['pos','cigar','mapq','seq','qual','other']:
        print n_diff[k], '(%1.2f %%) single hits differs by %s' % (float(n_diff[k]) / n *100, k)
      print
      print n_diff_nmulti1, ' diff hits from 1 is multi'
      print n_diff_nmulti2, ' diff hits from 2 is multi'
      print n_diff['nhits'], ' diff hits differs by number of nhits'
      break

    if len(as1) == 0 and len(as2) != 0:
      print as1[0].qname, "is not aligned in 1"
      break
  
    if len(as2) == 0 and len(as1) != 0:
      print as2[0].qname, "is not aligned in 2"
      break

    if as1[0].qname != as2[0].qname:
      print as1[0].qname, " in 1 is different from ", as2[0].qname, "in 2"
      print sam1.a.qname, ' is next in 1'
      print sam2.a.qname, ' is next in 2'
      break
  
    n += 1
    if ''.join([_.line for _ in as1]) == ''.join([_.line for _ in as2]):
      n_allmatch += 1
      continue
  
    # single hits
    if len(as1) == len(as2) == 1:
      a1 = as1[0]
      a2 = as2[0]
  
      majorsame = True
      if a1.cigar != a2.cigar:
        n_diff['cigar'] += 1
        majorsame = False
        if target == 'cigar':
          print '%s%s\n' % (a1,a2)
  
      if a1.chrm != a2.chrm or a1.pos != a2.pos:
        n_diff['pos'] += 1
        majorsame = False
        if target == 'pos':
          print '%s%s\n' % (a1,a2)
  
  
      if a1.mapq != a2.mapq:
        n_diff['mapq'] += 1
        majorsame = False
        if target == 'mapq':
          print '%s%s\n' % (a1,a2)
  
  
      if a1.seq != a2.seq:
        n_diff['seq'] += 1
        majorsame = False
        if target == 'seq':
          print '%s%s\n' % (a1,a2)
  
        
      if a1.qual != a2.qual:
        n_diff['qual'] += 1
        majorsame = False
        if target == 'qual':
          print '%s%s\n' % (a1,a2)
  
      if majorsame:
        n_diff['other'] += 1
        if target == 'other':
          print '%s%s\n' % (a1,a2)
  
  
    # multi-hits
    if len(as1) != len(as2):
      n_diff['nhits'] += 1
  
    if len(as1) > 1:
      n_diff_nmulti1 += 1
  
    if len(as2) > 1:
      n_diff_nmulti2 += 1
  
    # print as1[0].qname, as2[0].qname, len(as1), len(as2)
  
  
  
  
if __name__ == "__main__":
  try:
    main()
  except IOError:
    pass
  
  
