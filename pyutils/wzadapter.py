#!/usr/bin/env python
import re
import sys

def opengz(fn):
    
    if fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn)
    else:
        fh = open(fn)

    return fh

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

def reverse(seq):

    return ''.join(reversed(seq))

adaptor = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCACXXXXXXATCTCGTATGCCGTCTTCTGCTTG'

allcnt = [0]*len(adaptor)
C2Tcnt = [0]*len(adaptor)

pattern = re.compile(r'AGAT[CT]G')
seqfh = opengz(sys.argv[1])
adaptor_cnt = 0
for i, seq in enumerate(seqfh):
    if i % 4 == 1:
        for pm in pattern.finditer(seq):
            beg = pm.start()
            seqa = seq[beg:].strip()
            # filter out sequence too short to be identified as adaptor
            if len(seqa) < 10:
                continue
            adaptor_found = True
            for j in xrange(min(len(adaptor), len(seqa))):
                if adaptor[j] == 'X':
                    continue
                if adaptor[j] != seqa[j] and adaptor[j] != 'C':
                    adaptor_found = False
                    break
            if adaptor_found:
                for j in xrange(min(len(adaptor), len(seqa))):
                    if adaptor[j] == 'X':
                        continue
                    if adaptor[j] == 'C':
                        allcnt[j] += 1
                        if seqa[j] == 'T':
                            C2Tcnt[j] += 1
                adaptor_cnt += 1
                break

sys.stdout.write('Analyzed %d adaptors.\n' % adaptor_cnt)
allcnt_merged = 0
C2Tcnt_merged = 0
sys.stdout.write('pos\tall\tT\tfrac\n')
for i, b in enumerate(adaptor):
    if b == 'C':
        sys.stdout.write('%d\t%d\t%d\t%s\n' % (i+1, allcnt[i], C2Tcnt[i], 'NA' if allcnt[i] == 0 else '%1.3f%%' % (float(C2Tcnt[i])/float(allcnt[i])*100,)))
        allcnt_merged += allcnt[i]
        C2Tcnt_merged += C2Tcnt[i]
sys.stdout.write('all C: %d\n' % allcnt_merged)
sys.stdout.write('C2T: %d (%s)\n' % (C2Tcnt_merged, 'NA' if allcnt_merged == 0 else '%1.3f%%' % (float(C2Tcnt_merged)/float(allcnt_merged)*100,)))
