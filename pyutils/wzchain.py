#!/usr/bin/env python
import argparse
import sys, datetime
from wzcore import *

MCHC = '\033[34m'
MISC = '\033[33m'
ENDC = '\033[0m'

def psort(i,j):
    if i<j:
        return i,j
    else:
        return j,i
    
def process_chain(hchr, hstr, hpos, mchr, mstr, mpos, m, s, t):
    
    hp = [] # human chromosome coordinate on each position
    mp = [] # mouse chromosome coordinate on each position
    # ip = [] # index in chns on each position
    gaps = []
    # chns = []
    hp1 = hpos
    mp1 = mpos
    ngap = 0
    for i, m1 in enumerate(m):
        
        # match
        for j in xrange(m1):
            hp.append((hp1+j) if hstr == '+' else (hp1-j))
            mp.append((mp1+j) if mstr == '+' else (mp1-j))
            # ip.append((i,i))
            gaps.append(ngap)

        hp1 += m1 if hstr == '+' else -m1
        mp1 += m1 if mstr == '+' else -m1
        # hp2 = hp1 + m1 if hstr == '+' else -m1
        # mp2 = mp1 + m1 if mstr == '+' else -m1
        # chns.append((min(hp1,hp2),max(hp1,hp2), min(mp1,mp2), max(mp1,mp2)))
        # hp1 = hp2
        # mp1 = mp2

        # gap
        if i != len(m)-1:
            s1 = s[i]
            t1 = t[i]
            mst = max(s1,t1)
            hp2 = (hp1+s1) if hstr == '+' else (hp1-s1)
            mp2 = (mp1+t1) if mstr == '+' else (mp1-t1)
            for j in xrange(mst):
                hp.append(hp1 if j<mst/2 else hp2)
                mp.append(mp1 if j<mst/2 else mp2)
                gaps.append(ngap)
                # ip.append((i,i+1))
                ngap += 1
            hp1 = hp2
            mp1 = mp2
    
    sys.stdout.flush()
    stepsize = 300
    windowsize = 3000
    maxgap = 0.5
    for i in xrange(0,len(hp)-windowsize,stepsize):
        
        ngap = gaps[i+windowsize] - gaps[i]
        if ngap > windowsize*maxgap:
            continue
        h1,h2 = psort(hp[i], hp[i+windowsize])
        m1,m2 = psort(mp[i], mp[i+windowsize])
        # i1, i2 = ip[i][1], ip[i+windowsize][0]
        # for j in xrange(i1, i2+1):
        #     hp1, hp2, mp1, mp2 = chns[j]
        #     assert (max(hp1,h1) < min(hp2,h2)) and (max(mp1,m1) < min(mp2,m2))
        #     yield hchr, max(hp1, h1), min(hp2, h2), hstr, h1, h2, hpos, mchr, max(mp1, m1), min(mp2, m2), mstr, m1, m2, mpos
        yield hchr, h1, h2, mchr, m1, m2, '%d/%s/%d/%s' % (hpos, hstr, mpos, mstr)

def visualize_chain(hchr, hstr, hpos, hend, mchr, mstr, mpos, mend, m, s, t, m_genome, h_genome, fhvis):

    href = h_genome.fetch_sequence(hchr,min(hpos,hend), max(hpos,hend))
    mref = m_genome.fetch_sequence(mchr,min(mpos,mend), max(mpos,mend))
    if hstr == '-':
        href = reverse_complement(href)
    if mstr == '-':
        mref = reverse_complement(mref)
    haln = []
    maln = []
    hi = 0
    mi = 0
    for i, m1 in enumerate(m):
        haln.extend(href[hi:hi+m1])
        maln.extend(mref[mi:mi+m1])
        hi += m1
        mi += m1
        if i != len(m)-1:
            # mst = max(s[i], t[i])
            # haln.extend(['-']*mst)
            # maln.extend(['-']*mst)

            # save some space
            haln.append('-')
            maln.append('-')
            hi += s[i]
            mi += t[i]

    i = 0
    step = 80
    fhvis.write('%s\n' % '\t'.join(map(str, [hchr, hstr, hpos, mchr, mstr, mpos])))
    while True:
        if i!=0: fhvis.write('+\n')

        hs = []
        ms = []
        for j in xrange(step):
            if i+j < len(haln):
                if haln[i+j] == maln[i+j]:
                    hs.append(MCHC+haln[i+j]+ENDC)
                    ms.append(MCHC+maln[i+j]+ENDC)
                else:
                    hs.append(MISC+haln[i+j]+ENDC)
                    ms.append(MISC+maln[i+j]+ENDC)
        fhvis.write('%s\n' % ''.join(hs))
        fhvis.write('%s\n' % ''.join(ms))

        # fhvis.write('%s\n' % ''.join(haln[i:i+step]))
        # fhvis.write('%s\n' % ''.join(maln[i:i+step]))
        i += step
        if i > len(haln):
            break

def slice_chains(chain_fn, proc_func, data):

    cnt = 0
    for line in opengz(chain_fn):
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split()
        if line.startswith('chain'):
            hchr = fields[2]
            hchrlen = int(fields[3])
            hstr = fields[4]
            if hstr == '+':
                hpos, hend = int(fields[5])+1, int(fields[6]) # 0-based begin
            else:
                hpos, hend = hchrlen-int(fields[6])-1, hchrlen-int(fields[5])
            mchr = fields[7]
            mchrlen = int(fields[8])
            mstr = fields[9]
            if mstr == '+':
                mpos, mend = int(fields[10])+1, int(fields[11]) # 0-based begin
            else:
                mpos, mend = mchrlen-int(fields[11])-1, mchrlen-int(fields[10])
            m = [] # matches
            s = [] # human skips
            t = [] # mouse skips
        elif line == '\n':
            continue
        elif len(fields) == 3:
            m.append(int(fields[0]))
            s.append(int(fields[1]))
            t.append(int(fields[2]))
        elif len(fields) == 1:

            m.append(int(fields[0]))
            cnt += 1
            if cnt <=300:
                err_print('%s Processing chain of size %d/%d' %
                          (datetime.datetime.now(), abs(hend-hpos), abs(mend-mpos)))
                sys.stderr.flush()

            proc_func(hchr, hstr, hpos, hend, mchr, mstr, mpos, mend, m, s, t, data)
            # visualize_chain(hchr, hstr, hpos, hend, mchr, mstr, mpos, mend,
            #                 m, s, t, m_genome, h_genome, fhvis)
            # for chain in process_chain(hchr, hstr, hpos, mchr, mstr, mpos, m, s, t):
            #     fh.write('\t'.join(map(str, chain))+'\n')
            #     if chain[5] < 0:
            #         fh.flush()
            #     assert(chain[5] >= 0)

    fhvis.close()
    fh.close()

def main_toBed(args):


    def proc_func_toBed(hchr, hstr, hpos, hend, mchr, mstr, mpos, mend, m, s, t, data):

        htoken = '%s:%d_%d:%s' % (hchr, hpos, hend, hstr)
        mtoken = '%s:%d_%d:%s' % (mchr, mpos, mend, mstr)
        hp = []
        mp = []
        gaps = []
        hp1 = 0
        mp1 = 0

        hgenome, mgenome = data
        if hgenome:
            hseq = hgenome.fetch_sequence(hchr, hpos, hend, uppercase=True)
            if hstr == '-':
                hseq = reverse_complement(hseq)
        if mgenome:
            mseq = mgenome.fetch_sequence(mchr, mpos, mend, uppercase=True)
            if mstr == '-':
                mseq = reverse_complement(mseq)

        for i, m1 in enumerate(m):

            wprint('\t'.join(map(str, [
                hchr,
                hpos + hp1 if hstr == '+' else hend - hp1 - m1,
                hpos + hp1 + m1 if hstr == '+' else hend - hp1,
                htoken,
                mchr,
                mpos + mp1 if mstr == '+' else hend - mp1 - m1,
                mpos + mp1 + m1 if mstr == '+' else mend - mp1,
                mtoken])))

            if hseq and mseq:
                hs = []
                ms = []
                for j in xrange(m1):
                    hb = hseq[hp1+j]
                    mb = mseq[mp1+j]
                    if hb == mb:
                        hs.append(MCHC+hb+ENDC)
                        ms.append(MCHC+mb+ENDC)
                    else:
                        hs.append(MISC+hb+ENDC)
                        ms.append(MISC+mb+ENDC)
                wprint(''.join(hs))
                wprint(''.join(ms))
            hp1 += m1
            mp1 += m1

            if i != len(m) - 1:
                hp1 += s[i]
                mp1 += t[i]

    import faidx
    h_genome = faidx.RefGenome(args.g1) if args.g1 else None
    m_genome = faidx.RefGenome(args.g2) if args.g2 else None
    slice_chains(args.i, proc_func_toBed, (h_genome, m_genome))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='process chain file')
    subparsers = parser.add_subparsers()

    parser_toBed = subparsers.add_parser('toBed', help='convert a chain file to bed')
    parser_toBed.add_argument('-i', required=True, help='chain file')
    parser_toBed.add_argument('-g1', default=None, help='fasta of genome 1')
    parser_toBed.add_argument('-g2', default=None, help='fasta of genome 2')
    parser_toBed.set_defaults(func=main_toBed)

    args = parser.parse_args()

    try:
        args.func(args)
    except IOError as e:
        sys.exit()

    
