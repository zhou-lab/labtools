#!/RIS/HPC_apps/AMD/python/python-2.7.2/bin/python
import os, sys
from collections import deque
import gzip
from subprocess import Popen, PIPE, check_call
import re
from pbsgen import Job
from glob import glob

class VRecord:

    def __init__(self, chrm, pos):
        self.chrm = chrm
        self.pos = int(pos)

    def format_vcf(self, samples):
        return "%s\t%d\t%s\t%s\n" % (self.chrm, self.pos, "\t".join(self.text2to8), "\t".join([self.data[_] for _ in samples]))

class VCF:

    def __init__(self, fn):
        self.fn = fn

    def open(self):
        if self.fn.endswith(".gz"):
            self.fh = gzip.open(self.fn, "r")
        else:
            self.fh = open(self.fn, "r")

    def close(self):
        self.fh.close()

    def read_header(self):

        while (True):
            line = self.fh.readline()
            if not line:
                break
            if line.startswith("#CHROM"):
                self.samplelist = line.split()[9:]
                break

    def read1(self):

        while (True):
            line = self.fh.readline()
            if not line:
                break
            if line[0] == '#':
                continue
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.id = pair[2]
            r.data = dict(zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]]))
            yield r

    def fetch_region(self, chrm, beg, end):

        for line in Popen(["tabix", self.fn, "%s:%d-%d" % (chrm, beg, end)], stdout=PIPE).stdout:
            pair = line.split("\t")
            r = VRecord(pair[0], pair[1])
            r.text2to8 = pair[2:9]
            r.data = dict(zip(self.samplelist, [_.split(":")[0] for _ in pair[9:]]))
            yield r

def count_read(bam_fn, chrm, loc):

    out = Popen(['samtools', "view", "-c", bam_fn, "%s:%d-%d" % (chrm, loc, loc+1)], stdout=PIPE).stdout
    lines = out.readlines()
    if len(lines)!=1:
        stderr.write("[Warning] multiple lines %s\n", str(lines))
        return 0
    return int(lines[0].strip())

def is_overlap(rsnp, rdel):

    if rsnp.chrm == rdel.chrm and abs(rsnp.pos - rdel.pos)<10000:
        return True
    else:
        return False


def cmp_vcf():
    vcf1 = VCF(sys.argv[1])
    vcf2 = VCF(sys.argv[2])

    vcf1.open()
    vcf2.open()

    vcf1.read_header()
    vcf2.read_header()

    intersample = set(vcf2.samplelist) & set(vcf1.samplelist)
    print len(intersample), "overlapping samples"
    
    cs2g2 = {}
    for r in vcf2.read1():
        for s in intersample:
            cs2g2[(r.id, s)] = r.data[s]

    cs2g1 = {}
    for r in vcf1.read1():
        for s in intersample:
            cs2g1[(r.id, s)] = r.data[s]

    print len(cs2g1), "call-sample pairs in vcf1"
    print len(cs2g2), "call-sample pairs in vcf2"

    overlap = set(cs2g1.keys()) & set(cs2g2.keys())
    print len(overlap), "overlapping call-sample pairs between the two vcfs"

    vars = ["./.", "0/0", "0/1", "1/1"]
    print 'vcf1\\vcf2', '\t'.join(vars)
    for var1 in vars:
        sys.stdout.write(var1)
        for var2 in vars:
            sys.stdout.write('\t%d' % len([_ for _ in overlap if cs2g1[_] == var1 and cs2g2[_] == var2]))
        sys.stdout.write('\n')

cmp_vcf()
