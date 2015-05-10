import pandas as pd


def normalize_chrm(chrm):

    if chrm == '23' or chrm == 'chr23': chrm = 'X'
    if chrm == '24' or chrm == 'chr24': chrm = 'Y'
    if chrm == '25' or chrm == 'chr25': chrm = 'M'
    if chrm == 'MT' or chrm == 'chrMT': chrm = 'M'
    if chrm.isdigit() or chrm in ['X', 'Y', 'M']:
        # not chrm.startswith('chr')
        chrm = 'chr'+chrm

    return chrm

def _overlap(b1, e1, b2, e2):

    if b1 <= e2 and b2 <= e1:
        return min(e1,e2) - max(b1,b2)
    else:
        return 0

def overlap(r1, r2):

    if r1.chrm != r2.chrm:
        return 0
    return _overlap(r1.beg, r1.end, r2.beg, r2.end)


def load_bed(bed, sort=True):

    b = Bed()
    b.data = pd.read_table(bed, header=None, sep='\t')
    b.data.rename(columns={0:'chrm',1:'beg',2:'end'},inplace=True)
    if sort:
        b.data.sort(columns=['chrm','beg','end'], inplace=True)

    return b

def nextrec(it):

    res = next(it, None)
    if res is None:
        return None
    else:
        return res[1]

class Bed():

    def __init__(self, binsize=10000, issorted=False):
        self.data = None
        self.bin2record = {}
        self.binsize = binsize
        self.issorted = issorted

    def _bin_add_record(self, k, r):
        if k in self.bin2record:
            self.bin2record[k].append(r)
        else:
            self.bin2record[k] = [r]
        
    def index(self):
        for index, r in self.data.iterrows():
            for ki in xrange(r.beg/self.binsize, r.end/self.binsize+1):
                k = (r.chrm, ki)
                self._bin_add_record(k, r)

    def get(self, chrm, beg, end):

        """ get record overlapping beg and end """
        
        if not end: end = beg
        chrm = normalize_chrm(chrm)
        kbeg = int(beg) / self.binsize
        kend = int(end) / self.binsize
        records = []
        for ki in xrange(kbeg, kend+1):
            k = (chrm, ki)
            if k in self.bin2record:
                for r in self.bin2record[k]:
                    print beg, end, r.beg, r.end
                    if _overlap(r.beg, r.end, beg, end) > 0:
                        if r not in records:
                            records.append(r)

        return records

    def __repr__(self):
        return '<Bed file with %dx%d records>' % self.data.shape

    def sort(self, force=False):
        
        if self.issorted and not force:
            return
        self.data.sort(columns=['chrm','beg','end'], inplace=True)
        self.issorted = True

    def intersect(self, b):

        """ assume input beds are all sorted """

        self.sort()
        b.sort()

        selfvals = self.data.columns[3:]
        bvals = b.data.columns[3:]

        iter2 = b.data.iterrows()
        next_r2 = nextrec(iter2)

        bout_rowlist = []
        _r2s = []
        for index1, r1 in self.data.iterrows():

            # remove all r2 with end before r1.beg
            _r2s2 = []
            for r2 in _r2s:
                if r1.chrm == r2.chrm and r2.end >= r1.beg:
                    _r2s2.append(r2)
            _r2s = _r2s2

            # flush additional r2
            while next_r2 is not None and next_r2.chrm < r1.chrm:
                next_r2 = nextrec(iter2)

            # get all r2 with beg before t1.end
            while (next_r2 is not None and
                   next_r2.chrm == r1.chrm and
                   next_r2.beg <= r1.end):

                if next_r2.end >= r1.beg:
                    _r2s.append(next_r2)

                next_r2 = nextrec(iter2)

            for r2 in _r2s:
                if r1.chrm == r2.chrm and r2.beg > r1.end:
                    break
                if overlap(r1, r2) > 0:
                    bout_rowlist.append([r1.chrm, max(r1.beg, r2.beg), min(r1.end, r2.end), overlap(r1, r2)]
                                        +[r1.beg, r1.end]+[r1[v] for v in selfvals]
                                        +[r2.beg, r2.end]+[r2[v] for v in bvals])

        b_out = Bed()
        b_out.data = pd.DataFrame(bout_rowlist)
        b_out.data.rename(columns={0:'chrm', 1:'beg', 2:'end'}, inplace=True)

        return b_out

    def exclude(self, b):

        """ assume input beds are all sorted """

        self.sort()
        b.sort()

        selfvals = self.data.columns[3:]
        bvals = b.data.columns[3:]

        iter2 = b.data.iterrows()
        next_r2 = nextrec(iter2)

        bout_rowlist = []
        _r2s = []
        for index1, r1 in self.data.iterrows():

            # remove all r2 with end before r1.beg
            _r2s2 = []
            for r2 in _r2s:
                if r1.chrm == r2.chrm and r2.end >= r1.beg:
                    _r2s2.append(r2)
            _r2s = _r2s2

            while next_r2 is not None and next_r2.chrm < r1.chrm:
                next_r2 = nextrec(iter2)

            # get all r2 with beg before t1.end
            while (next_r2 is not None and
                   next_r2.chrm == r1.chrm and
                   next_r2.beg <= r1.end):

                if next_r2.end >= r1.beg:
                    _r2s.append(next_r2)

                next_r2 = nextrec(iter2)

            intersect = False
            for r2 in _r2s:
                if r1.chrm == r2.chrm and r2.beg > r1.end:
                    break
                if overlap(r1, r2) > 0:
                    intersect = True

            if not intersect:
                bout_rowlist.append(r1)

        b_out = Bed()
        b_out.data = pd.DataFrame(bout_rowlist)
        b_out.data.rename(columns={0:'chrm', 1:'beg', 2:'end'}, inplace=True)

        return b_out
    
    def complement(self, fai_path):

        chrm2len = {}
        with open(fai_path) as fh:
            for line in fh:
                fields = line.strip().split('\t')
                chrm2len[fields[0]] = int(fields[1])

        self.sort()

        pbeg = 0
        pchrm = None
        bc_rowlist = []
        chrmseen = []
        for index, r in self.data.iterrows():

            if r.chrm not in chrmseen:
                chrmseen.append(r.chrm)

            if r.chrm != pchrm:
                if pchrm is not None:
                    if pbeg < chrm2len[pchrm]:
                        bc_rowlist.append([pchrm, pbeg, chrm2len[pchrm]])
                pchrm = r.chrm
                pbeg = 0
            if r.beg > pbeg:
                bc_rowlist.append([r.chrm, pbeg, r.beg])
                pbeg = r.end
        if pchrm is not None:
            if pbeg < chrm2len[pchrm]:
                bc_rowlist.append([pchrm, pbeg, chrm2len[pchrm]])

        for chrm, chrmlen in chrm2len.iteritems():
            if chrm not in chrmseen:
                bc_rowlist.append([chrm, 0, chrmlen])

        bc = Bed()
        bc.data = pd.DataFrame(bc_rowlist, columns=['chrm', 'beg', 'end'])

        return bc
