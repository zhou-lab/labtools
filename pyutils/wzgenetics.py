import wzcore
import faidx

refgenomehg19 = faidx.RefGenome('/Users/wandingzhou/references/hg19/hg19.fa')

class GeneticElement(object):

    def __init__(self, flank=1000, flank1=None, flank2=None):
        self.flank1 = flank     # upstream (gene order)
        self.flank2 = flank     # downstream (gene order)
        if flank1:
            self.flank1 = flank1
        if flank2:
            self.flank2 = flank2
        
    def __len__(self):
        return self.end - self.beg
    
    def __repr__(self):
        
        return '<Element: %s:%d-%d>' % (self.chrm, self.beg, self.end)

    def rep_coord(self):

        return '%s:%d-%d' % (self.chrm, self.beg, self.end)

    def seq_beg(self):
        """ note that beg, end, seq_beg(), seq_end() are all in reference coordinates """
        if self.strand == '+':
            return self.beg-self.flank1
        else:
            return int(round(self.end-self.flank2))
        
    def seq_end(self):
        """ note that beg, end, seq_beg(), seq_end() are all in reference coordinates """
        if self.strand == '+':
            return int(round(self.beg+self.flank2))
        else:
            return self.end + self.flank1

    def ownseq(self):
        return self.seq[self.flank1:self.flank1+len(self)]

    def cpgcnt(self):
        return self.ownseq().count('CG')
    
    def cpgdensity(self):
        return float(self.cpgcnt()) / float(self.__len__())

    def upstream(self,l):
        return self.seq[self.flank1-l:self.flank1]

    def dwstream(self,l):
        return self.seq[self.flank1+len(self):self.flank1+len(self)+l]
    
    def flanking_cpgdensity(self, flank=1000):

        if flank>=self.flank1 or flank+len(self)>=self.flank2:
            raise Exception('need to get larger sequence')

        seq_up = self.upstream(flank)
        seq_dw = self.dwstream(flank)
        if len(seq_up) + len(seq_dw) > 0:
            return float(seq_up.count('CG') + seq_dw.count('CG')) / (len(seq_up)+len(seq_dw))
        else:
            return np.nan

    def area_cpgdensities(self, step=50, window=100, area_flank1=1000, area_flank2=1300):

        cgd = []
        window_begs = []
        for start in xrange(self.flank1-area_flank1, self.flank1+area_flank2-window, step):
            window_begs.append(start-self.flank1)
            cgd.append(self.seq[start:start+window].count('CG'))
        
        return window_begs, cgd

class TE(GeneticElement):

    def __init__(self, flank=1000):
        super(TE, self).__init__(flank=flank)
        self.seq = None
        self.tetype = 'RepElement'

    def __repr__(self):
        return '<%s: %s:%d-%d>' % (self.tetype, self.chrm, self.beg, self.end)

    def load_seq(self, refgenome=refgenomehg19, flank=1000, flank1=None, flank2=None):
        self.flank1 = flank1 if flank1 else flank
        self.flank2 = flank2 if flank2 else flank

        if self.strand == '+':
            self.seq = refgenome.fetch_sequence(self.chrm, self.seq_beg(), self.seq_end()).upper()
        else:
            self.seq = wzcore.reverse_complement(refgenome.fetch_sequence(self.chrm, self.seq_beg(), self.seq_end()).upper())

        if len(self.seq) != self.seq_end()-self.seq_beg()+1:
            self.seq = None    # set to None if sequence retrieval error

    def get_microenv(self, d=100):

        """
        self.rmskbed must be tabix indexed, still this takes 14hr to load all TE environment
        this is pretty slow
        """
        
        import tabix
        te_fh = tabix.open(self.rmskbed)

        self.upstream = []
        winend = self.beg
        beg = self.beg
        while True:
            hits = list(te_fh.query(self.chrm, winend-1000, winend))
            winend -= 1000
            disconnected = False if hits else True
            for fields in reversed(hits):
                _chrm, _beg, _end, _strand, _tetype, _tetype2, _tetype3 = fields
                te = TE()
                te.chrm = _chrm
                te.beg = int(_beg)
                te.end = int(_end)
                te.strand = _strand
                te.tetype = _tetype
                te.tetype2 = _tetype2
                te.tetype3 = _tetype3
                te.dist = te.end - beg
                if -te.dist > d:
                    disconnected = True
                    break
                beg = te.beg
                self.upstream.append(te)
                
            if disconnected:
                break

        self.downstream = []
        winbeg = self.end+1
        end = self.end+1
        while True:
            hits = list(te_fh.query(self.chrm, winbeg, winbeg+1000))
            winbeg += 1000
            disconnected = False if hits else True
            for fields in hits:
                _chrm, _beg, _end, _strand, _tetype, _tetype2, _tetype3 = fields
                te = TE()
                te.chrm = _chrm
                te.beg = int(_beg)
                te.end = int(_end)
                te.strand = _strand
                te.tetype = _tetype
                te.tetype2 = _tetype2
                te.tetype3 = _tetype3
                te.dist = te.beg - end
                if te.dist > d:
                    disconnected = True
                    break
                end = te.end
                self.downstream.append(te)

            if disconnected:
                break

def _te_load_seqs(refgenome, te):

    if te.strand == '+':
        te.seq = refgenome.fetch_sequence(te.chrm, te.seq_beg(), te.seq_end()).upper()
    else:
        te.seq = wzcore.reverse_complement(refgenome.fetch_sequence(te.chrm, te.seq_beg(), te.seq_end()).upper())
    if len(te.seq) != te.seq_end()-te.seq_beg()+1:
        raise IndexError()

def load_te_and_seqs(rmskbed='/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/rmsk.bed.gz',
                     load_seq=False, tetype=None, tetype2=None, tetype3=None):

    refgenome = faidx.RefGenome('/Users/wandingzhou/references/hg19/hg19.fa')
    tes = {}
    wzcore.err_print_sig()
    for i,line in enumerate(wzcore.opengz(rmskbed)):
        if i%100000 == 0:
            wzcore.err_print_m(' %d' % i)

        fields = line.strip().split('\t')
        te = TE()
        te.chrm = fields[0]
        if te.chrm.find('_')>0:
            continue
        te.beg = int(fields[1])
        te.end = int(fields[2])
        te.rmskbed = rmskbed
        te.strand = fields[3]
        te.tetype = fields[4]
        te.tetype2 = fields[5]
        te.tetype3 = fields[6]

        if tetype is not None and te.tetype != tetype:
            continue

        if tetype2 is not None and te.tetype2 != tetype2:
            continue

        if tetype3 is not None and te.tetype3 != tetype3:
            continue
        
        if load_seq:
            try:
                _te_load_seqs(refgenome, te)
            except IndexError:      # TE at chromosome boundaries, ignore
                # te.seq == None
                pass

        tes[(te.chrm,te.beg,te.end)] = te

    wzcore.err_print_m('\n')
    wzcore.err_print('Loaded %d TEs' % len(tes))
    return tes

def te_load_seqs(tes):

    refgenome = faidx.RefGenome('/Users/wandingzhou/references/hg19/hg19.fa')

    tes2 = {}
    if isinstance(tes, dict):
        it = tes.itervalues()
    else:
        it = iter(tes)
    for te in it:
        try:
            _te_load_seqs(refgenome, te)
        except IndexError:      # TE at chromosome boundaries, ignore
            continue

        tes2[(te.chrm,te.beg,te.end)] = te

    return tes2

class CGI(GeneticElement):

    def __init__(self, flank=1000):

        super(CGI, self).__init__(flank=flank)

    def __repr__(self):
        
        return '<CGI %s: %d bp>' % (self.cgitype, self.__len__())
    
    def seq_beg(self):
        return self.beg - self.flank
        
    def seq_end(self):
        return int(round(self.beg+1.5*self.flank))

def load_cgi_and_seqs():

    cgis = []
    refgenome = faidx.RefGenome('/Users/wandingzhou/references/hg19/hg19.fa')
    for line in wzcore.opengz('/Users/wandingzhou/projects/hs-tcga/data/2015_03_24_cpg_island/TakaiJones/takai.jones.strict.bed.gz'):
        fields = line.strip().split('\t')
        cgi = CGI()
        cgi.chrm = fields[0]
        cgi.beg = int(fields[1])
        cgi.end = int(fields[2])
        cgi.cgitype = fields[3]
        cgi.seq = refgenome.fetch_sequence(cgi.chrm, cgi.seq_beg(), cgi.seq_end()).upper()
        cgis.append(cgi)
        
    wzcore.err_print('Loaded %d CGIs' % len(cgis))
    return cgis

def load_tss_and_seqs(tssfn='/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/hg19_mRNA_tss', cgifn='/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/hg19_mRNA_tss_1k1k_cgi_uniq'):

    import pandas as pd
    tss_table = pd.read_table(tssfn, header=None)
    tss_table.columns = ['chrm','tss','strand','gene','_cnt','transname']
    tss_table.index = tss_table['chrm']+":"+tss_table['tss'].map(str)
    tss_table = tss_table.groupby(level=0).first()
    if cgifn is not None:
        cgi_anno = pd.read_table(cgifn,sep='\t', index_col='chrmpos')
        tss_table['cgi'] = cgi_anno.loc[tss_table.index, 'CGI']

    return tss_table

def tss_load_te_state(tss_table, usetype=1):

    import tabix
    import numpy as np
    import pandas as pd
    te_fh = tabix.open('/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/rmsk.bed.gz')

    d_upstream = 2500
    d_downstream = 2500
    poses = range(-d_upstream, d_downstream+1)

    j = 0
    _allposes = []
    for tss_index, fields1 in tss_table.iterrows():
        j+=1
        # if j%1000 ==0:
        #     print j
        chrm, tss, strand, gene, _cnt, transnames, cgi = fields1
        chrm = fields1[0]
        tss = fields1[1]
        strand = fields1[2]
        gene = fields1[3]
        _cnt = fields1[4]
        transnames = fields1[5]
        if len(fields1)>6:
            cgi = fields1[6]
        
        _poses = ['NA']*len(poses)
        if strand == "+":
            d1 = d_upstream
            d2 = d_downstream
        else:
            d1 = d_downstream
            d2 = d_upstream

        for fields2 in te_fh.query(chrm, tss - d1, tss + d2):
            _chrm, _beg, _end, _strand, _te_type1, _te_type2, _te_type3 = fields2
            te_type = [_te_type1,_te_type2,_te_type3][usetype-1]
            _beg = int(_beg)
            _end = int(_end)
            for i in range(max(0, _beg - tss + d1), min(_end - tss, d2) + d1 + 1):
                if strand == '+':
                    _poses[i] = te_type
                else:
                    _poses[d1+d2 - i] = te_type

        _allposes.append(_poses)

    te_state = pd.DataFrame(_allposes, index=tss_table.index, columns=poses)

    return te_state

def tss_load_cpgdensity(tss_table):

    d_upstream = 2500
    d_downstream = 2500
    poses = range(-d_upstream, d_downstream+1)

    import faidx
    import numpy as np
    import pandas as pd
    refgenome = faidx.RefGenome('/Users/wandingzhou/references/hg19/hg19.fa')
    _CpG_density = []
    j = 0
    for k, f in tss_table.iterrows():
        chrm = f[0]
        tss = f[1]
        strand = f[2]

        j+=1
        # if j%1000 ==0:
        # print j
        if strand == "+":
            d1 = d_upstream
            d2 = d_downstream
        else:
            d1 = d_downstream
            d2 = d_upstream

        _cg = [np.nan]*(d1+d2+1)
        genome_seq = refgenome.fetch_sequence(chrm, tss-d1-200, tss+d2+200).upper()
        for i in xrange(d1+d2+1):
            seq = genome_seq[200+i-50:200+i+50]
            if strand == '+':
                _cg[i] = seq.count('CG')
            else:
                _cg[d1+d2-i] = seq.count('CG')

        _CpG_density.append(_cg)
    CpG_density = pd.DataFrame(_CpG_density, index=tss_table.index, columns=poses)

    return CpG_density
