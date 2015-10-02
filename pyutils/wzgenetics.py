import wzcore
import faidx

class GeneticElement(object):

    def __init__(self, flank=1000):
        self.flank = flank
        
    def __len__(self):
        return self.end - self.beg
    
    def __repr__(self):
        
        return '<%s: %d bp>' % (self.alutype, self.__len__())
    
    def cgcnt(self):
        return self.seq().count('CG')
    
    def cgdensity(self):
        return float(self.cgcnt()) / float(self.__len__())

    def seq(self):
        return self.area[self.flank:self.flank+len(self)]

    def area_cgdensities(self, step=20, window=100):

        cgd = []
        for start in xrange(0, len(self.area)-window, step):
            cgd.append(self.area[start:start+window].count('CG'))
        
        return cgd

class Alu(GeneticElement):

    def __init__(self, flank=1000):

        super(Alu, self).__init__(flank=flank)

    def area_beg(self):
        if self.strand == '+':
            return self.beg - self.flank
        else:
            return int(round(self.end-1.5*self.flank))
        
    def area_end(self):
        if self.strand == '+':
            return int(round(self.beg+1.5*self.flank))
        else:
            return self.end + self.flank

def load_alu_and_seqs():

    refgenome = faidx.RefGenome('/Users/wandingzhou/references/hg19/hg19.fa')
    alus = {}
    wzcore.err_print_sig()
    for i,line in enumerate(wzcore.opengz('/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/alu.bed.gz')):
        if i%100000 == 0:
            wzcore.err_print_m(' %d' % i)

        fields = line.strip().split('\t')
        alu = Alu()
        alu.chrm = fields[0]
        if alu.chrm.find('_')>0:
            continue
        alu.beg = int(fields[1])
        alu.end = int(fields[2])
        alu.strand = fields[3]
        alu.alutype = fields[4]
        try:
            if alu.strand == '+':
                alu.area = refgenome.fetch_sequence(alu.chrm, alu.area_beg(), alu.area_end()).upper()
            else:
                alu.area = wzcore.reverse(refgenome.fetch_sequence(alu.chrm, alu.area_beg(), alu.area_end()).upper())
        except IndexError:      # Alu at chromosome boundaries, ignore
            continue
        alus[(alu.chrm,alu.beg,alu.end)] = alu

    wzcore.err_print_m('\n')
    wzcore.err_print('Loaded %d Alus' % len(alus))
    return alus

class CGI(GeneticElement):

    def __init__(self, flank=1000):

        super(CGI, self).__init__(flank=flank)

    def __repr__(self):
        
        return '<CGI %s: %d bp>' % (self.cgitype, self.__len__())
    
    def area_beg(self):
        return self.beg - self.flank
        
    def area_end(self):
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
        cgi.area = refgenome.fetch_sequence(cgi.chrm, cgi.area_beg(), cgi.area_end()).upper()
        cgis.append(cgi)
        
    wzcore.err_print('Loaded %d CGIs' % len(cgis))
    return cgis

def load_tss_and_seqs():

    import pandas as pd
    tss_table = pd.read_table('/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/hg19_mRNA_tss', header=None)
    tss_table.columns = ['chrm','tss','strand','gene','_cnt','transname']
    tss_table.index = tss_table['chrm']+":"+tss_table['tss'].map(str)
    tss_table = tss_table.groupby(level=0).first()
    cgi_anno = pd.read_table('/Users/wandingzhou/projects/pj-mm/2015-04-23-alu/hg19_mRNA_tss_1k1k_cgi_uniq',sep='\t', index_col='chrmpos')
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
        _poses = [np.nan]*len(poses)
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
