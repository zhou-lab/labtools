import wzcore
import faidx

class GeneticElement():

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

        super(self).__init__(flank=flank)

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
    alus = []
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
        alus.append(alu)

    wzcore.err_print_m('\n')
    wzcore.err_print('Loaded %d Alus' % len(alus))
    return alus

class CGI(GeneticElement):

    def __init__(self, flank=1000):

        super(self).__init__(flank=flank)

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

