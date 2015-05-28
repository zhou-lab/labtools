import pandas as pd

# def cluster_sample_by_variable_probes(betas):

def clean_450k(df, probe_fn='/Users/wandingzhou/projects/hs-tcga/data/2015_03_05_TCGA_450/450k_probes'):

    probe_loc = pd.read_table(probe_fn,index_col=3,header=None, names=['chrm','beg','end','gene'])
    # remove X,Y chromosome
    df = df[(~probe_loc.chrm.isin(['chrX','chrY']))[df.index]]
    # remove NA
    df = df.dropna()
    # remove SNP
    df = df[df.index.map(lambda x: not x.startswith('rs'))]

    print "There are %d probes and %d samples" % df.shape
    return df

""" select probes that are methylated in some samples and unmethylated in others """
def polarized(df, upthres=0.7,dwthres=0.3, kind='any', minsupp=1):

    def _isvar(row):
        hi = (row > upthres).sum()
        lo = (row < dwthres).sum()
        if hi >= minsupp and lo >= minsupp:
            if kind == 'any':
                return True
            elif kind == 'all' and hi + lo == row.shape[0]:
                return True
            else:
                return False
        else:
            return False

    dfv = df[df.apply(_isvar, axis=1)]
    print 'Selected %d variable probes from %d samples' % dfv.shape
    return dfv

def select_std(df, topn=1000, std_thres=None):

    std = df.std(axis=1)
    if std_thres is None:
        std.sort(ascending=False)
        dfv = df[std[:topn]].copy()
    else:
        dfv = df[std > std_thres].copy()

    print 'Selected %d variable probes from %d samples' % dfv.shape
    return dfv
    
def dichotomize(df):
    pass

def uniformly_methylated(df, thres=0.8):

    dfm = df[df.apply(lambda x: (x>thres).all())]

    print 'Selected %d methylated probes (>%1.3f) from %d samples' % (dfm.shape[0], thres, dfm.shape[1])
    return dfm

def uniformly_unmethylated(df, thres=0.2):

    dfu = df[df.apply(lambda x: (x<thres).all())]

    print 'Selected %d methylated probes (<%1.3f) from %d samples' % (dfu.shape[0], thres, dfu.shape[1])
    return dfu

""" filter uniformly methylated and uniformly unmethylated probes """
def filter_uniform(df, maxbeta=0.7, minbeta=0.3, maxsupp=0.9):

    maxsuppn = float(df.shape[1]) * maxsupp
    def _is_nonuniform(row):
        return not (((row > maxbeta).sum() > maxsuppn) or ((row < minbeta).sum() > maxsuppn))

    dfv = df[df.apply(_is_nonuniform, axis=1)]
    print 'Filtered %d uniform probes from %d' % (df.shape[0] - dfv.shape[0], df.shape[0])

    return dfv

def split_tumor_normal(df):

    dft = df.loc[:,df.columns.map(lambda x:x[13]=='0')]
    dfn = df.loc[:,df.columns.map(lambda x:x[13]=='1')]

    return dft, dfn

# class BloodTest:

#     def __init__(self):

#         self.sample = 
        
def discretize(df, upthres=0.8, dwthres=0.2):

    df_bin = pd.DataFrame()
    df_bin = df.copy()
    def _discretize(v):
        if v >= upthres:
            return 1
        elif v <= dwthres:
            return 0
        elif v <= 0.5:
            return 2
        else:
            return 3
    df_bin = df_bin.applymap(_discretize)
    return df_bin
