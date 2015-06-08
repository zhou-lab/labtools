import pandas as pd
import numpy as np
import wzcore
# def cluster_sample_by_variable_probes(betas):

def check_probe(probe, probe_fn='/Users/wandingzhou/projects/hs-tcga/data/2015_03_05_TCGA_450/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz'):

    fh = wzcore.opengz(probe_fn)
    # skip 7 lines
    for i in xrange(7):
        fh.readline()

    item = fh.readline().strip().split(',')
    for line in fh:
        fields = line.strip().split(',')
        # p2info[fields[0]] = fields
        if fields[0] == probe:
            for i, field in enumerate(fields):
                if len(field)>0:
                    print item[i]+':', field
    
    return

def clean_450k(df, nahow='strong', probe_fn='/Users/wandingzhou/projects/hs-tcga/data/2015_03_05_TCGA_450/450k_probes', verbose=True):

    probe_loc = pd.read_table(probe_fn,index_col=3,header=None, names=['chrm','beg','end','gene'])

    # remove X,Y chromosome
    df = df[(~probe_loc.chrm.isin(['chrX','chrY']))[df.index]]
    # remove NA
    if nahow == 'strong':
        df = df.dropna(how='any')
    if nahow == 'weak':
        df = df.dropna(how='all')
    # remove SNP
    df = df[df.index.map(lambda x: not x.startswith('rs'))]

    if verbose:
        if len(df.shape) == 1:
            wzcore.err_print("Kept %d probes after remov" % df.shape[0])
        else:
            wzcore.err_print("Kept %d probes and %d samples" % df.shape)
        
    return df

def mask_snp_probes(df, mask_fn='/Users/wandingzhou/projects/hs-tcga/data/2015_03_05_TCGA_450/toMask.tsv', remove=False):

    probetomask = pd.read_table(mask_fn, index_col=0)
    df.loc[df.index.isin(probetomask.index)] = np.nan
    # wzcore.err_print('Masked %d probes from %d probes' % (df.shape[0] - dfm.shape[0], df.shape[0]))
    return

def polarized(df, upthres=0.7,dwthres=0.3, kind='any', minsupp=1):
    """ select probes that are methylated in some samples and unmethylated in others """

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
    wzcore.err_print('Selected %d variable probes from %d samples' % dfv.shape)
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

    dfm = df[df.apply(lambda x: (x>thres).all(), axis=1)]

    print 'Selected %d uniformly methylated probes (>%1.3f) from %d samples' % (dfm.shape[0], thres, dfm.shape[1])
    return dfm

def uniformly_unmethylated(df, thres=0.2):

    dfu = df[df.apply(lambda x: (x<thres).all(), axis=1)]

    print 'Selected %d uniformly unmethylated probes (<%1.3f) from %d samples' % (dfu.shape[0], thres, dfu.shape[1])
    return dfu

""" filter uniformly methylated and uniformly unmethylated probes """
def nonuniform(df, maxbeta=0.7, minbeta=0.3, maxsupp=0.95):

    maxsuppn = float(df.shape[1]) * maxsupp
    def _is_nonuniform(row):
        return not (((row > maxbeta).sum() > maxsuppn) or ((row < minbeta).sum() > maxsuppn))

    dfv = df[df.apply(_is_nonuniform, axis=1)]
    wzcore.err_print('Selected %d nonuniform probes from %d' % (dfv.shape[0], df.shape[0]))

    return dfv

def split_tumor_normal(df):

    dft = df.loc[:,df.columns.map(lambda x:x[13]=='0')]
    dfn = df.loc[:,df.columns.map(lambda x:x[13]=='1')]
    dfc = df.loc[:,df.columns.map(lambda x:x[13]=='2')]
    
    wzcore.err_print('Found %d tumor, %d normal, %d cell line and %d others.' % (dft.shape[1], dfn.shape[1], dfc.shape[1], df.shape[1]-dft.shape[1]-dfn.shape[1]-dfc.shape[1]))
    return dft, dfn, dfc

def get_tumor(df):
    dft, dfn, dfc = split_tumor_normal(df)
    return dft

def get_normal(df):
    
    dft, dfn, dfc = split_tumor_normal(df)
    return dfn

def get_cellline(df):
    
    dft, dfn, dfc = split_tumor_normal(df)
    return dfc


# class BloodTest:

#     def __init__(self):

#         self.sample =

def test_blood(target_df, figfn, blood_fn='/Users/wandingzhou/projects/hs-tcga/data/2015_04_10_sorted_cell_population/blood_beta.tsv', max_meth=0.85, min_meth=0.15):
    import pandas as pd
    blood_betas = pd.read_table(blood_fn)
    import matplotlib.pyplot as plt
    plt.figure()
    plt.subplot(2,1,1)
    blood_methylated = uniformly_methylated(blood_betas, thres=max_meth)
    plt.hist(target_df.loc[blood_methylated.index.intersection(target_df.index)].dropna(), bins=30)
    plt.yscale('log', nonposy='clip')
    plt.xlim(0,1)

    plt.subplot(2,1,2)
    blood_unmethylated = uniformly_unmethylated(blood_betas, thres=min_meth)
    plt.hist(target_df.loc[blood_unmethylated.index.intersection(target_df.index)].dropna(), bins=30)
    plt.yscale('log', nonposy='clip')
    plt.xlim(0,1)
    plt.xlabel('betas')

    plt.savefig(figfn)
    return

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


def data_load_WGBS_betas_common_PMD():

    commons = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_05_04_pmd/stringent_common_probes',sep='\t',header=None, index_col=3)
    commons.rename(columns={5:'segments', 4:'MD',0:'chrm',1:'beg',2:'end'}, inplace=True)

    betas, cancer_types = data_load_WGBS_betas(commons.index)
    return betas, cancer_types, commons

def data_load_WGBS_betas(probes=None):

    WGBS_samples = ["BLCA", "BRCA", "COAD", "GBM", "LUAD", "LUSC", "READ", "STAD", "UCEC", 'LAML']
    
    return data_load_samples(WGBS_samples, probes)

def data_load_pancan12(probes=None):

    pancan12_samples = ['OV', 'BRCA', 'GBM', 'KIRC', 'LAML', 'COAD', 'READ', 'HNSC', 'LUAD', 'LUSC', 'BLCA', 'UCEC']
    
    return data_load_samples(pancan12_samples, probes)

def data_load_samples(samples, probes=None):

    if not isinstance(samples, list):
        samples = [samples]

    data_home = '/Users/wandingzhou/projects/hs-tcga/data/2015_03_05_TCGA_450/dat/'
    _betases = []
    cancer_types = pd.Series()

    wzcore.err_print_sig()
    for cancer_type in samples:
        wzcore.err_print_m(' '+cancer_type)
        dat_fn = cancer_type+'.pkl'
        cancer_type = dat_fn.strip('.pkl')
        _betas = pd.read_pickle(data_home+'/'+dat_fn)
        if probes is None:
            _betases.append(_betas) # = pd.concat([betas, _betas], axis=1)
        else:
            _betases.append(_betas.loc[probes,])  # _betas = pd.concat([betas, _betas.loc[probes,]], axis=1)
        cancer_types = cancer_types.append(pd.Series([cancer_type]*_betas.shape[1], index=_betas.columns))

    betas = pd.concat(_betases, axis=1)
    wzcore.err_print_m('\n')

    # some cell line sample belong to multiple cancer type, choose the last cancer
    cancer_types = cancer_types.groupby(level=0).last()
    wzcore.err_print('Loaded %d probes and %d samples' % betas.shape)

    return betas, cancer_types

def filter_by_purity(betas, min_purity=0.8, keep_normal=True):

    Hui_annot = pd.read_csv('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/sampleAnnotSubWB20130619.tsv', index_col=1, error_bad_lines=False, sep='\t')
    # print Hui_annot.columns[Hui_annot.columns.map(lambda x: np.bool(re.search('purity', x)))]
    # there are two purities, ABSOLUTE.purity and abs.purity, ABSOLUTE.purity has more value (4776 vs 4468)
    # the two values agree on overlap
    # I hereby use ABSOLUTE.purity
    # 4239 samples with purity estimates
    purity_anno = Hui_annot['ABSOLUTE.purity']
    purity_anno = purity_anno[(purity_anno.notnull()) & (purity_anno != 1.0)]

    sample_is_pure = betas.columns.map(lambda x: (x in purity_anno and purity_anno[x] > min_purity) or x[13]=='1')
    betasv = betas.loc[:,sample_is_pure].copy()
    wzcore.err_print('Selected %d pure samples (>%1.2f) from %d samples.' % (betasv.shape[1], min_purity, betas.shape[1]))
    return betasv

def positive_diff(b1, b2, thres=0.3, verbose=True):

    select = (b1 - b2) > thres
    if verbose:
        wzcore.err_print('Selected %d probes' % select.sum())
    
    return b1.loc[select,], b2.loc[select,]

def negative_diff(b1, b2, thres=-0.3, verbose=True):

    select = (b1 - b2) < thres
    if verbose:
        wzcore.err_print('Selected %d probes' % select.sum())

    return b1.loc[select,], b2.loc[select,]


def local_maxima(x):

    n = len(x)
    maxima = []
    if x[0] >= x[1]:
        maxima.append(0)
    for i in xrange(1,n-1):
        if x[i] >= x[i-1] and x[i] >= x[i+1]:
            maxima.append(i)
    if x[n-1] >= x[n-2]:
        maxima.append(n-1)

    return maxima

def is_broad(x, i, span):

    n = len(x)
    for j in xrange(max(0,i-span), min(n-1,i+span)):
        if x[j] > x[i]:
            return False
    return True

def broad_local_maxima(x, span=20):

    broad_maxima = []
    for i in local_maxima(x):
        if is_broad(x, i, span):
            broad_maxima.append(i)
    return broad_maxima

def DMKSM_test(betas_t, betas_c, mindelta=0.2, plotfn=None, verbose=True):
    """ DMKSM: Differential Methylation of Kernel Smoothing Max """

    import scipy.stats as stats
    import numpy as np

    diffp = clean_450k(positive_diff(betas_c, betas_t, thres=mindelta, verbose=verbose)[1], verbose=verbose)
    if len(diffp) <= 1000:
        err_print('Too few differential probes (p). Consider lowering min delta?')
        return -1, -1, -1

    diffn = clean_450k(negative_diff(betas_c, betas_t, thres=-mindelta, verbose=verbose)[1], verbose=verbose)
    if len(diffn) <= 1000:
        err_print('Too few differential probes (n). Consider lowering min delta?')
        return -1, -1, -1

    kernelp = stats.gaussian_kde(diffp)
    kerneln = stats.gaussian_kde(diffn)

    grid = np.linspace(0,1,1000)

    maxps = broad_local_maxima(kernelp(grid))
    if maxps:
        maxp = grid[maxps[0]]
    else:
        maxp = max(grid, key=kernelp)

    maxns = broad_local_maxima(kerneln(grid))
    if maxns:
        maxn = grid[maxns[-1]]
    else:
        maxn = max(grid, key=kerneln)

    if verbose:
        err_print('Estimate from cm/tu probes: %1.3f' % (maxp, ))
        err_print('Estimate from cu/tm probes: %1.3f' % (1.0-maxn,))

    if plotfn is not None:
        import matplotlib.pyplot as plt
        xvals = np.linspace(0,1,100)

        plt.figure(figsize=(8,4))
        plt.subplot(1,2,1)
        plt.title('contaminant methylated')
        plt.hist(diffp,normed=True,color='b',alpha=0.5)
        plt.plot(xvals, kernelp(xvals),color='r',lw=2)

        plt.subplot(1,2,2)
        plt.title('contaminant unmethylated')
        plt.hist(diffn,normed=True,color='b',alpha=0.5)
        plt.plot(xvals, kerneln(xvals),color='r',lw=2)
        
        plt.savefig(plotfn)

    estimate = (maxp + 1-maxn)/2.0
    if estimate + mindelta > 0.85:
        err_print('Warning, min delta is too large.')

    if (maxp - 1 + maxn) / maxp > 0.1:
        err_print('Warning, estimate error greater than 10%')

    return estimate, maxp, 1-maxn


def DMKSM_delta_scan(betas_t, betas_c, plot_fn):
    import numpy as np
    import matplotlib.pyplot as plt
    
    xs = []
    allests = []
    pests = []
    nests = []
    for x in np.arange(0.01,0.8,0.02):
        est, ep, en = DMKSM_test(betas_t, betas_c, mindelta=x, verbose=False)
        if est>0:
            xs.append(x)
            allests.append(est)
            pests.append(ep)
            nests.append(en)

    plt.figure()
    plt.plot(xs, allests, lw=3, alpha=0.5)
    plt.plot(xs, pests, lw=3, alpha=0.5)
    plt.plot(xs, nests, lw=3, alpha=0.5)
    plt.plot([0,1],[1,0], ls='dashed')
    plt.xlabel('min delta')
    plt.ylabel('contamination estimate')
    plt.savefig(plot_fn)
