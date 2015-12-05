import pandas as pd
import numpy as np
import wzcore
# def cluster_sample_by_variable_probes(betas):

def makedict(df, k, v):

    d = df[v].copy()
    d.index = df[k]
    d = d.groupby(level=0).last()
    return d

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

def check_sample(sample):

    Hui_annot = pd.read_csv('/Users/wandingzhou/projects/hs-tcga/data/2015_03_23_Hui_annotation/sampleAnnotSubWB20130619.tsv', index_col=1, error_bad_lines=False, sep='\t')

    with pd.option_context('display.max_rows', 999):
        print Hui_annot.loc[sample]


def remove_snps(df):

    df = df[df.index.map(lambda x: not x.startswith('rs'))].copy()

    return df

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
    """ select probes that are methylated in some samples and unmethylated in others
    "any" means more than minsupp probes have hypo and hyper meth
    "all" all probes have either hypo or hyper meth
    """

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

def take_segment_mean(df, probe2seg, min_support=10, quiet=False):

    _df_seg = df.copy()
    if isinstance(_df_seg, pd.Series):
        _df_seg = _df_seg.to_frame()
    _df_seg['seg'] = probe2seg.loc[df.index]
    _df_seg_gb = _df_seg.groupby('seg')
    df_seg_count = _df_seg_gb.count().iloc[:,0]
    df_seg_mean = _df_seg_gb.mean()[df_seg_count >= min_support]
    if not quiet:
        wzcore.err_print('There are %d segments well supported.' % df_seg_mean.shape[0])

    return df_seg_mean

class MutData():

    def __init__(self):

        self.sample2muts = {}
        self.genes = set()
        with open('/Users/wandingzhou/projects/hs-tcga/data/2015_06_17_TCGA_mutations/merged_maf') as fh:
            for line in fh:
                fields = line.strip().split()
                if fields[2] in ['Silent', 'RNA']:
                    continue
                sample = fields[4][:12]
                if sample not in self.sample2muts:
                    self.sample2muts[sample] = []
                mut = fields[1]
                self.sample2muts[sample].append((mut, fields[2], fields[3]))
                self.genes.add(mut)

        wzcore.err_print('Loaded %d genes and %d samples' % (len(self.genes), len(self.sample2muts)))

        return

    def samples(self):

        return sample2muts.keys()

    def mutstat(self, samples, gene, verbose=True):

        t = []
        cnt = 0
        cnt1 = 0
        if gene not in self.genes:
            wzcore.err_print('Gene ID not identified %s' % gene)
            return None
        
        for s in samples:
            s = s[:12]
            if s in self.sample2muts:
                mut = False
                for g, mt, mt1 in self.sample2muts[s]:
                    if g == gene:
                        mut = True
                        cnt1 += 1
                        break
                t.append(mut)
                cnt += 1
            else:
                t.append(np.nan)

        if verbose:
            wzcore.err_print('Identified info from %d/%d samples (%d muts).' % (cnt, len(samples), cnt1))

        return pd.Series(t, index=samples)

    def itergene(self, samples, min_muts=5, verbose=False):

        for g in self.genes:
            s = self.mutstat(samples, g, verbose=verbose)
            if s is not None and s.sum() >= min_muts:
                yield g, s

    def associate_continuous(self, ts, prefix):

        import scipy.stats as stats
        self.ts = ts
        self.pvals = []
        self.foldcs = []
        self.genes_select = []
        for g, ss in self.itergene(ts.index):
            U, pval = stats.mannwhitneyu(ts[ss == True],ts[ss == False])
            foldc = (float(ts[ss == True].median()) / ts[ss == False].median())
            self.pvals.append(-np.log2(pval))
            self.foldcs.append(np.log2(foldc))
            self.genes_select.append(g)
            # if len(self.genes) > 100:
            # break

        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(self.foldcs, self.pvals, edgecolor='none', alpha=0.5)
        plt.xlabel('mutation / nonmutation')
        plt.ylabel('p-value')
        plt.savefig(prefix+'_volcano.png')

        self.associate_continuous_outliers(prefix+'_outlier.png')

    def associate_continuous_outliers(self, outlier_fn, fw=20.0, printn=10):
        
        balanced_up = [((pval+foldc*fw)/2.0,g,pval,foldc)
                       for pval, foldc, g in zip(self.pvals, self.foldcs, self.genes_select)
                       if (not pd.isnull(pval)) and (not pd.isnull(foldc))]
        balanced_dw = [((pval-foldc*fw)/2.0,g,pval,foldc)
                       for pval, foldc, g in zip(self.pvals, self.foldcs, self.genes_select)
                       if (not pd.isnull(pval)) and (not pd.isnull(foldc))]

        import wzplotlib

        tss = self.ts.copy()
        tss.sort()
        cbs = [wzplotlib.WZCbar(tss, continuous=True)]
        balanced_up.sort(reverse=True)
        wzcore.err_print('upside:')
        for i in xrange(min(printn, len(self.genes_select))):
            b, g, pval, foldc = balanced_up[i]
            cbs.append(wzplotlib.WZCbar(self.mutstat(tss.index, g, verbose=False)[tss.index], title=g+' [up]'))
            wzcore.err_print('%s\tpval:%1.2f\tfoldc:%1.2f\tbalanced:%1.2f' % (g, pval, foldc, b))

        balanced_dw.sort(reverse=True)
        wzcore.err_print('\ndownside:')
        for i in xrange(min(printn, len(self.genes_select))):
            b, g, pval, foldc = balanced_dw[i]
            cbs.append(wzplotlib.WZCbar(self.mutstat(tss.index, g, verbose=False)[tss.index], title=g+' [down]'))
            wzcore.err_print('%s\tpval:%1.2f\tfoldc:%1.2f\tbalanced:%1.2f' % (g, pval, foldc, b))

        wzplotlib.row_stack_layout(cbs, figfile=outlier_fn)

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


def associate_sample_continuous(df, s, mut=True, cnv=True):

    cbs = [s]

    if mut:
        for iname, idat in df.iteritems():
            if iname.startswith('mut_'):
                cbs.append(idat[s.index])

    if cnv:
        for iname, idat in df.iteritems():
            if iname.startswith('cnv_'):
                cbs.append(idat[s.index])

    wzplotlib.row_stack_layout(cbs)

class ExpData():
    
    def __init__(self, genes=None):

        self.rsems, self.cancer_types = data_load_rnaseq_all(genes=genes)

    def _expstat(self, ts, dgene):
        sample2v = {}
        for sample, v in dgene.iteritems():
            if not np.isnan(v):
                sample2v[sample[:12]] = v
        ts2v = {}
        for sample, v in ts.iteritems():
            if not np.isnan(v):
                ts2v[sample[:12]] = v
        _dgene = []
        _ts = []
        for s in sample2v:
            if s in ts2v:
                _dgene.append(sample2v[s])
                _ts.append(ts2v[s])
        return _dgene, _ts

    def expstat(self, ts, gene):

        return self._expstat(ts, self.rsems.loc[gene])

    def itergene(self, ts):

        for g, dgene in self.rsems.iterrows():
            _dgene, _ts = self._expstat(ts, dgene)
            yield g, _dgene, _ts
        
    def associate_continuous(self, ts, prefix):

        import scipy.stats as stats
        self.rhos = []
        self.pvals = []
        self.genes_select = []
        self.ts = ts

        for g, _dgene, _ts in self.itergene(ts):
            rho, pval = stats.spearmanr(_dgene, _ts)
            self.rhos.append(rho)
            self.pvals.append(-np.log2(pval))
            self.genes_select.append(g)
            # if len(self.rhos) > 100:
            # break
        
        import matplotlib.pyplot as plt
        plt.figure()
        plt.scatter(self.rhos, self.pvals, edgecolor='none', alpha=0.5)
        plt.xlabel("spearman's rho")
        plt.ylabel('p-value')
        plt.savefig(prefix+'_volcano.png')

        self.associate_continuous_outliers(prefix+'_outlier.png')

    def associate_continuous_outliers(self, fn, fw=10, printn=10):

        import matplotlib.pyplot as plt
        plt.figure(figsize=(10,10))
        plt.subplots_adjust(hspace=1)
        
        toplist = sorted(zip(self.pvals, self.rhos, self.genes_select), reverse=True)
        n = min(printn, len(toplist))
        for i in xrange(n):
            pval, rho, gene = toplist[i]
            wzcore.err_print('%s\tpval:%1.2f\trho:%1.2f' % (gene, pval, rho))
            plt.subplot(5,2,i+1)
            _dgene, _ts = self.expstat(self.ts, gene)
            plt.scatter(np.log2(_dgene), _ts, edgecolor='none', alpha=0.5, s=4)
            plt.xlabel('log2('+gene+' %1.2f)' % rho)

        plt.savefig(fn, bbox_inches='tight')

"""" data load subroutines """

def data_load_tissue(source):

    if source == 'Laird':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_06_03_AML_normal_sorted/GSE49618_betas.tsv')
        samples = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_06_03_AML_normal_sorted/samples', index_col='barcode')
        betas.columns = betas.columns.map(lambda x: 'blood_'+samples.loc[x.split('_')[0],'name'])

    if source == 'Encode':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE40699_ENCODE/GSE40699_betas.tsv')
        samples = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE40699_ENCODE/samples', index_col='barcode')
        betas.columns = samples.loc[betas.columns.map(lambda x:x.split('_',2)[2]), 'short']+'_'+samples.loc[betas.columns.map(lambda x:x.split('_',2)[2]), 'cellline']

    if source == 'Bonder':      # muscle and fat

        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/Bonder2014_BMCGenomics/GSE61454_severely_obsese/betas.tsv')
        samples = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/Bonder2014_BMCGenomics/GSE61454_severely_obsese/samples',header=None,index_col=0,names=['barcode','sample'])
        betas.columns = samples.loc[betas.columns.map(lambda x:x.split('_',1)[0]),'sample']+"_"+betas.columns.map(lambda x:x.split('_',1)[0])
        betas = betas.loc[:,~betas.columns.map(lambda x:x.startswith('Liver'))] # exclude liver

    if source == 'Slieker':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE48472_Slieker_2013_EpigeneticsAndChromatin/betas.tsv')
        samples = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE48472_Slieker_2013_EpigeneticsAndChromatin/samples',header=None,index_col=0,names=['barcode','sample'])
        betas.columns = samples.loc[betas.columns.map(lambda x:x.split('_',1)[0]),'sample']+"_"+betas.columns.map(lambda x:x.split('_',1)[0])

    if source == 'Wong':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_06_03_AlexWong_skin/AlexWong_skin_betas.tsv')
        samples = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_06_03_AlexWong_skin/samples', index_col='barcode')
        betas.columns = 'skin_'+betas.columns.to_series()

    if source == 'Lawlor1133':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_06_03_Lawlor_tumor_lungfibroblast/s1133_betas.tsv')
        names = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_06_03_Lawlor_tumor_lungfibroblast/1133samples.csv.unix', index_col='Complete_Barcode', sep='\t')
        betas.columns = names.loc[betas.columns,'GROUP_NAME']+'_'+betas.columns
        
    if source == 'Guintivano':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE41826_Brain/betas.tsv')
        mask_snp_probes(betas)
        names = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE41826_Brain/samples.csv',index_col='barcode')
        betas.columns = 'brain_'+names.loc[betas.columns,'sample'].map(str)+'_'+betas.columns

    if source == 'Wagner':
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_03_18_tumor_purity/GSE52025_Wagner_fibroblast/GSE52025_betas.tsv',sep='\t')
        betas.columns = 'fibroblast_'+betas.columns.map(str)

    if source == 'Reinus':      # blood
        betas = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_04_10_sorted_cell_population/blood_beta.tsv', sep='\t')
        names = pd.read_table('/Users/wandingzhou/projects/hs-tcga/data/2015_04_10_sorted_cell_population/Sorted_Blood/sample_sheet_IDAT.csv.unix.tsv',index_col='barcode')
        betas.columns = names.loc[betas.columns,'Type']+'-'+betas.columns
        
    wzcore.err_print("Loaded %d samples" % betas.shape[1])
    
    return betas

def data_create_tissue_sample():
    samples = pd.read_table('/Users/wandingzhou/projects/hs-tcga/2015_07_14_purity/2015_07_14_sample_select.txt', index_col="sample")
    betas = pd.DataFrame()
    for k, v in samples.sort('source').groupby('source'):
        _betas = data_load_tissue(v['source'][0])
        _betas2 = _betas[_betas.columns.intersection(v.index)]
        betas = pd.concat([betas, _betas2], axis=1)
        wzcore.err_print("From %s chose %d samples" % (v['source'][0], _betas2.shape[1]))

    wzcore.err_print("Loaded %d samples" % betas.shape[1])

    return betas, samples
    

def data_load_commonPMD(commonfn='/Users/wandingzhou/projects/hs-tcga/2015_10_12_Huy_segments/common_probes'):
    # commonfn='/Users/wandingzhou/projects/hs-tcga/2015_05_04_pmd/stringent_common_probes'):

    commons = pd.read_table(commonfn, sep='\t',header=None, index_col=3)
    commons.rename(columns={5:'domain', 4:'domaintype', 0:'chrm',1:'beg',2:'end'},
                   inplace=True)
    seg2MD = makedict(commons, 'domain', 'domaintype')
    return commons, seg2MD

def data_load_WGBS_betas_common_PMD(commonfn='/Users/wandingzhou/projects/hs-tcga/2015_10_12_Huy_segments/common_probes'):
    # commonfn='/Users/wandingzhou/projects/hs-tcga/2015_05_04_pmd/stringent_common_probes'):

    commons, seg2MD = data_load_commonPMD(commonfn)
    betas, cancer_types = data_load_WGBS_betas(commons.index)
    return betas, cancer_types, commons, seg2MD

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
    betas = betas.groupby(level=0, axis=1).mean()
    wzcore.err_print('Loaded %d probes and %d samples' % betas.shape)

    return betas, cancer_types

def data_load_rnaseq_all(genes=None):

    data_home = '/Users/wandingzhou/projects/hs-tcga/data/2015_04_30_TCGA_rnaseq/dat/'
    import os
    cancer_types = [pkl.strip('.pkl') for pkl in os.listdir(data_home)]
    return data_load_rnaseq_samples(cancer_types, genes)

def data_load_rnaseq_samples(_cancer_types, genes=None):

    if not isinstance(_cancer_types, list):
        _cancer_types = [_cancer_types]

    data_home = '/Users/wandingzhou/projects/hs-tcga/data/2015_04_30_TCGA_rnaseq/dat/'
    _rsems = []
    cancer_types = pd.Series()
    
    wzcore.err_print_sig()
    for cancer_type in _cancer_types:
        wzcore.err_print_m(' '+cancer_type)
        dat_fn = cancer_type+'.pkl'
        _rsem = pd.read_pickle(data_home+'/'+dat_fn)
        if genes is None:
            _rsems.append(_rsem)
        else:
            _rsems.append(_rsem.loc[genes,])
        cancer_types = cancer_types.append(pd.Series([cancer_type]*_rsem.shape[1], index=_rsem.columns))

    rsems = pd.concat(_rsems, axis=1)
    wzcore.err_print_m('\n')
            
    cancer_types = cancer_types.groupby(level=0).last()
    wzcore.err_print('Loaded %d genes and %d samples' % rsems.shape)
    
    return rsems, cancer_types

def probe_select1(betas, v, upq=0.75, loq=0.25, mindelta=0.1):

    inind = betas.columns.isin(v.index)
    _select1 = []
    _select2 = []

    hyperrows = []
    hyporows = []
    for i, row in betas.iterrows():
        c = row[inind]
        cb = row[~inind]
        hyperrows.append((c.min()-cb.quantile(upq),i))
        hyporows.append((cb.quantile(loq)-c.max(),i))

    _select1 = [(i,j) for i,j in sorted(hyporows, reverse=True)[:100] if i>mindelta]
    _select2 = [(i,j) for i,j in sorted(hyperrows, reverse=True)[:100] if i>mindelta]

    wzcore.err_print("selected %d hypo and %d hyper probes." % (len(_select1), len(_select2)))

    return list(set([j for i,j in _select1]) | set([j for i,j in _select2]))

def probe_select_pairwise(betas, v1, v2, mindelta=0.5, upq=0.75, loq=0.25):

    inind1 = betas.columns.isin(v1.index)
    inind2 = betas.columns.isin(v2.index)

    select1 = []
    select2 = []
    for i, row in betas.iterrows():
        c1 = row[inind1]
        c2 = row[inind2]
        if c2.quantile(loq) - c1.quantile(upq) > mindelta: # hypo
            select1.append(i)
        if c1.quantile(loq) - c2.quantile(upq) > mindelta: # hyper
            select2.append(i)

    wzcore.err_print("selected %d hypo and %d hyper." % (len(select1), len(select2)))
    return select1, select2

def mark_alternate(series):
    segmarks = []
    prev = None
    mark = True
    for i, v in series.iteritems(): # uniseg['segments'][betas_uniseg_tumor.index].iteritems():
        if v != prev:
            mark = not mark
            prev = v
        segmarks.append(mark)

    return segmarks
