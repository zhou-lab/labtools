import inspect
import sys

class WZData(object):
    pass

def tprint(array, out=None):
    if out is None:
        out = sys.stdout

    out.write('\t'.join(map(str, array))+'\n')

def wprint(msg):

    try:
        print(msg)
    except IOError as e:
        sys.exit()

def err_print(msg):
    fn = inspect.stack()[1][3]
    sys.stderr.write('[%s] %s\n' % (fn, str(msg)))

def err_print_sig():

    fn = inspect.stack()[1][3]
    sys.stderr.write('[%s]' % (fn))

def err_print_m(msg):

    sys.stderr.write(msg)
    sys.stderr.flush()

def pd_print_n(x, n):
    import pandas as pd
    pd.set_option('display.max_rows', n)
    print(x)
    pd.reset_option('display.max_rows')
    
def pd_print_full(x):
    import pandas as pd
    pd.set_option('display.max_rows', len(x))
    print(x)
    pd.reset_option('display.max_rows')

def opengz(fn,m='r'):
    
    if fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn,m)
    else:
        fh = open(fn,m)

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

def uniquify(seq):

    seq2 = []
    item2ind = {}
    for i in seq:
        if i not in item2ind:
            seq2.append(i)
            item2ind[i] = 1
        else:
            seq2.append(i+str(item2ind[i]))
            item2ind[i] += 1

    return seq2

def mode(l, dmin=None, dmax=None, bw=None, bins=100):

    import numpy as np
    minl = np.min(l)
    maxl = np.max(l)
    from scipy.stats import gaussian_kde
    if dmin is None:
        dmin = np.min(l)
    if dmax is None:
        dmax = np.max(l)
    if bw is None:
        bw = (dmax - dmin) / 20.

    kde = gaussian_kde(l, bw)

    support = np.linspace(dmin, dmax, bins)
    density = kde.evaluate(support)

    return support, density
