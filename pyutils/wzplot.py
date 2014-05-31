#!/RIS/HPC_apps/AMD/python/python-2.7.2/bin/python

import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
import matplotlib.cm as cm
import itertools
import numpy as np

colorcycle = ['DarkSlateGray', 'Brown', 'Burlywood', 'yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'y', 'k', 'c', 'g']
ccycle = itertools.cycle(colorcycle)
mpl.rcParams['axes.color_cycle'] = colorcycle

def __prelude__(args):

    if args.figsize:
        plt.figure(figsize=map(int, args.figsize.split(',')))

    plt.figure(frameon=False)

def __core__(args, bbox_inches='tight'):

    if args.hline:
        ax = plt.gca()
        ax.axhline(y=args.hline, ls='--')

    if args.vline:
        ax = plt.gca()
        ax.axvline(x=args.vline, ls='--', color='k')

    if args.xlabel:
        plt.xlabel(args.xlabel)
        
    if args.ylabel:
        plt.ylabel(args.ylabel)

    if args.xmin:
        ax = plt.gca()
        ax.set_xlim(left=args.xmin)
    
    if args.xmax:
        ax = plt.gca()
        ax.set_xlim(right=args.xmax)

    if args.ymin:
        ax = plt.gca()
        ax.set_ylim(bottom=args.ymin)

    if args.ymax:
        ax = plt.gca()
        ax.set_ylim(top=args.ymax)

    plt.savefig(args.outfig, bbox_inches=bbox_inches, dpi=args.dpi)

def main_scatter(args):

    __prelude__(args)

    if args.skipheader:
        args.table.readline()

    if args.skipline > 0:
        for i in xrange(args.skipline):
            args.table.readline()

    x = []
    y = []
    colors = []
    for i, line in enumerate(args.table):
        if i > args.maxline:
            break
        pair = line.strip().split(args.delim)
        x.append(float(pair[args.x-1]))
        y.append(float(pair[args.y-1]))
        if args.color>=0:
            colors.append(pair[args.color-1])

    if args.s:
        __ms = args.s
    else:
        __ms = args.ms

    if args.color>=0:
        cats = set(colors)
        for catind, cat in enumerate(cats):
            plt.scatter([x[i] for i, _ in enumerate(colors) if _ == cat], [y[i] for i, _ in enumerate(colors) if _ == cat], alpha=args.alpha, linewidth=0, s=__ms, color=next(ccycle))
    else:
        plt.scatter(x, y, alpha=args.alpha, linewidth=0, s=__ms, color = 'DarkSlateGray')
    if args.xlog:
        ax = plt.gca()
        ax.set_xscale("log")
    if args.ylog:
        ax = plt.gca()
        ax.set_yscale("log")


    __core__(args)


def main_hist(args):

    __prelude__(args)
    if args.skipheader:
        args.table.readline()

    if args.skipline > 0:
        for i in xrange(args.skipline):
            args.table.readline()

    data = []
    for i, line in enumerate(args.table):
        if i > args.maxline:
            break
        pair = line.strip().split(args.delim)
        data.append(float(pair[args.c-1]))

    if args.binseq:
        bins = map(int, args.binseq.split(','))
        hts, locs = np.histogram(data, bins=bins)
        coords = range(1,len(locs))
        plt.bar(coords, hts, edgecolor='none', linewidth=0, width=0.95,
                facecolor=next(ccycle), alpha=args.alpha)
        plt.xticks(coords, locs)
    elif args.xlog:
        bins=np.logspace(np.log10(min(data)), np.log10(max(data)), args.bins)
        hts, locs = np.histogram(data, bins=bins)
        coords = range(1,len(locs))
        plt.bar(coords, hts, edgecolor='none', linewidth=0, width=0.90,
                facecolor=next(ccycle), alpha=args.alpha)
        if len(coords)>8:
            inds = map(int, np.linspace(0, len(coords)-1, 8))
            inds = sorted(set(inds))
            plt.xticks([coords[_] for _ in inds], [int(locs[_]) for _ in inds])
        else:
            plt.xticks(coords, locs)
    else:
        plt.hist(data, bins=args.bins, log=args.log, linewidth=0, edgecolor='none', alpha=args.alpha, rwidth=0.9, facecolor=next(ccycle))
    __core__(args)

def main_box(args):

    def drawmean(starts, ends, means, color):
        for start, end, mean in zip(starts, ends, means):
            plt.hlines(y=mean, xmin=start, xmax=end, color=color)

        return

    def setcolor(ob, col):
        for item in ['medians', 'fliers', 'whiskers', 'boxes', 'caps']:
            plt.setp(ob[item], color=col, lw=1.3)
            
        plt.setp(ob["medians"], lw=0)

        return

    __prelude__(args)

    if args.dim1>=0 and args.dim2>=0:
        if args.skipheader:
            args.table.readline()

        if args.skipline > 0:
            for i in xrange(args.skipline):
                args.table.readline()

        cat2key2vals = {}
        for i, line in enumerate(args.table):
            if i > args.maxline:
                break
            pair = line.strip().split(args.delim)
            cat = pair[args.dim2-1]
            key = pair[args.dim1-1]
            if cat not in cat2key2vals:
                cat2key2vals[cat] = {}
            if key not in cat2key2vals[cat]:
                cat2key2vals[cat][key] = []
            cat2key2vals[cat][key].append(float(pair[args.c-1]))

        keys = set([val for _ in cat2key2vals.values() for val in _.keys()])
        keys = sorted(list(keys), key=lambda x: int(x))

        cats = cat2key2vals.keys()
        xmin = 1
        curr = xmin
        base_poses = []
        for p in range(1, len(keys)+1):
            base_poses.append(curr)
            curr += len(cats) * (args.width+0.1) + args.catsep

        for i, cat in enumerate(cats):
            key2vals = cat2key2vals[cat]
            ob = plt.boxplot([key2vals[_] if _ in key2vals else [] for _ in keys], positions=[_+i*(args.width+0.1) for _ in base_poses], widths=args.width)
            setcolor(ob, next(ccycle))
            
        ax=plt.gca()
        ax.set_xlim(left=xmin-0.5-args.width/2)
        ax.set_xticks(base_poses)
        ax.set_xticklabels(keys, size=args.xtlsize, rotation=args.xtlrotat)
            
        __core__(args)

def main_bar(args):

    __prelude__(args)

    rows = []
    catlabels = args.table.readline().strip().split('\t')
    print catlabels
    for line in args.table:
        pair = line.strip().split('\t')
        rows.append([pair[0]]+map(float, pair[1:]))

    cols = zip(*rows)

    ncats = len(cols)-1
    xdim = len(cols[0])
    sep = 0.2
    lefts = np.arange(1, (ncats+sep)*xdim+1, ncats+sep)
    mids = np.arange(1+(ncats+sep)/2, (ncats+sep)*xdim+1, ncats+sep)
    for col in cols[1:]:
        plt.bar(lefts, [_+0.1 for _ in col], width=0.98, facecolor=next(ccycle), linewidth=0, alpha=args.alpha)
        lefts += 1

    plt.axis('equal')
    plt.legend(catlabels, loc=args.legloc)
    plt.xticks(mids, cols[0])
    plt.ylim(bottom = -1)
    plt.xlim(left=0.5, right=(ncats+sep)*xdim+1)

    # # stacked
    # bottom = np.zeros(len(cols[0]))
    # print bottom
    # for col in cols[1:]:
    #     print col, bottom
    #     plt.bar(range(1,len(col)+1), height=col, bottom = bottom)
    #     bottom += col

    __core__(args)
        
def main_pie(args):

    __prelude__(args)

    fracs = []
    labels = []
    for line in args.table:
        pair = line.strip().split('\t')
        fracs.append(float(pair[args.frac-1]))
        labels.append(pair[args.label-1])

    (fracs,labels) = zip(*sorted(zip(fracs,labels)))
    # cmap = plt.cm.prism
    # colors = cmap(np.linspace(0., 1., len(fracs)))
    pie_wedges = plt.pie(fracs, labeldistance=1.05,
                         autopct='%1.1f%%', colors=colorcycle)
    for pie_wedge in pie_wedges[0]:
        pie_wedge.set_edgecolor('white')
        pie_wedge.set_alpha(args.alpha)

    plt.axis('equal')
    plt.legend(labels, loc=args.legloc, shadow=True, prop={'size':11})

    __core__(args)

def add_std_options(psr):
    psr.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    psr.add_argument('--delim', default="\t", 
                     help="table delimiter [\\t]")
    psr.add_argument('-o', dest='outfig', default="tmp.png",
                     help='output figure file name [tmp.png]')
    psr.add_argument('--xlabel', help='x label')
    psr.add_argument('--ylabel', help='y label')
    psr.add_argument('--skipheader', action='store_true', help='skip header')
    psr.add_argument('--skipline', default=0, type=int, help='skip line from table')
    psr.add_argument('--maxline', default=100000, type=int, help='number of lines to plot [100000]')

    psr.add_argument('--vline', default=None, type=float, help="draw an extra vertical line at given x position")
    psr.add_argument('--hline', default=None, type=float, help="draw an extra horizontal line at given y position")
    psr.add_argument('--xmin', default=None, type=float, help="minimum x [auto]")
    psr.add_argument('--xmax', default=None, type=float, help="maximum x [auto]")
    psr.add_argument('--ymin', default=None, type=float, help="minimum y [auto]")
    psr.add_argument('--ymax', default=None, type=float, help="maximum y [auto]")
    psr.add_argument('--alpha', default=0.9, type=float, help="alpha")
    psr.add_argument('--dpi', default=100, type=int, help="resolution in dpi")
    psr.add_argument('--xtlsize', default='large', help="size of x tick label")
    psr.add_argument('--xtlrotat', default='horizontal', help="rotation of x tick label. {horizontal, vertical, angle in degrees} [horizontal]")
    psr.add_argument('--figsize', default=None, help="figure size in tuple. e.g., 10,8 ([width],[height]).")
    psr.add_argument('--legloc', default=2, type=int, help='legend location')
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='plotter')
    subparsers = parser.add_subparsers()

    # x,y scatter plot
    psr_scatter = subparsers.add_parser("scatter", help=""" x-y scatter plot """)
    psr_scatter.add_argument('-x', type=int, help="column of x (1-based)")
    psr_scatter.add_argument('-y', type=int, help="column of y (1-based)")
    psr_scatter.add_argument('--color', type=int, default=-1, help="column of color (1-based)")
    psr_scatter.add_argument('-s', type=int, default=None, help="column to plot s (1-based, size of the markers)")
    psr_scatter.add_argument('--ms', type=float, default=10, help="fixed marker size [10]")
    psr_scatter.add_argument('--xlog', action="store_true", help="x axis log-scaled")
    psr_scatter.add_argument('--ylog', action="store_true", help="y axis log-scaled")
    add_std_options(psr_scatter)
    psr_scatter.set_defaults(func=main_scatter)

    # boxplot
    psr_box = subparsers.add_parser("box", help=""" boxplot """)
    psr_box.add_argument('-c', type=int, help="column of value to plot [required]")
    psr_box.add_argument('--dim1', type=int, default=-1, help="column of dimention 1 category (1-based)")
    psr_box.add_argument('--dim2', type=int, default=-1, help="column of dimention 2 category (1-based)")
    psr_box.add_argument('--width', type=int, default=0.5, help="width of each boxplot [0.5]")
    psr_box.add_argument('--catsep', type=int, default=0.5, help="separation between categories [0.1]")
    add_std_options(psr_box)
    psr_box.set_defaults(func=main_box)

    # histogram
    psr_hist = subparsers.add_parser("hist", help=""" histogram """)
    psr_hist.add_argument('-c', type=int, help="column to plot (1-based)")
    psr_hist.add_argument('--bins', type=int,
                          default=20, help="number of bins [20]")
    psr_hist.add_argument('--binseq', default=None, help='bin sequences. e.g., 1,20,100,400,1000 means [1,20),[20,100)...')
    psr_hist.add_argument('--log', action='store_true', help='log scaled x and y')
    psr_hist.add_argument('--xlog', action='store_true', help='log scaled x only')
    add_std_options(psr_hist)
    psr_hist.set_defaults(func=main_hist)

    # pie
    psr_pie = subparsers.add_parser("pie", help=""" pie chart """)
    psr_pie.add_argument('-c', type=int, help="column to plot (1-based)")
    psr_pie.add_argument('--frac', type=int, default=1, help='index of fraction column')
    psr_pie.add_argument('--label', type=int, default=2, help='index of label column')
    add_std_options(psr_pie)
    psr_pie.set_defaults(func=main_pie)

    # bar
    psr_bar = subparsers.add_parser("bar", help=""" bar chart """)
    psr_bar.add_argument('-l', type=int, help="x dimensional column index (1-based)")
    psr_bar.add_argument('--cat', help='column for category (stacked or juxtaposed)')
    psr_bar.add_argument('--catlabels', help='category labels')
    # psr_bar.add_argument('--frac', type=int, default=1, help='index of fraction column')
    # psr_bar.add_argument('--label', type=int, default=2, help='index of label column')
    add_std_options(psr_bar)
    psr_bar.set_defaults(func=main_bar)


    args = parser.parse_args()
    args.func(args)
