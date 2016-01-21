#!/usr/bin/env python

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import argparse
import numpy as np
import matplotlib.cm as cm
import itertools
import numpy as np
import operator
import re

class Indices:

    def __init__(self):
        self.spans = []

    def extend(self, start, end):
        self.spans.append((start, end))

    def extract(self, lst):
        result = []
        for start, end in self.spans:
            if not end:
                end = len(lst)
            result.extend([lst[_] for _ in xrange(start, end)])

        return result

def parse_indices(indstr):
    indices = Indices()
    if not indstr: return indices
    rgs = indstr.split(',')
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = None
            indices.extend(int(pair[0])-1 if pair[0] else 0,
                           int(pair[1]) if pair[1] else None)
        else:
            indices.extend(int(rg)-1, int(rg))

    return indices


colorcycle = ['k', 'DarkSlateGray', 'Brown', 'Burlywood', 'DarkSlateBlue', 'yellowgreen', 'gold', 'lightskyblue', 'lightcoral', 'y', 'k', 'c', 'g']
ccycle = itertools.cycle(colorcycle)
mpl.rcParams['axes.color_cycle'] = colorcycle

def __prelude__(args):

    if args.figsize:
        plt.figure(figsize=tuple(map(int, args.figsize.split(','))), frameon=False)
    else:
        plt.figure(frameon=False)

def __core__(args, bbox_inches='tight'):

    if args.hline:
        ax = plt.gca()
        for hval in args.hline:
            ax.axhline(y=hval, ls='--', color='k')

    if args.vline:
        ax = plt.gca()
        for vval in args.vline:
            ax.axvline(x=vval, ls='--', color='k')

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

    if args.outfig:
        plt.savefig(args.outfig, bbox_inches=bbox_inches, dpi=args.dpi)
    else:
        plt.show()

def main_scatter(args):

    __prelude__(args)

    if args.skipheader:
        headerfields = args.table.readline().strip().split(args.delim)
        if not args.xlabel:
            args.xlabel = headerfields[args.x-1]
        if not args.ylabel:
            args.ylabel = headerfields[args.y-1]

    if args.skipline > 0:
        for i in xrange(args.skipline):
            args.table.readline()

    x = []
    y = []
    colors = []
    sizes = []
    for i, line in enumerate(args.table):
        if i > args.maxline:
            break
        field = line.strip().split(args.delim)
        x.append(float(field[args.x-1]))
        y.append(float(field[args.y-1]))
        if args.color>=0:
            colors.append(field[args.color-1])
        if args.s:
            sizes.append(float(field[args.s-1]))

    if sizes:                   # scale size
        minsize = args.minsize; maxsize = args.maxsize
        sizerange = maxsize - minsize
        _minsize = min(sizes); _maxsize = max(sizes);
        _sizerange = _maxsize - _minsize
        sizes = [minsize + float(_ - _minsize) / _sizerange * sizerange for _ in sizes]

    if args.s:
        __ms = sizes
    else:
        __ms = args.ms

    if args.beta != 0:
        minx = min(x)
        maxx = max(x)
        plt.plot([minx, maxx],
                 [args.gamma,
                  args.gamma+args.beta*(maxx-minx)],
                 linestyle='dashed', color='k')

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

    if args.color >= 0:
        plt.legend(cats, loc=args.legloc)

    __core__(args)


def main_ezscatter_table(args):

    """
    plot x vs multiple y with each y in a separate column (table format)
    assumes the reader line contains the label for each category
    e.g.,
    x y_apple y_orange y_banana
    1 1       1        1
    2 2       2        2
    ...
    """

    __prelude__(args)

    yinds = parse_indices(args.ys)
    labels = args.table.readline().strip().split('\t')
    # need to make xs because there might exist 'NA' for some y
    xs = [[] for _ in yinds.extract(labels)]
    ys = [[] for _ in xs]
    for i, line in enumerate(args.table):
        field = line.strip().split(args.delim)
        for j, yval in enumerate(yinds.extract(field)):
            if yval != 'NA':
                xs[j].append(float(field[args.x-1]))
                ys[j].append(float(yval))

    for j in xrange(len(xs)):
        plt.scatter(xs[j], ys[j], color=colorcycle[j])

    plt.legend(yinds.extract(labels))

    __core__(args)

                
def main_ezscatter_factor(args):

    """
    plot x vs multiple y in with a factor column
    e.g.,
    x  y  category
    1  1  apple
    2  2  orange
    ...
    """
    pass


def adjust_spines(ax,spines):
    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward',10)) # outward by 10 points
            spine.set_smart_bounds(True)
        else:
            spine.set_color('none') # don't draw spine

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

def _hist_multiread(args):
    

    cat2data = {}
    for i, line in enumerate(args.table):
        fields = line.strip().split(args.delim)
        cat = fields[args.cat-1]
        datastr = fields[args.c-1]
        if datastr == 'NA':
            continue
        data = float(datastr)

        # find min, max to determine range
        if i==0:
            dmin = data
            dmax = data
        if dmin > data:
            dmin = data
        if dmax < data:
            dmax = data

        if cat in cat2data:
            cat2data[cat].append(data)
        else:
            cat2data[cat] = [data]

    return cat2data, dmin, dmax


def _hist_multiplot(args):

    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'

    cat2data, dmin, dmax = _hist_multiread(args)

    if args.xrange:
        dmin, dmax = eval(args.xrange)

    fig = plt.figure()
    ncats = len(cat2data)
    for i, (cat, data_seq) in enumerate(cat2data.iteritems()):
        ax = fig.add_subplot(ncats, 1, i+1)
        color = next(ccycle)
        ax.hist(data_seq, range=(dmin, 45), alpha=args.alpha, ec=color, fc=color)
        adjust_spines(ax, ["left", "bottom"])
        ax.set_ylabel(str(cat))

    ax=plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    return

def main_ezhist_stack(args):

    """
    this plots histogram in which different categories stack on top of each other
    e.g.,
    hist1
    hist2
    hist3
    """

    pass
    

def main_ezhist_sidebyside(args):

    """
    this plots histograms in which different categories stand side by side
    input format:
    category  count
    cat1      12
    cat1      11
    ...
    cat2      13
    cat3      18
    
    plots
    cat1,cat2,cat3  cat1,cat2,cat3
         bin1            bin2
    """

    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    
    __prelude__(args)
    cat2data, dmin, dmax = _hist_multiread(args)
    ncats = len(cat2data)

    barwidth = 0.8
    spacing = 0.3
    binwidth = barwidth * ncats + spacing

    cats = []
    obs  = []
    for i, (cat,data_seq) in enumerate(cat2data.iteritems()):
        hts, edges = np.histogram(data_seq, range=(dmin, dmax), bins=args.bins)
        binlefts = [1+j*binwidth for j in xrange(len(hts))]
        barlefts = [k+i*barwidth for k in binlefts]
        ob = plt.bar(barlefts, hts, edgecolor="none", linewidth=0, width=barwidth, facecolor=next(ccycle), alpha=args.alpha)
        obs.append(ob)
        cats.append(cat)
    plt.xticks(binlefts, map(lambda a: '%1.2g' % a, edges[:-1]), rotation=args.xtlrotat)

    plt.legend(obs, cats)

    __core__(args)

def main_ezhist_substack(args):

    """
    this plots histogram in different subplots stacking on top of each other
    """
    pass



def _hist_oneplot(args):

    if args.skipheader:
        headerline = args.table.readline()

    if args.skipline > 0:
        for i in xrange(args.skipline):
            args.table.readline()

    data = []
    dmax = None
    dmin = None
    for i, line in enumerate(args.table):
        if i > args.maxline:
            break
        pair = line.strip().split(args.delim)
        if pair[args.c-1] == 'NA':
            continue
        d = float(pair[args.c-1])

        if args.xmax and d > args.xmax:
            continue

        if args.xmin and d < args.xmin:
            continue
        
        data.append(d)
        
        if dmax is None or dmax < d:
            dmax = d
        if dmin is None or dmin > d:
            dmin = d

    if args.binseq:
        bins = map(int, args.binseq.split(','))
        hts, locs = np.histogram(data, bins=bins)
        coords = range(1,len(locs))
        plt.bar(coords, hts, edgecolor='none', linewidth=0, width=0.95,
                facecolor=next(ccycle), alpha=args.alpha)
        plt.xticks(coords, locs)
    elif args.xlog:
        bins=np.logspace(np.log10(max(0.1, min(data))), np.log10(max(data)), args.bins)
        hts, locs = np.histogram(data, bins=bins)
        coords = range(1,len(locs))
        plt.bar(coords, hts, edgecolor='none', linewidth=0, width=0.90,
                facecolor=next(ccycle), alpha=args.alpha)
        if len(coords)>8:
            inds = map(int, np.linspace(0, len(coords)-1, 8))
            inds = sorted(set(inds))
            plt.xticks([coords[_] for _ in inds], [int(locs[_]) for _ in inds], rotation=args.xtlrotat)
        else:
            plt.xticks(coords, locs, rotation=args.xtlrotat)
    else:
        plt.hist(data, bins=args.bins, log=args.log, linewidth=0, edgecolor='none', alpha=args.alpha, rwidth=0.9, facecolor=next(ccycle))
        plt.xticks(rotation=args.xtlrotat)

    if not args.xlabel and args.skipheader:
        args.xlabel = headerline.strip().split(args.delim)[args.c-1]

    return

def main_hist(args):

    __prelude__(args)

    if args.cat >= 0:

        _hist_multiplot(args)

    else:

        _hist_oneplot(args)

    __core__(args)

def main_cumhist(args):
    
    """ cumulative stepwise histogram with increment 1 """
    data = []
    dmax = None
    dmin = None
    for i, line in enumerate(args.table):
        pair = line.strip().split(args.delim)
        if pair[args.c-1] == 'NA':
            continue
        d = float(pair[args.c-1])

        if args.xmax and d > args.xmax:
            continue

        if args.xmin and d < args.xmin:
            continue
        
        data.append(d)
        
        if dmax is None or dmax < d:
            dmax = d
        if dmin is None or dmin > d:
            dmin = d

    plt.hist(data, bins=range(int(np.percentile(data, 99))+1), histtype='step', cumulative=-1)

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

    if args.dim1>=0 and args.dim2>=0: # two levels of category
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

    elif args.dim1 >= 0:        # one level of category

        k2vs = {}
        if args.inputraw:
            for i, line in enumerate(args.table):
                if i > args.maxline:
                    break
                fields = line.strip().split(args.delim)
                k = fields[args.dim1-1]
                vs = eval(fields[args.c-1])
                if k not in vs:
                    k2vs[k] = []
                k2vs[k].extend(vs)

        ks, vss = zip(*k2vs.items())

        lefts = range(1,len(ks)+1)
        mids = [_+args.width/2 for _ in lefts]
        plt.boxplot(vss, positions=lefts, widths=args.width)
        ax=plt.gca()
        ax.set_xlim(left=lefts[0]-0.5-args.width/2)
        ax.set_xticks(mids)
        ax.set_xticklabels(ks, size=args.xtlsize, rotation=args.xtlrotat)

        __core__(args)



def main_boxbin(args):

    """ x is continuous and binned 
    y is plotted
    """

    __prelude__(args)

    def drawmean(starts, ends, means, color):
        for start, end, mean in zip(starts, ends, means):
            plt.hlines(y=mean, xmin=start, xmax=end, color=color, lw=3)

        return

    def setcolor(ob, col):
        for item in ['medians', 'fliers', 'whiskers', 'boxes', 'caps']:
            plt.setp(ob[item], color=col, lw=1.3)
            
        plt.setp(ob["medians"], lw=1.3)

        return
    
    yinds = parse_indices(args.ys)
    labels = args.table.readline().strip().split(args.delim)

    ylabels = yinds.extract(labels)
    xs = [[] for _ in ylabels]
    ys = [[] for _ in xs]
    x0 = []
    for i, line in enumerate(args.table):
        fields = line.strip().split(args.delim)

        if fields[args.x-1] == 'NA':
            continue

        x = float(fields[args.x-1])
        x0.append(x)
        for j, yval in enumerate(yinds.extract(fields)):
            if yval == 'NA':
                continue
            y = float(yval)
            xs[j].append(x)
            ys[j].append(float(yval))

    grids = np.linspace(min(x0), max(x0), args.bins)
    obs = []
    for j, y in enumerate(ys):
        dataseq = [[y[i] for i, x in enumerate(xs[j]) if x >= grids[k] and x < grids[k+1]] for k in xrange(len(grids)-1)]
        centerpoints = range(1, len(dataseq)+1)
        ob = plt.boxplot(dataseq, sym=args.symb, positions=centerpoints, widths=0.2)
        if args.mean:
            drawmean([_-0.1 for _ in centerpoints], [_+0.1 for _ in centerpoints], map(np.mean, dataseq), colorcycle[j])
        setcolor(ob, colorcycle[j])
        obs.append(ob['boxes'][0])

    if not args.xlabel:
        args.xlabel = labels[args.x-1]

    if not args.ylabel and len(ys) == 1:
        args.ylabel = yinds.extract(labels)[0]

    plt.xticks(range(len(grids)), map(lambda a: '%1.2f' % a, grids))
        
    if len(ys) > 1:
        plt.legend(obs, yinds.extract(labels), loc=args.legloc)
    
    __core__(args)

def main_bar(args):

    """
    simple bar plot

    input:
    xlabel   ylabel
    bar1name bar1height
    bar2name bar2height
    bar3name bar3height
    """

    __prelude__(args)

    fields = args.table.readline().strip().split('\t')
    if not args.xlabel:
        args.xlabel = fields[args.x-1]
    if not args.ylabel:
        args.ylabel = fields[args.y-1]
    
    x = []
    y = []
    if (re.match(r'^[\-\.0-9]*$', fields[args.y-1])):
        y.append(float(fields[args.y-1]))
        x.append(fields[args.x-1])
        if not args.ylabel:
            args.ylabel = "column %d" % args.y
        if not args.xlabel:
            args.xlabel = "column %d" % args.x
    else:
        if not args.ylabel:
            args.ylabel = "column %d" % args.y
        if not args.xlabel:
            args.xlabel = "column %d" % args.x

    for line in args.table:
        fields = line.strip().split('\t')
        y.append(float(fields[args.y-1]))
        x.append(fields[args.x-1])

    # print x
    # print y

    barwidth = 1.0
    barsep = 0.2
    barstep = barwidth + barsep
    
    lefts = [1+_*barstep for _ in xrange(len(x))]
    mids = [_+barwidth/2 for _ in lefts]
    # print len(x)
    # print lefts
    # print mids
    plt.bar(lefts, y, width=barwidth, facecolor=next(ccycle), linewidth=0, alpha=args.alpha)

    if not args.xmin:
        args.xmin = 1-barwidth
    if not args.xmax:
        args.xmax = lefts[-1]+2*barwidth
    # rows = []
    # catlabels = args.table.readline().strip().split('\t')
    # print catlabels
    # for line in args.table:
    #     pair = line.strip().split('\t')
    #     rows.append([pair[0]]+map(float, pair[1:]))

    # cols = zip(*rows)
    # ncats = len(cols)-1
    # xdim = len(cols[0])
    # sep = 0.2
    # lefts = np.arange(1, (ncats+sep)*xdim+1, ncats+sep)
    # mids = np.arange(1+(ncats+sep)/2, (ncats+sep)*xdim+1, ncats+sep)
    # for col in cols[1:]:
    #     plt.bar(lefts, [_+0.1 for _ in col], width=0.98, facecolor=next(ccycle), linewidth=0, alpha=args.alpha)
    #     lefts += 1

    # plt.axis('equal')
    # plt.legend(catlabels, loc=args.legloc)
    plt.xticks(mids, x, rotation=args.xtlrotat)
    # plt.ylim(bottom = -1)
    # plt.xlim(left=0.5, right=(ncats+sep)*xdim+1)

    # # stacked
    # bottom = np.zeros(len(cols[0]))
    # print bottom
    # for col in cols[1:]:
    #     print col, bottom
    #     plt.bar(range(1,len(col)+1), height=col, bottom = bottom)
    #     bottom += col

    __core__(args)


def _bar_multiread_table(args):

    indices = parse_indices(args.hs)
    headerline = args.table.readline()
    headerfields = headerline.strip().split(args.delim)
    xlabel = headerfields[args.x-1]
    cats = indices.extract(headerfields)
    
    data = [[] for j in xrange(len(cats))]
    xticklabels = []
    for line in args.table:
        fields = line.strip().split(args.delim)
        xticklabels.append(fields[args.x-1])
        for j, n in enumerate(indices.extract(fields)):
            data[j].append(float(n))

    return xlabel, xticklabels, cats, data

def main_ezbar_sidebyside(args):

    """ input format:

    header cat1 cat2
    x1     n11   n12
    x2     n21   n22
    x3     n31   n32
    ...

    results:


    bar(n11) bar(n12)   bar(n21) bar(n22)   bar(n31) bar(n32)
    ==========================================================
           x1                   x2                 x3

    """
    
    xlabels, xticklabels, cats, data = _bar_multiread_table(args)

    __prelude__(args)
    ncats = len(cats)
    barwidth = args.barwidth
    spacing = args.spacing
    binwidth = barwidth * ncats + spacing

    obs = []
    for i, bardata in enumerate(data):
        binlefts = [1+j*binwidth for j in xrange(len(bardata))]
        barlefts = [k+i*barwidth for k in binlefts]
        ob = plt.bar(barlefts, bardata, edgecolor="none", linewidth=0, width=barwidth, facecolor=next(ccycle), alpha=args.alpha)
        obs.append(ob)

    binrights = [k+barwidth*ncats+spacing for k in binlefts]
    bincenter = map(lambda x,y: (x+y)/2.0, binlefts, binrights)
    plt.xticks(bincenter, xticklabels, rotation=args.xtlrotat, horizontalalignment='center')
    plt.legend(obs, cats)

    if not args.xmin:
        args.xmin = binlefts[0]-binwidth
        
    if not args.xmax:
        args.xmax = binrights[-1]+binwidth

    __core__(args)
        
    

def main_ezbar_stack(args):


    """ input format:
    
    header cat1 cat2
    x1     n11  n12
    x2     n21  n22
    x3     n31  n32
    ...

    bar(n12)    bar(n22)   bar(n32)
    bar(n11)    bar(n21)   bar(n31)
    ===============================
       x1          x2         x3
    """
    xlabels, xticklabels, cats, data = _bar_multiread_table(args)

    ncats = len(cats)
    barwidth = 2
    spacing = 0.4
    binwidth = barwidth + spacing

    obs = []
    for i, bardata in enumerate(data):
        barlefts = [1+j*binwidth for j in xrange(len(bardata))]
        if i == 0:
            bottom = [0 for _ in bardata]
        ob = plt.bar(barlefts, bardata, bottom=bottom, edgecolor="none", linewidth=0, width=barwidth, facecolor=next(ccycle), alpha=args.alpha)
        obs.append(ob)
        bottom = map(operator.add, bottom, bardata)
        bottom = map(lambda x: x+args.vspacing, bottom)

    plt.xticks(barlefts, xticklabels, rotation=args.xtlrotat)
    plt.legend(obs, cats)

    __core__(args)

def main_ezbar_substack(args):

    """ same as stack but in different subplots so to allow different scale """
    xlabels, xticklabels, cats, data = _bar_multiread_table(args)

    ncats = len(cats)
    barwidth = 2
    spacing = 0.4
    binwidth = barwidth + spacing

    obs = []
    nrows = len(data)
    plt.subplots_adjust(hspace=0)
    for i, bardata in enumerate(data):
        plt.subplot(nrows, 1, i+1)
        barlefts = [1+j*binwidth for j in xrange(len(bardata))]
        ob = plt.bar(barlefts, bardata, edgecolor="none", linewidth=0, width=barwidth, facecolor=next(ccycle), alpha=args.alpha)
        obs.append(ob)
        plt.ylabel(cats[i])
        if i != 0:
            yticks, yticklabels = plt.yticks()
            plt.yticks(yticks[:-2])
        if i != nrows - 1:
            plt.xticks([])

        barmids = [_+barwidth/2 for _ in barlefts]
        plt.xticks(barmids, xticklabels, horizontalalignment='center', rotation=args.xtlrotat)
        
        if not args.xmin:
            plt.gca().set_xlim(left=1-barwidth)
        if not args.xmax:
            plt.gca().set_xlim(right=barlefts[-1]+barwidth*2)

    __core__(args)


### horizontal bar
# #!/usr/bin/env python
# # -*- coding: utf-8 -*-
# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches

# # ----------
# # Data to be represented

# products = ['Vendor A - Product A', 'Vendor A - Product B', 'Vendor A - Product C',
#             'Vendor B - Product A', 'Vendor B - Product B', 'Vendor B - Product C',
#             'Vendor C - Product A', 'Vendor C - Product B', 'Vendor C - Product C']

# values = np.random.uniform(10,60,len(products))

# # ----------

# # Choose some nice colors
# matplotlib.rc('axes', facecolor = '#6E838A')
# matplotlib.rc('axes', edgecolor = '#737373')
# matplotlib.rc('axes', linewidth = 1)
# matplotlib.rc('ytick', direction='out')
# matplotlib.rc('xtick', direction='out')
# matplotlib.rc('figure.subplot', left=0.25)

# # Make figure background the same colors as axes
# fig = plt.figure(figsize=(8,6), facecolor='#6E838A')

# # Remove left and top axes spines
# axes = plt.subplot(1,1,1)
# axes.spines['right'].set_color('none')
# axes.spines['top'].set_color('none')
# axes.xaxis.set_ticks_position('bottom')
# axes.yaxis.set_ticks_position('left')

# # Adjust yticks to the number of products
# plt.yticks(np.arange(len(products)+1), [])

# # Set tick labels color to white
# for label in axes.get_xticklabels()+axes.get_yticklabels():
#     label.set_color('white')

# # Set tick labels line width to 1
# for line in axes.get_xticklines() + axes.get_yticklines():
#     line.set_markeredgewidth(1)

# # Set axes limits
# ymin, ymax = 0, len(products)
# xmin, xmax = 0, 60
# plt.xlim(xmin,xmax)
# plt.ylim(ymin,ymax)

# # Start with blue colormap
# cmap = plt.cm.Blues

# for i, label in enumerate(products):

#     # Alternate band of light background
#     if not i%2:
#         p = patches.Rectangle(
#             (0, i), xmax, 1, fill=True, transform=axes.transData,
#             lw=0, facecolor='w', alpha=.1)
#         axes.add_patch(p)

#     # Product name left to the axes
#     plt.text(-.5, i+0.5, label, color="white", size=10,
#              horizontalalignment='right', verticalalignment='center')

#     # Plot the bar with gradient (1 to .65)
#     value = values[i]
#     X = np.array([1,.65]).reshape((1,2))
#     axes.imshow(X,extent=(0,value,i+.25,i+.75),cmap=cmap, vmin=0, vmax=1)
#     plt.text(value-0.5, i+0.5, '%.1f' % value, color="white", size=10,
#              horizontalalignment='right', verticalalignment='center')

#     # Change colormap every 3 values
#     if i >= 2: cmap = plt.cm.Greens
#     if i >= 5: cmap = plt.cm.Reds

# # Set a nice figure aspect
# axes.set_aspect(4.5)

# # Write some title & subtitle
# plt.text(1, 10.0, "Vendor benchmarks", color="1.0", fontsize=14)
# plt.text(1,  9.7, "(higher is better)", color="0.75", fontsize=10)

# # Done
# matplotlib.rc('savefig', facecolor = '#6E838A')
# plt.show()


# polar plot

# # Data to be represented
# # ----------
# labels = ['January', 'Feburary', 'March', 'April', 'May', 'June',
#           'July', 'August', 'September', 'October', 'November', 'December']
# n = len(labels)
# data = np.random.uniform(0,1,n)
# # ----------

# # Make figure square and background the same colors as axes (white)
# fig = plt.figure(figsize=(8,6), facecolor='white')

# # Make a new polar axis
# axes = plt.subplot(111, polar=True, axisbelow=True)

# # Put labels on outer 
# T = np.arange(np.pi/n, 2*np.pi, 2*np.pi/n)
# R = np.ones(n)*10
# width = 2*np.pi/n

# # Label background
# bars  = axes.bar(T, R, width=width, bottom=9,
#                  linewidth = 2, facecolor = '0.9', edgecolor='1.00')
# # Labels
# for i in range(T.size):
#     theta = T[n-1-i]+np.pi/n + np.pi/2
#     plt.text(theta, 9.5, labels[i], rotation=180*theta/np.pi-90,
#              family='Helvetica Neue', size=7,
#              horizontalalignment="center", verticalalignment="center")

# # Data
# R = 1 + data*6
# bars = axes.bar(T, R, width=width, bottom=2,
#                 linewidth=1, facecolor = '0.75', edgecolor='1.00')
# for i,bar in enumerate(bars):
#     bar.set_facecolor(plt.cm.hot(R[i]/10))

# # Text i the center
# plt.text(1*np.pi/2, 0.05, "2012",
#          size=16, family='Helvetica Neue Light',
#          horizontalalignment="center", verticalalignment="bottom")
# plt.text(3*np.pi/2, 0.05, "some levels", color="0.50",
#          size=8, family='Helvetica Neue Light',
#          horizontalalignment="center", verticalalignment="top")

# # Set ticks, tick labels and grid
# plt.ylim(0,10)
# plt.xticks(T)
# plt.yticks(np.arange(2,9))
# axes.grid(which='major', axis='y', linestyle='-', color='0.75')
# axes.grid(which='major', axis='x', linestyle='-', color='1.00')
# for theta in T:
#     axes.plot([theta,theta], [4,9], color='w', zorder=2, lw=1)
# axes.set_xticklabels([])
# axes.set_yticklabels([])

# plt.show()

# cluster and dendrogram

# import scipy
# import pylab
# import scipy.cluster.hierarchy as sch

# # Generate features and distance matrix.
# x = scipy.rand(40)
# D = scipy.zeros([40,40])
# for i in range(40):
#     for j in range(40):
#         D[i,j] = abs(x[i] - x[j])

# # Compute and plot dendrogram.
# fig = pylab.figure()
# axdendro = fig.add_axes([0.09,0.1,0.2,0.8])
# Y = sch.linkage(D, method='centroid')
# Z = sch.dendrogram(Y, orientation='right')
# axdendro.set_xticks([])
# axdendro.set_yticks([])

# # Plot distance matrix.
# axmatrix = fig.add_axes([0.3,0.1,0.6,0.8])
# index = Z['leaves']
# D = D[index,:]
# D = D[:,index]
# im = axmatrix.matshow(D, aspect='auto', origin='lower')
# axmatrix.set_xticks([])
# axmatrix.set_yticks([])

# # Plot colorbar.
# axcolor = fig.add_axes([0.91,0.1,0.02,0.8])
# pylab.colorbar(im, cax=axcolor)

# # Display and save figure.
# fig.show()
# fig.savefig('dendrogram.png')

def main_pie(args):

    __prelude__(args)

    if args.skipheader:
        args.table.readline()

    fracs = []
    labels = []
    for line in args.table:
        pair = line.strip().split('\t')
        fracs.append(float(pair[args.frac-1]))
        labels.append(pair[args.label-1])

    (fracs,labels) = zip(*sorted(zip(fracs,labels), reverse=True))
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


def main_venn(args):

    __prelude__(args)

    cnts = eval(args.c)
    labels = args.l.split(',')

    if args.n == 2:
        from matplotlib_venn import venn2
        venn2(subsets=cnts, set_labels=labels)
    elif args.n == 3:
        from matplotlib_venn import venn3
        venn3(subsets=cnts, set_labels=labels)

    __core__(args)


def main_hexbin(args):

    __prelude__(args)

    x = []
    y = []
    for line in args.table:
        fields = line.strip().split(args.delim)
        x.append(float(fields[args.x-1]))
        y.append(float(fields[args.y-1]))

    bins = (args.bins, args.ybins) if args.ybins>0 else args.bins
    if args.nolog:
        plt.hexbin(x,y, gridsize=bins, cmap=plt.cm.YlGn)
    else:
        plt.hexbin(x,y, gridsize=bins, bins='log', cmap=plt.cm.YlGn)
    if args.scatter:
        plt.scatter(x,y,s=0.5,alpha=0.5,edgecolor='none')

    # plt.xlim(0,20)
    __core__(args)

def main_dhist(args):

    from collections import Counter
    cs = []
    for line in args.i:
        fields = line.strip('\n').split('\t')
        cs.append(fields[args.c-1])

    cter = Counter(cs)
    items = sorted(cter.iteritems(),key=lambda k:k[1],reverse=True)
    names, vals = zip(*items)
    ps = range(len(names))
    plt.bar(ps, vals)
    plt.xticks([_+0.5 for _ in ps], names, rotation=90)
    __core__(args)


def add_std_options_table(psr):

    psr.add_argument('-t', '--table', help="data table", type = argparse.FileType('r'), default='-')
    psr.add_argument('--delim', default="\t", 
                     help="table delimiter [\\t]")
    psr.add_argument('--skipheader', action='store_true', help='skip header')
    psr.add_argument('--skipline', default=0, type=int, help='skip line from table')
    psr.add_argument('--maxline', default=1000000, type=int, help='number of lines to plot [1000000]')


def add_std_options_out(psr):

    psr.add_argument('-o', dest='outfig', default=None,
                     help='output figure file name [tmp.png]')
    psr.add_argument('--xlabel', help='x label')
    psr.add_argument('--ylabel', help='y label')

    psr.add_argument('--vline', action="append", type=float, help="draw an extra vertical line at given x position")
    psr.add_argument('--hline', action="append", type=float, help="draw an extra horizontal line at given y position")
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

def add_std_options(psr):

    add_std_options_table(psr)
    add_std_options_out(psr)
    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='plotter')
    subparsers = parser.add_subparsers()

    # x,y scatter plot
    psr_scatter = subparsers.add_parser("scatter", help=""" x-y scatter plot """)
    psr_scatter.add_argument('-x', type=int, help="column of x (1-based)")
    psr_scatter.add_argument('-y', type=int, help="column of y (1-based)")
    psr_scatter.add_argument('--color', type=int, default=-1, help="column of color (1-based)")
    psr_scatter.add_argument('-s', type=int, default=None, help="column to plot s (1-based, size of the markers)")
    psr_scatter.add_argument('--minsize', type=float, default=3.0, help='minimum size')
    psr_scatter.add_argument('--maxsize', type=float, default=20.0, help='maximum size')
    psr_scatter.add_argument('--ms', type=float, default=10, help="fixed marker size [10]")
    psr_scatter.add_argument('--xlog', action="store_true", help="x axis log-scaled")
    psr_scatter.add_argument('--ylog', action="store_true", help="y axis log-scaled")
    psr_scatter.add_argument('--beta', type=float, default=0, help='tilt angle of a dashed straight line (gamma + beta*(maxx-minx))')
    psr_scatter.add_argument('--gamma', type=float, default=0.0, help='intercept of the dashed straight line (gamma + beta*(maxx-minx))')
    add_std_options(psr_scatter)
    psr_scatter.set_defaults(func=main_scatter)

    ezscatter_parsers = subparsers.add_parser('ezscatter', help='easy scatter plots')
    ezscatter_subparsers = ezscatter_parsers.add_subparsers()
    psr_ezs_table = ezscatter_subparsers.add_parser('table', help="plot x vs multiple y with each y in a separate column (table format)")
    psr_ezs_table.add_argument('-ys', default='', help='column of ys (1-based)')
    psr_ezs_table.add_argument('-x', type=int, help="column of x (1-based)")
    add_std_options(psr_ezs_table)
    psr_ezs_table.set_defaults(func=main_ezscatter_table)
    
    # hexbin
    psr_hexbin = subparsers.add_parser("hexbin", help=""" hexbin plot """)
    psr_hexbin.add_argument('-x', type=int, help="column of x (1-based)")
    psr_hexbin.add_argument('-y', type=int, help="column of y (1-based)")
    psr_hexbin.add_argument('--nolog', action='store_true')
    psr_hexbin.add_argument('--bins', type=int, default=10, help='#bins along dimension')
    psr_hexbin.add_argument('--ybins', type=int, default=-1, help='if specified, used as number of bins in the y-axis')
    psr_hexbin.add_argument('--scatter', action='store_true', help='overlay scatter plot on top of hexbin')
    add_std_options(psr_hexbin)
    psr_hexbin.set_defaults(func=main_hexbin)

    # boxplot
    psr_box = subparsers.add_parser("box", help="box plot with discrete x-value")
    psr_box.add_argument('-c', type=int, help="column of value to plot [required]")
    psr_box.add_argument('-d', '--dim1', type=int, default=-1, help="column of dimention 1 category (1-based)")
    psr_box.add_argument('-d2', '--dim2', type=int, default=-1, help="column of dimention 2 category (1-based)")
    psr_box.add_argument('--width', type=int, default=0.5, help="width of each boxplot [0.5]")
    psr_box.add_argument('--catsep', type=int, default=0.5, help="separation between categories [0.1]")
    psr_box.add_argument('--inputraw', action='store_true', help='input format is raw "[1,2,3]" for column fields[args.c-1]')
    add_std_options(psr_box)
    psr_box.set_defaults(func=main_box)

    # boxbin
    psr_boxbin = subparsers.add_parser("boxbin", help="box plot binned by continuous value rather than from a given category")
    psr_boxbin.add_argument('-ys', help='column')
    psr_boxbin.add_argument('-x', type=int, help='value to be binned')
    psr_boxbin.add_argument('--mean', action='store_true', help='draw mean on each box')
    psr_boxbin.add_argument('--bins', type=int, default=15, help="number of bins [15]")
    psr_boxbin.add_argument('--catsep', type=int, default=0.5, help="separation between categories [0.1]")
    psr_boxbin.add_argument('--symb', default='+', help='symbol for outliers')
    add_std_options(psr_boxbin)
    psr_boxbin.set_defaults(func=main_boxbin)

    # histogram
    psr_hist = subparsers.add_parser("hist", help=""" histogram """)
    psr_hist.add_argument('-c', type=int, default=1, help="column to plot (1-based)")
    psr_hist.add_argument('--bins', type=int, default=20, help="number of bins [20]")
    psr_hist.add_argument('--cat', type=int, help="column to indicate category (1-based)")
    psr_hist.add_argument('--binseq', default=None, help='bin sequences. e.g., 1,20,100,400,1000 means [1,20),[20,100)...')
    psr_hist.add_argument('--log', action='store_true', help='log scaled x and y')
    psr_hist.add_argument('--xlog', action='store_true', help='log scaled x only')
    psr_hist.add_argument('--xrange', default=None, help="range of x axes. E.g., (xmin, xmax)")
    add_std_options(psr_hist)
    psr_hist.set_defaults(func=main_hist)

    parser_cumhist = subparsers.add_parser('cumhist', help=' cumulative histogram ')
    parser_cumhist.add_argument('-c', type=int, default=1, help="column to plot (1-based)")
    parser_cumhist.add_argument('--bins', type=int, default=20, help="number of bins [20]")
    add_std_options(parser_cumhist)
    parser_cumhist.set_defaults(func=main_cumhist)

    parser_dhist = subparsers.add_parser('dhist', help='discrete histogram')
    parser_dhist.add_argument('-i', type=argparse.FileType('r'), default='-', help='input table')
    parser_dhist.add_argument('-c', type=int, default=1, help="column to plot (1-based)")
    add_std_options(parser_dhist)
    parser_dhist.set_defaults(func=main_dhist)

    ezhist_parsers = subparsers.add_parser('ezhist', help='easy histogram plots')
    ezhist_subparsers = ezhist_parsers.add_subparsers()
    psr_ezh_sidebyside = ezhist_subparsers.add_parser('sidebyside')
    psr_ezh_sidebyside.add_argument('-c', type=int, help="column to plot (1-based)")
    psr_ezh_sidebyside.add_argument('--bins', type=int, default=10, help="number of bins [10]")
    psr_ezh_sidebyside.add_argument('--cat', type=int, help="column to indicate category (1-based)")
    add_std_options(psr_ezh_sidebyside)
    psr_ezh_sidebyside.set_defaults(func=main_ezhist_sidebyside)

    # pie
    psr_pie = subparsers.add_parser("pie", help=""" pie chart """)
    psr_pie.add_argument('-c', type=int, help="column to plot (1-based)")
    psr_pie.add_argument('--frac', type=int, default=1, help='index of fraction column')
    psr_pie.add_argument('--label', type=int, default=2, help='index of label column')
    add_std_options(psr_pie)
    psr_pie.set_defaults(func=main_pie)

    # bar
    psr_bar = subparsers.add_parser("bar", help=""" bar chart """)
    psr_bar.add_argument('-x', type=int, help="x axis column (1-based)")
    psr_bar.add_argument('-y', type=int, help='y axis column (1-based)')
    add_std_options(psr_bar)
    psr_bar.set_defaults(func=main_bar)

    ezbar_parsers = subparsers.add_parser('ezbar', help='easy bar plots')
    ezbar_subparsers = ezbar_parsers.add_subparsers()
    psr_ezb_sidebyside = ezbar_subparsers.add_parser('sidebyside')
    psr_ezb_sidebyside.add_argument('-x', type=int, help="column to plot bar labels (1-based), should contain unique text labels")
    psr_ezb_sidebyside.add_argument('-hs', help="column to plot bar heights (1-based), if only one column is given the plot is a simple bar plot, if more than one columns are given, the plot stacks multiple columns next to each other")
    psr_ezb_sidebyside.add_argument('--spacing', type=float, default=0.3, help='default[0.3]')
    psr_ezb_sidebyside.add_argument('--barwidth', type=float, default=1.5, help='default[1.5]')
    add_std_options(psr_ezb_sidebyside)
    psr_ezb_sidebyside.set_defaults(func=main_ezbar_sidebyside)

    psr_ezb_stack = ezbar_subparsers.add_parser('stack')
    psr_ezb_stack.add_argument('-x', type=int, help="column to plot bar labels (1-based), should contain unique text labels")
    psr_ezb_stack.add_argument('-hs', help="column to plot bar heights (1-based), if only one column is given the plot is a simple bar plot, if more than one columns are given, the plot stacks multiple columns on top of another")
    psr_ezb_stack.add_argument('--vspacing', type=float, default=0, help="vertical spacing")
    add_std_options(psr_ezb_stack)
    psr_ezb_stack.set_defaults(func=main_ezbar_stack)

    psr_ezb_substack = ezbar_subparsers.add_parser('substack')
    psr_ezb_substack.add_argument('-x', type=int, help="column to plot bar labels (1-based), should contain unique text labels")
    psr_ezb_substack.add_argument('-hs', help="column to plot bar heights (1-based), if only one column is given the plot is a simple bar plot, if more than one columns are given, the plot substacks multiple columns on top of another")
    add_std_options(psr_ezb_substack)
    psr_ezb_substack.set_defaults(func=main_ezbar_substack)


    # venn diagram
    psr_venn = subparsers.add_parser("venn", help=""" Venn's diagram (need matplotlib venn) """)
    psr_venn.add_argument('-n', type=int, default=2, help="number of ways of comparison. n=2 or 3")
    psr_venn.add_argument('-c', help="counts of each category. For 2-way venn, category Ab, aB, AB. For 3-way venn, category Abc, aBc, ABc, abC, AbC, aBC, ABC")
    psr_venn.add_argument('-l', default=None, help="labels of each category")
    add_std_options_out(psr_venn)
    psr_venn.set_defaults(func=main_venn)

    args = parser.parse_args()
    args.func(args)
