import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as mcolls
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
import numpy as np
import itertools
import wzcolors
import wzcore
import os
import pandas as pd

reload(wzcolors)
"""
low level plotting:
* plot_heatmap_discrete
* plot_heatmap
* 

high level plotting:
* WZHmap
* WZCbar
* WZCbarGroup
* WZCbarCollection
"""

def plot_heatmap_discrete(data, dim=None, fig=None, xlabels=None, ylabels=None, label_fontsize=8, level2color=None):

    """ plot heatmap for data frame with categorical value """

    # set defaults
    if dim is None:
        dim = [0.1,0.1,0.85,0.85]

    if fig is None:
        fig = plt.figure()

    if isinstance(fig, tuple):
        fig = plt.figure(figsize=fig)
        
    # map levels to color
    levels = sorted(list(set(data.values.flatten())))
    if level2color is None:
        colors = wzcolors.get_distinct_colors_rgb(len(levels))
        level2color = dict(zip(levels, colors))
    else:
        colors = [level2color[level] for level in levels]
    
    cmap = mcolors.ListedColormap(colors)
    bounds = np.arange(0.5, len(levels)+0.5, 1.0)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    centers = [_+0.5 for _ in bounds]
    level2center = dict(zip(levels, centers))

    dataplot = data.applymap(lambda x: level2center[x])

    ax = fig.add_axes(dim, frameon=False)
    ax.pcolor(data, cmap=cmap)
    # ax.imshow(dataplot, aspect='auto', origin='lower', cmap=cmap, interpolation='none')
    ax.set_xticks([])
    ax.set_yticks([])

    return ax, level2color

def plot_heatmap(data, dim=None, fig=None, xlabels=None, ylabels=None, label_fontsize=8,
                 axhlines=None, cmap=None, interpolation=None):

    """ plot heatmap for data frame with continuous value """
    
    # data = np.ma.array(data.as_matrix(), mask=data.isnull())

    # set defaults
    if dim is None:
        dim = [0.1,0.1,0.85,0.85]

    if fig is None:
        fig = plt.figure()

    if isinstance(fig, tuple):
        fig = plt.figure(figsize=fig)
        
    if cmap is None:
        cmap = cm.jet
    if isinstance(cmap, str):
        cmap = cm.get_cmap(cmap)
    cmap.set_bad('#E6E6E6', 1)

    ax = fig.add_axes(dim, frameon=False)
    # ax.pcolor(data, cmap=cmap)
    ax.imshow(data, aspect='auto', origin='lower', cmap=cmap, interpolation=interpolation)

    # ax.imshow(data, extent=[0,data.shape[0],0,data.shape[1]], aspect='equal', cmap=cmap, interpolation=interpolation)
    
    if axhlines is not None:
        for y in axhlines:
            ax.axhline(y, color='w', lw=0.2)

    # print ax.get_ylim()
    # print data.shape[1]
    # ax.set_ylim((0, data.shape[1]+1))
    ax.set_xticks([])
    ax.set_yticks([])

    return ax, cmap

def discrete_colorshow(data, dim, fig, orientation='horizontal',
                       greyscale=False, greyscale_range=(0.1,0.9), level2color=None):
    
    """ color bar of discrete values """

    n = len(data)
    if orientation == 'horizontal' or orientation == 'h':
        beg = '(i,0)'
        xmax = n
        ymax = 1
    elif orientation == 'vertical' or orientation == 'v':
        beg = '(0,i)'
        xmax = 1
        ymax = n
    else:
        raise Exception("unknown orientation")

    import wzcolors
    ax = fig.add_axes(dim, frameon=False)
    if level2color is None:     # if level2color is not given
        colors, level2color = wzcolors.map_distinct_colors_hex(
            data, other2grey=True, greyscale=greyscale, greyscale_range=greyscale_range)
    else:                       # if level2color is given
        colors = [level2color[datum] if datum in level2color else '#E6E6E6' for datum in data]
    color2patches = {}
    for i, color in enumerate(colors):
        if color not in color2patches:
            color2patches[color] = []
        color2patches[color].append(mpatches.Rectangle(eval(beg),1,1))

    for color, patches in color2patches.iteritems():
        ax.add_collection(mcolls.PatchCollection(patches, facecolor=color, edgecolor='none'))

    ax.set_xlim(0,xmax)
    ax.set_ylim(0,ymax)
    ax.set_xticks([])
    ax.set_yticks([])

    return ax, level2color

def colorshow_legend(l2c, (left,bottom,width,height), fig, label_fontsize=7,
                     horizontalspace=0.06, title=None, title_fontsize=7, alpha=1):

    """ legend of discrete color bar """

    n = len(l2c)
    ax = fig.add_axes((left, bottom, width, height), frameon=False)
    ax.set_xlim(0,1)
    ax.set_ylim(0,n)
    levels, colors = zip(*l2c.items())
    for i, color in enumerate(colors):
        ax.add_patch(mpatches.Rectangle((0,i),1,1,facecolor=color,edgecolor='white',alpha=alpha))
        ax.text(1+horizontalspace,i+0.5,str(levels[i]),verticalalignment="center", fontsize=label_fontsize, fontweight='light')
    # ax.set_yticks([i+0.5 for i in xrange(n)])
    # ax.set_yticklabels(levels, fontsize=fontsize, fontweight='light')
    # ax.yaxis.tick_right()
    ax.set_yticks([])
    ax.set_xticks([])
    if title:
        ax.set_title(title, fontsize=title_fontsize)

    return ax

def continuous_array_colorshow(data, dim=[0.1,0.1,0.85,0.85], fig=None, cmap='jet', orientation='horizontal'):

    """ color bar of continuous value """
    
    dmin = min(data)
    dmax = max(data)

    if cmap is None:
        colormap = cm.jet
    else:
        colormap = cm.get_cmap(cmap)
    colormap.set_bad('#E6E6E6', 1)
    
    norm = mcolors.Normalize(vmin=dmin, vmax=dmax)
    data_mat = np.array(data)
    if orientation == 'horizontal' or orientation == 'h':
        data_mat.shape = (1, len(data))
    if orientation == 'vertical' or orientation == 'v':
        data_mat.shape = (len(data), 1)

    if fig == None:
        fig = plt.figure()

    if isinstance(fig, tuple):
        fig = plt.figure(figsize=fig)


    ax = fig.add_axes(dim, frameon=False)
    ax.imshow(data_mat, cmap=colormap, aspect='auto', norm=norm, interpolation='none')
    ax.set_xticks([])
    ax.set_yticks([])
    
    return ax, norm, colormap

def continuous_colorshow_legend(norm, colormap, dim=[0.1,0.1,0.85,0.85], fig=None, title=None, fontsize=7, title_fontsize=7):

    """ plot legend of color bar with continuous values """
    if fig == None:
        fig = plt.figure()

    if isinstance(fig, tuple):
        fig = plt.figure(figsize=fig)

    ax = fig.add_axes(dim, frameon=False)
    cb = mcolorbar.ColorbarBase(ax, cmap=colormap, norm=norm)

    cb.ax.tick_params(labelsize=fontsize)
    # cb.outline.set_color('white') # this would wipe out the plot
    cb.outline.set_linewidth(0.2)
    # cb.dividers.set_color('white')
    # cb.dividers.set_linewidth(2)
    
    if title:
        ax.set_title(title, fontsize=title_fontsize)

    return ax

def fig_add_line(fig, xcoor, ycoor, color='k', linewidth=.5):
    fig.lines.append(mlines.Line2D(xcoor, ycoor, transform=fig.transFigure,
                                   figure=fig, color=color, linewidth=linewidth))

class WZHmap():

    def __init__(self, data,

                 continuous = True,

                 # discrete heat map
                 label2color = None,
                 greyscale = False,
                 greyscale_range = (0.1,0.9),

                 # continuous heat map
                 cmap = 'jet',
                 norm = None,
                 interpolation = 'none',
                 
                 # tick label on x axis
                 xticklabels = None,
                 xticklabel_fontsize = 5,
                 xticklabel_side = 'bottom',
                 xticklabel_rotat = 90,

                 yticklabels = None,
                 yticklabel_fontsize = 5,
                 yticklabel_side = 'l',

                 # title
                 label = None,
                 labelside = 'r',
                 labelspacing = 0.001,
                 labelfontsize = 8,
                 labelfontweight = 'light',

                 # legend
                 legend_title = None,
                 legend_title_fontsize = 8,
                 legend_label_fontsize = 8,
             ):
        
        self.ax = None
        self.data = data
        import inspect
        ia = inspect.getargspec(WZHmap.__init__)
        # print ia
        for a in ia.args[-len(ia.defaults):]:
            # print a, locals()[a]
            setattr(self, a, locals()[a])

        # print self.xticklabels

    def plot(self, dim=None, fig=None):

        # heat map
        if self.continuous:
            self.ax, self.colormap = plot_heatmap(self.data, dim, fig, interpolation=self.interpolation)
        else:
            self.ax, self.label2color = plot_heatmap_discrete(self.data, dim, fig, level2color=self.label2color)

        height = self.data.shape[0]
        # xtick labels
        if self.xticklabels is not None:
            if type(self.xticklabels) == bool:
                self.xticklabels = self.data.columns.format()
            for i in xrange(self.data.shape[1]):
                if self.xticklabel_side in ['bottom', 'b']:
                    self.ax.text(i, -1, self.xticklabels[i], rotation=self.xticklabel_rotat, horizontalalignment='center',
                                 verticalalignment='top', fontsize=self.xticklabel_fontsize)
                elif self.xticklabel_side in ['top', 't']:
                    self.ax.text(i, height*1.01, self.xticklabels[i], rotation=self.xticklabel_rotat, horizontalalignment='center',
                                 verticalalignment='bottom', fontsize=self.xticklabel_fontsize)

        if self.yticklabels is not None:
            if type(self.yticklabels) == bool:
                self.yticklabels = self.data.index.format()
            for i in xrange(self.data.shape[0]):
                self.ax.text(-1, i, self.yticklabels[i], horizontalalignment='right', fontsize=self.yticklabel_fontsize)

    def plot_legend(self, dim=[0.1,0.1,0.03,0.4], fig=None, unitheight=0.015):

        if fig is None:
            fig = plt.figure()

        if isinstance(fig, tuple):
            fig = plt.figure(figsize=fig)

        kwargs = {}
        if self.legend_title is not None:
            kwargs['title'] = self.legend_title
        if self.legend_title_fontsize is not None:
            kwargs['title_fontsize'] = self.legend_title_fontsize
        if self.legend_label_fontsize is not None:
            kwargs['label_fontsize'] = self.legend_label_fontsize

        if self.continuous:
            continuous_colorshow_legend(None, self.colormap, dim, fig, **kwargs)
        else:
            dim[-1] = unitheight*len(self.label2color)
            colorshow_legend(self.label2color, dim, fig, **kwargs)


def text_reconcile(pos1, interval, error=0.01):
    pos1 = sorted(pos1)
    # first resolve interval
    pos2 = [0] * len(pos1)
    for j in xrange(len(pos1)):
        if j == 0:
            pos2[0] = pos1[0]
            continue
        if pos1[j] - pos2[j-1] < interval:
            pos2[j] = pos2[j-1] + interval
        else:
            pos2[j] = pos1[j]
        
    concord = False
    ccc = 0
    # print 'pos1', pos1
    # print 'pos2 initial', pos2
    while not concord:
        concord = True
        
        # left shifts
        j = 0
        while j < len(pos1):
            group = [j]
            while j < len(pos1)-1 and pos2[j+1] <= pos2[j]+interval:
                j += 1
                group.append(j)
            j += 1
            mean1 = np.mean([pos1[_] for _ in group])
            mean2 = np.mean([pos2[_] for _ in group])
            if mean1 < mean2: # need to left shift
                delta_shift = mean2-mean1
                if group[0] != 0:
                    delta_shift = min(delta_shift, pos2[group[0]]-pos2[group[0]-1]-interval)
                if delta_shift < error:
                    continue
                for _j in group: # shift the entire group left
                    pos2[_j] -= delta_shift
                    concord = False

        #print 'left pos2', pos2
        # right shifts
        j = len(pos1)-1
        while j >= 0:
            group = [j]
            while j > 0 and pos2[j-1]+interval >= pos2[j]:
                j -= 1
                group.append(j)
            j -= 1
            mean1 = np.mean([pos1[_] for _ in group])
            mean2 = np.mean([pos2[_] for _ in group])
            if mean1 > mean2: # need to right shift
                delta_shift = mean1 - mean2
                if group[-1] != len(pos1)-1:
                    delta_shift = min(delta_shift, pos2[group[-1]+1]-pos2[group[-1]]-interval)
                if delta_shift < error:
                    continue
                for _j in group:  # shift the entire group right
                    pos2[_j] += delta_shift
                    concord = False
               
        #print 'right pos2', pos2
        ccc += 1
        if ccc > 3:
            break
    return pos2

class WZCbar(object):

    def __init__(self, data,
                 
                 continuous=False,

                 # discrete color bar
                 label2color = None,
                 greyscale = False,
                 greyscale_range = (0.1,0.9),

                 # continuous color bar
                 cmap = 'jet',
                 norm = None,

                 # title
                 title = None,
                 title_side = 'l',
                 title_spacing = 0.001,
                 title_fontsize = 8,
                 title_fontweight = 'light',

                 # slanted annotation
                 lineanno = None,
                 annhei = None,
                 annlft = None,
                 annlen = 0.03,
                 anntan = 1,

                 # legend
                 legend_title = None,
                 legend_title_fontsize = 8,
                 legend_label_fontsize = 8,

                 # perpendicular tick label
                 pticklabels = None,
                 pticklabel_fontsize = 4,
                 pticklabel_pad = 0.01,
             ):

        self.ax = None
        if isinstance(data, pd.Series):
            self.data = data
        else:
            self.data = pd.Series(data)
        import inspect
        ia = inspect.getargspec(WZCbar.__init__)
        for a in ia.args[-len(ia.defaults):]:
            setattr(self, a, locals()[a])
            
        if self.legend_title is None:
            self.legend_title = self.title

    def __len__(self):

        return 1

    def plot(self, dim, fig, orientation='horizontal'):

        left, bottom, width, height = dim
        top = bottom + height

        if self.continuous:
            self.ax, self.norm, self.colormap = continuous_array_colorshow(
                self.data, dim, fig, cmap=self.cmap)
        else:
            self.ax, self.label2color = discrete_colorshow(
                self.data, dim, fig, orientation=orientation, level2color=self.label2color)

        if self.lineanno == 'topleft' and self.title is not None:
            if self.annhei is None and self.annlft is None:
                self.annlft = left - 0.05
                self.annhei = top + 0.05
            elif self.annhei is None:
                self.annhei = top + (left + width / 2. - self.annlft) * self.anntan
            elif self.annlft is None:
                self.annlft = left + width / 2. - (self.annhei - top)

            anno_lw = 0.5
            # slant line
            fig_add_line(fig, [self.annlft, left+width/2.], [self.annhei, top], linewidth=anno_lw)
            # horizontal line
            fig_add_line(fig, [self.annlft-self.annlen, self.annlft],
                         [self.annhei, self.annhei], linewidth=anno_lw)
            fig.text(self.annlft, self.annhei, self.title,
                     fontsize = self.title_fontsize, fontweight = self.title_fontweight,
                     horizontalalignment='right', verticalalignment='bottom')

        elif self.title is not None:
            if self.title_side == 'right' or self.title_side == 'r':
                fig.text(left + width + self.title_spacing,
                         bottom + height / 2., self.title,
                         fontsize = self.title_fontsize, fontweight = self.title_fontweight,
                         horizontalalignment = 'left', verticalalignment = 'center')
            elif self.title_side == 'left' or self.title_side == 'l':
                fig.text(left - self.title_spacing, bottom + height / 2., self.title,
                         fontsize = self.title_fontsize, fontweight = self.title_fontweight,
                         horizontalalignment = 'right', verticalalignment = 'center')
            else:
                raise Exception('Unacceptable title side %s' % self.title_side)

        if self.pticklabels is not None:

            pl_inds = [_[0] for _ in self.pticklabels]
            pl_plots = text_reconcile(pl_inds, 10)
            ax = fig.add_axes([left-self.pticklabel_pad, bottom, self.pticklabel_pad, height], frameon=False)
            for _real, _plot in zip(pl_inds, pl_plots):
                ax.plot([1.0,0.6,0.4,0.], [_real, _real, _plot, _plot], color='k', lw=0.5)
            ax.set_ylim(0, len(self.data))
            ax.set_xlim(0,1)
            ax.set_yticks(pl_plots)
            ax.xaxis.set_tick_params(direction='in',width=0)
            ax.yaxis.set_tick_params(direction='in',width=0, pad=0)
            ax.set_yticklabels([_[1] for _ in self.pticklabels], fontsize=self.pticklabel_fontsize)
            ax.set_xticks([])

        return

    def plot_legend(self, dim=[0.1,0.1,0.03,0.4], fig=None, unitheight=0.015):

        """ plot legend of the color bar """

        if fig is None:
            fig = plt.figure()

        if isinstance(fig, tuple):
            fig = plt.figure(figsize=fig)

        kwargs = {}
        if self.legend_title is not None:
            kwargs['title'] = self.legend_title
        if self.legend_title_fontsize is not None:
            kwargs['title_fontsize'] = self.legend_title_fontsize
        if self.legend_label_fontsize is not None:
            kwargs['label_fontsize'] = self.legend_label_fontsize

        if self.continuous:
            continuous_colorshow_legend(self.norm, self.colormap, dim, fig,
                                        title = self.legend_title,
                                        title_fontsize = self.legend_title_fontsize)
        else:
            dim[-1] = unitheight*len(self.label2color)
            colorshow_legend(self.label2color, dim, fig, **kwargs)

        return

def cast_WZCbar(cbs):
    _cbs = []
    for cb in cbs:
        if not isinstance(cb, WZCbar):
            cb = WZCbar(cb)
        _cbs.append(cb)
    return _cbs


""" a group of color bars with the same color map """
class WZCbarGroup(object):

    def __init__(self,

                 continuous = False,
                 
                 # title
                 title = None,
                 title_side = 'r',
                 title_spacing = 0.001,
                 title_fontsize = 5,
                 title_fontweight = 'light',

                 # legend
                 legend_title = None,
                 legend_title_fontsize = 8,
                 legend_label_fontsize = 8,
             ):

        """ label2color is set at plot time """
        self.cbars = []
        self.label2color = None
        self.uid2cbar = {}


        import inspect
        ia = inspect.getargspec(WZCbarGroup.__init__)
        for a in ia.args[-len(ia.defaults):]:
            setattr(self, a, locals()[a])
            
        if self.legend_title is None:
            self.legend_title = self.title

    def getcat(self, catid):
        g2 = WZCbarGroup(title=self.title)
        g2.cbars = [c for c in self.cbars if c.cat==catid]
        g2.label2color = self.label2color
        return g2

    def __len__(self):
        
        return len(self.cbars)
    
    def add_cbar(self, data, uid=None, cat=1, **kwargs):

        cbar = WZCbar(data, **kwargs)
        self.cbars.append(cbar)
        cbar.cat = cat
        if uid is not None:
            self.uid2cbar[uid] = cbar

    def finalize_label2color(self, **kwargs):

        if self.label2color is None:
            levels = set()
            for cbar in self.cbars:
                levels |= set(cbar.data)
            self.label2color = wzcolors.map_level2color(levels, **kwargs)
            for cbar in self.cbars:
                cbar.label2color = self.label2color

    def finalize_sublabel2color(self, pg):

        pdata = []
        cdata = []
        for pcbar, ccbar in zip(pg.cbars, self.cbars):
            pdata.extend(pcbar.data)
            cdata.extend(ccbar.data)

        p2c = {}
        for _c, _p in zip(cdata, pdata):
            if _p in p2c:
                p2c[_p].add(_c)
            else:
                p2c[_p] = set()

        c2colors = wzcolors.map2sub_alpha(p2c, pg.label2color)
        self.label2color = c2colors
        for cbar in self.cbars:
            cbar.label2color = c2colors
        return


    def get(self, uid):

        if uid in self.uid2cbar:
            return self.uid2cbar[uid]
        else:
            return None
        
    # def plot(self, *args, **kwargs):
        
    #     """ plot all color bars """
    #     self.finalize_label2color(**kwargs)
    #     kwargs['label2color'] = self.label2color

    #     space = kwargs['space'] if 'space' in kwargs else 0
    #     if 'orientation' in kwargs and kwargs['orientation'] in ['v', 'vertical']:
    #         left, bottom, width, height = args[0]
    #         fig = args[1]
    #         for i, cbar in enumerate(self.wcbars):
    #             dim = (left + i*(width+space), bottom, width, height)
    #             cbar.plot(dim, fig, **kwargs)

    def plot_cbar(self, label, *args, **kwargs):
        cbar = self.label2cbar[label]
        self.finalize_label2color(**kwargs)
        kwargs['label2color'] = self.label2color
        cbar.plot(*args, **kwargs)

    def plot_legend(self, dim=[0.1,0.1,0.03,0.4], fig=None, unitheight=0.015):

        if fig is None:
            fig = plt.figure()

        if isinstance(fig, tuple):
            fig = plt.figure(figsize=fig)

        kwargs = {}
        if self.legend_title is not None:
            kwargs['title'] = self.legend_title
        if self.legend_title_fontsize is not None:
            kwargs['title_fontsize'] = self.legend_title_fontsize
        if self.legend_label_fontsize is not None:
            kwargs['label_fontsize'] = self.legend_label_fontsize

        if self.continuous:
            continuous_colorshow_legend(self.norm, self.colormap, dim, fig,
                                        title = self.legend_title,
                                        title_fontsize = self.legend_title_fontsize)
        else:
            dim[-1] = unitheight*len(self.label2color)
            colorshow_legend(self.label2color, dim, fig, **kwargs)


            
        # height = len(self.label2color)*patchheight

        # colorshow_legend(self.label2color, (left, bottom, width, height), fig, title=self.title)

        # self.leghei = height
        return

class WZCbarCollection():

    """ a collection of color bars or color bar groups """

    def __init__(self, colls=[]):

        self.colls = colls

    def __len__(self):

        return sum([len(coll) for coll in self.colls])

    def add_cbar(self, cbar):
        self.colls.append(cbar)

    def cbars(self):

        for cb in self.colls:
            if isinstance(cb, WZCbarGroup):
                for _cb in cb.cbars:
                    yield _cb
            elif isinstance(cb, WZCbar):
                yield cb

def scatter_matrix(frame, alpha=0.5, figsize=None, ax=None, grid=False,
                   diagonal='hist', marker='.', density_kwds=None,
                   hist_kwds=None, range_padding=0.05, rot=0, **kwds):
    """
    Draw a matrix of scatter plots.

    Parameters
    ----------
    frame : DataFrame
    alpha : float, optional
        amount of transparency applied
    figsize : (float,float), optional
        a tuple (width, height) in inches
    ax : Matplotlib axis object, optional
    grid : bool, optional
        setting this to True will show the grid
    diagonal : {'hist', 'kde'}
        pick between 'kde' and 'hist' for
        either Kernel Density Estimation or Histogram
        plot in the diagonal
    marker : str, optional
        Matplotlib marker type, default '.'
    hist_kwds : other plotting keyword arguments
        To be passed to hist function
    density_kwds : other plotting keyword arguments
        To be passed to kernel density estimate plot
    range_padding : float, optional
        relative extension of axis range in x and y
        with respect to (x_max - x_min) or (y_max - y_min),
        default 0.05
    kwds : other plotting keyword arguments
        To be passed to scatter function

    Examples
    --------
    >>> df = DataFrame(np.random.randn(1000, 4), columns=['A','B','C','D'])
    >>> scatter_matrix(df, alpha=0.2)
    """
    import matplotlib.pyplot as plt
    from matplotlib.artist import setp
    import pandas.tools.plotting as pdtp
    import pandas.core.common as com
    import numpy as np
    import pandas.compat as pdcompat

    df = frame._get_numeric_data()
    n = df.columns.size
    naxes = n * n
    fig, axes = pdtp._subplots(naxes=naxes, figsize=figsize, ax=ax,
                          squeeze=False)

    # no gaps between subplots
    fig.subplots_adjust(wspace=0, hspace=0)

    mask = com.notnull(df)

    marker = pdtp._get_marker_compat(marker)

    hist_kwds = hist_kwds or {}
    density_kwds = density_kwds or {}

    # workaround because `c='b'` is hardcoded in matplotlibs scatter method
    kwds.setdefault('c', plt.rcParams['patch.facecolor'])

    boundaries_list = []
    for a in df.columns:
        values = df[a].values[mask[a].values]
        rmin_, rmax_ = np.min(values), np.max(values)
        rdelta_ext = (rmax_ - rmin_) * range_padding / 2.
        boundaries_list.append((rmin_ - rdelta_ext, rmax_+ rdelta_ext))

    for i, a in zip(pdcompat.lrange(n), df.columns):
        for j, b in zip(pdcompat.lrange(n), df.columns):
            ax = axes[i, j]

            if i == j:
                values = df[a].values[mask[a].values]

                # Deal with the diagonal by drawing a histogram there.
                if diagonal == 'hist':
                    ax.hist(values, **hist_kwds)

                elif diagonal in ('kde', 'density'):
                    from scipy.stats import gaussian_kde
                    y = values
                    gkde = gaussian_kde(y)
                    ind = np.linspace(y.min(), y.max(), 1000)
                    ax.plot(ind, gkde.evaluate(ind), **density_kwds)

                ax.set_xlim(boundaries_list[i])

            else:
                common = (mask[a] & mask[b]).values

                ax.scatter(df[b][common], df[a][common],
                           marker=marker, alpha=alpha, **kwds)

                ax.set_xlim(boundaries_list[j])
                ax.set_ylim(boundaries_list[i])

            ax.set_xlabel('')
            ax.set_ylabel('')

            _label_axis(ax, kind='x', label=b, position='bottom', rotate=True, labelrot=rot)

            _label_axis(ax, kind='y', label=a, position='left', labelrot=rot)

            if j!= 0:
                ax.yaxis.set_visible(False)
            if i != n-1:
                ax.xaxis.set_visible(False)

    for ax in axes.flat:
        setp(ax.get_xticklabels(), fontsize=8)
        setp(ax.get_yticklabels(), fontsize=8)

    return axes

def _label_axis(ax, kind='x', label='', position='top',
                ticks=True, rotate=False, labelrot=0):

    from matplotlib.artist import setp
    if kind == 'x':
        ax.set_xlabel(label, visible=True, rotation=labelrot)
        ax.xaxis.set_visible(True)
        ax.xaxis.set_ticks_position(position)
        ax.xaxis.set_label_position(position)
        if rotate:
            setp(ax.get_xticklabels(), rotation=90)
    elif kind == 'y':
        ax.yaxis.set_visible(True)
        ax.set_ylabel(label, visible=True, rotation=90-labelrot)
        # ax.set_ylabel(a)
        ax.yaxis.set_ticks_position(position)
        ax.yaxis.set_label_position(position)
        
    return

def cb_reorder_generic(cb, order):

    """ utility function for reorder data within color bar """

    if isinstance(cb, WZCbar):
        if order is not None:
            cb.data = cb.data.iloc[order] # in-situ modification of cb
    elif isinstance(cb, pd.Series):
        if order is not None:
            cb = cb.iloc[order]
    else:
        cb = pd.Series(cb)
        if order is not None:
            cb = cb.iloc[order]

    return cb

def hmap_cluster_generic(hmap, fast=True, mode='b'):

    import wzhierarchy
    if mode == 'b':
        clustfun = wzhierarchy.ez_cluster
    elif mode == 'r':
        clustfun = wzhierarchy.ez_cluster_row
    elif mode == 'c':
        clustfun = wzhierarchy.ez_cluster_column
    
    if isinstance(hmap, WZHmap):
        cd = clustfun(hmap.data, fast=fast)
        hmap.data = cd.df
    elif isinstance(hmap, pd.DataFrame):
        cd = clustfun(hmap, fast=fast)
        hmap = WZHmap(cd.df)
    else:                       # right now nothing different
        cd = clustfun(hmap, fast=fast)
        hmap = WZHmap(cd.df)

    return hmap, cd

def subset_kwargs(kwargs, keys):

    kwargs2 = {}
    for key in keys:
        if key in kwargs:
            kwargs2[key] = kwargs[key]

    return kwargs2

def single_cluster_layout(hmap, lcbs=[], tcbs=[], fast=True, mode='b', **kwargs):

    hmap, cd = hmap_cluster_generic(hmap, fast=fast, mode=mode)
    if 'td' in kwargs and kwargs['td']:
        kwargs['td'] = cd.D_top
    single_hmap_layout(hmap,
                       lcbs=[cb_reorder_generic(cb, cd.lft_order()) for cb in lcbs],
                       tcbs=[cb_reorder_generic(cb, cd.top_order()) for cb in tcbs], **kwargs)

def row_stack_layout(cbs, figwid=10, fighei=10, figfile=None,
                     pad=0.002, hei=0.01,
                     lft=0.05, btm=0.05, wid=0.9):

    cbs = cast_WZCbar(cbs)
    fig = plt.figure(figsize=(figwid, fighei))
    for i, cb in enumerate(cbs):
        cb.plot([lft, btm + i*(hei+pad), wid, hei], fig)

    if figfile is not None:
        wzcore.err_print("Saving to %s." % figfile)
        fig.savefig(figfile, bbox_inches='tight', dpi=150)
        
def single_hmap_layout(hmap, lcbs=[], tcbs=[], td=None, ld=None,
                       figwid=10, fighei=10, figfile=None,
                       tdhei = 0.1, tdpad = 0.005, ldwid = 0.1, ldpad = 0.003,
                       lcwid = 0.01, lcpad = 0.001, tchei = 0.02, tcpad = 0.001,
                       maheipad = 0.001,
                       **kwargs):
    
    fig = plt.figure(figsize=(figwid, fighei))
    fig.patch.set_facecolor('white')

    ## plot margin
    mar = 0.05

    ## top dendrogram
    if td is None:
        tdhei = 0.1
        tdpad = 0
    tduni = tdhei + tdpad

    ## left dendrogram
    if ld is None:
        ldwid = 0
        ldpad = 0
    lduni = ldwid + ldpad

    ## left colorbar
    lcuni = lcwid + lcpad
    lcwid_all = len(lcbs) * lcwid + (((len(lcbs) - 1) * lcpad) if len(lcbs)>1 else 0)

    ## top color bar
    tcuni = tchei + tcpad
    tchei_all = len(tcbs) * tchei + (((len(tcbs) - 1) * tcpad) if len(tcbs)>1 else 0)

    ## matrix
    mawid = 1 - 2*mar - lduni - lcwid_all
    mahei = 1 - 2*mar - tduni - tchei_all
    maheiuni = mahei + maheipad

    assert(mawid>0)
    assert(mahei>0)

    if ld is not None:
        wzcore.err_print("left dendro")
        ld.plot(fig, [mar, mar, ldwid, mahei], orientation="left")
        lclft = mar + lduni
    else:
        lclft = mar

    if len(lcbs) > 0:
        wzcore.err_print("left color")
        lcb = cast_WZCbar(lcbs)
        for i, cbar in enumerate(lcbs):
            cbar.lineanno = 'topleft'
            cbar.annlft = lclft
            cbar.plot([lclft+i*lcuni, mar, lcwid, mahei], fig, orientation='v')
    malft = lclft + len(lcbs) * lcuni

    wzcore.err_print("heatmap")
    if not isinstance(hmap, WZHmap):
        hmap = WZHmap(hmap)
    hmap.plot(dim=[malft, mar, mawid, mahei], fig=fig)

    tcbtm = mar + maheiuni

    if len(tcbs) > 0:
        wzcore.err_print("top color")
        tcbs = cast_WZCbar(tcbs)
        for i, cbar in enumerate(tcbs):
            cbar.plot([malft, tcbtm+i*tcuni, mawid, tchei], fig)
        tdbtm = tcbtm + len(tcbs)*tcuni
    else:
        tdbtm = tcbtm

    if td is not None:
        wzcore.err_print("top dendro")
        td.plot(fig, [malft, tdbtm, mawid, tdhei], orientation='top')

    if figfile is not None:
        wzcore.err_print("Saving to %s." % figfile)
        fig.savefig(figfile, bbox_inches='tight', dpi=150)

    # plot legend separately
    if lcbs or tcbs and figfile is not None:
        wzcore.err_print('Legend for color bars')
        legfigwid = 10
        legfighei = 8
        legbtm = 0.1
        leglft = 0.1
        legpad = 1.0 / legfigwid
        legwid = 0.6 / legfigwid

        fig = plt.figure(figsize=(10,8))
        for cbar in itertools.chain(lcbs, tcbs):
            cbar.plot_legend([leglft, legbtm, legwid, 0.4], fig)
            leglft += legwid + legpad

        ff = os.path.splitext(figfile)
        legfigfile = ff[0]+'_legend'+ff[1]
        fig.savefig(legfigfile, bbox_inches='tight', dpi=150)
        wzcore.err_print('Saving legend..')
        wzcore.err_print("Done.")

    return

def normal_tumor_layout(cd_tumor, cd_normal=None, top_cbars=WZCbarCollection(), left_cbars=WZCbarCollection(),
                        top_normal_cbars=WZCbarCollection(), 
                        figwid=18, fighei=8, figfile="tmp.png", **kwargs):

    fig = plt.figure(figsize=(figwid, fighei))
    fig.patch.set_facecolor('white')

    mar = 0.05
    tdhei = 0.1                 # top dendro height
    tdpad = 0.005
    tduni = tdhei + tdpad

    ldwid = 0.1                 # left dendro width
    ldpad = 0.003
    lduni = ldwid + ldpad

    nmwid = 0.1 if cd_normal is not None else 0
    nmpad = 0.001 if cd_normal is not None else 0

    nmuni = nmwid + nmpad
    lcwid = 0.01                # left colorbar width
    lcpad = 0.001               # left colorbar spacing
    lcuni = lcwid + lcpad       # left colorbar unit width

    tchei = 0.02  # top colorbar height
    if 'tcpad' in kwargs and kwargs['tcpad']:
        tcpad = kwargs['tcpad']
    else:
        tcpad = 0.001 # top colorbar spacing
    tcuni = tchei + tcpad       # top colorbar unit height

    # print left_cbars, top_cbars
    mawid = (1 - 2*mar - ldpad - ldwid - nmwid - 2*nmpad - ((len(left_cbars)-1)*lcpad if len(left_cbars)>1 else 0) - len(left_cbars)*lcwid)
    mahei = (1 - 2*mar - tdpad - tdhei - ((len(top_cbars)-1)*lcpad if len(top_cbars)>1 else 0) - len(left_cbars)*lcwid)
    maheipad = 0.001
    maheiuni = mahei + maheipad

    figfile = figfile

    assert(mawid>0)
    assert(mahei>0)

    print "left dendro"
    cd_tumor.D_lft.plot(fig, [mar, mar, ldwid, mahei], orientation="left")

    print "left color"
    lclft = mar+lduni
    for i, cbar in enumerate(left_cbars.cbars()):
        cbar.lineanno = 'topleft'
        cbar.annlft = lclft
        cbar.plot([lclft+i*lcuni, mar, lcwid, mahei], fig, orientation='v')

    print "normal heatmap"
    nmlft = lclft + len(left_cbars)*lcuni
    if cd_normal is not None:
        plot_heatmap(cd_normal.df, [nmlft, mar, nmwid, mahei], fig, interpolation='none')

    print "heatmap"
    malft = nmlft + nmuni
    plot_heatmap(cd_tumor.df, [malft, mar, mawid, mahei], fig, interpolation='none')

    print "top normal color"
    tcbtm = mar + maheiuni
    if top_normal_cbars.colls:
        for i, top_normal_cbar in enumerate(top_normal_cbars.cbars()):
            top_normal_cbar.plot([nmlft, tcbtm+i*tcuni, nmwid, tchei], fig)

    print "top color"
    if top_cbars.colls:
        for i, cbar in enumerate(top_cbars.cbars()):
            cbar.plot([malft, tcbtm+i*tcuni, mawid, tchei], fig)

    print "top dendro"
    tdbtm = tcbtm + len(top_cbars)*tcuni
    cd_tumor.D_top.plot(fig, [malft, tdbtm, mawid, tdhei], orientation='top')

    print "Saving.."
    fig.savefig(figfile, bbox_inches='tight', dpi=150)

    # plot legend separately
    print 'Legend'
    legfigwid = 10
    legfighei = 8
    legbtm = 0.1
    leglft = 0.1
    legpad = 1.0 / legfigwid
    legwid = 0.6 / legfigwid
    fig = plt.figure(figsize=(10,8))
    for cbar in left_cbars.colls+top_cbars.colls:
        cbar.plot_legend([leglft, legbtm, legwid, 0.4], fig)
        leglft += legwid + legpad
    ff = os.path.splitext(figfile)
    legfigfile = ff[0]+'_legend'+ff[1]
    fig.savefig(legfigfile, bbox_inches='tight', dpi=150)
    print 'Saving legend..'
    print "Done."

def normal_tumor_layout_simple(df_tumor, df_normal=None, top_cbars=WZCbarCollection(), left_cbars=WZCbarCollection(),
                               top_normal_cbars=WZCbarCollection(), tumor_axhlines=None, dpi=150,
                               figwid=18, fighei=8, figfile="tmp.png", **kwargs):

    """ version without dendrogram """
    fig = plt.figure(figsize=(figwid, fighei))
    fig.patch.set_facecolor('white')

    left_cbars = left_cbars
    top_cbars = top_cbars
    top_normal_cbars = top_normal_cbars

    mar = 0.05

    nmwid = 0.1 if df_normal is not None else 0
    nmpad = 0.001 if df_normal is not None else 0

    nmuni = nmwid + nmpad
    lcwid = 0.01                # left colorbar width
    lcpad = 0.001               # left colorbar spacing
    lcuni = lcwid + lcpad       # left colorbar unit width

    tchei = 0.02  # top colorbar height
    if 'tcpad' in kwargs and kwargs['tcpad']:
        tcpad = kwargs['tcpad']
    else:
        tcpad = 0.001 # top colorbar spacing
    tcuni = tchei + tcpad       # top colorbar unit height

    # print left_cbars, top_cbars
    mawid = (1 - 2*mar - nmwid - 2*nmpad - ((len(left_cbars)-1)*lcpad if len(left_cbars)>1 else 0) - len(left_cbars)*lcwid)
    mahei = (1 - 2*mar - ((len(top_cbars)-1)*lcpad if len(top_cbars)>1 else 0) - len(left_cbars)*lcwid)
    maheipad = 0.001
    maheiuni = mahei + maheipad

    figfile = figfile

    assert(mawid>0)
    assert(mahei>0)

    print "left color"
    lclft = mar
    for i, cbar in enumerate(left_cbars.cbars()):
        cbar.lineanno = 'topleft'
        cbar.annlft = lclft
        cbar.plot([lclft+i*lcuni, mar, lcwid, mahei], fig, orientation='v')

    print "normal heatmap"
    nmlft = lclft + len(left_cbars)*lcuni
    if df_normal is not None:
        plot_heatmap(df_normal, [nmlft, mar, nmwid, mahei], fig, axhlines=tumor_axhlines, interpolation='none')

    print "heatmap"
    malft = nmlft + nmuni
    plot_heatmap(df_tumor, [malft, mar, mawid, mahei], fig, axhlines=tumor_axhlines, interpolation='none')
                 # interpolation='nearest')

    print "top normal color"
    tcbtm = mar + maheiuni
    if top_normal_cbars.colls:
        for i, top_normal_cbar in enumerate(top_normal_cbars.cbars()):
            top_normal_cbar.plot([nmlft, tcbtm+i*tcuni, nmwid, tchei], fig)

    print "top color"
    if top_cbars.colls:
        for i, cbar in enumerate(top_cbars.cbars()):
            cbar.plot([malft, tcbtm+i*tcuni, mawid, tchei], fig)


    print "Saving.."
    fig.savefig(figfile, bbox_inches='tight', dpi=dpi)

    # plot legend separately
    print 'Legend'
    legfigwid = 10
    legfighei = 8
    legbtm = 0.1
    leglft = 0.1
    legpad = 1.0 / legfigwid
    legwid = 0.6 / legfigwid
    fig = plt.figure(figsize=(10,8))
    for cbar in left_cbars.colls+top_cbars.colls:
        cbar.plot_legend([leglft, legbtm, legwid, 0.4], fig)
        leglft += legwid + legpad
    ff = os.path.splitext(figfile)
    legfigfile = ff[0]+'_legend'+ff[1]
    fig.savefig(legfigfile, bbox_inches='tight', dpi=150)
    print 'Saving legend..'
    print "Done."


def normal_tumor_with_contrast_layout(cd_tumor, cd_normal=None, cd_tumor_contrast=None, cd_normal_contrast=None,
                                      top_cbars=[], left_cbars=[], left_cbars_contrast=[], top_normal_cbars=[],
                                      contrast_hei_ratio = 0.5, figwid=18, fighei=10, figfile="tmp.png", **kwargs):

    fig = plt.figure(figsize=(figwid, fighei))
    fig.patch.set_facecolor('white')

    left_cbars = left_cbars
    top_cbars = top_cbars
    top_normal_cbars = top_normal_cbars
    cd_tumor = cd_tumor
    cd_normal = cd_normal

    mar = 0.05
    tdhei = 0.1                 # top dendro height
    tdpad = 0.005
    tduni = tdhei + tdpad

    ldwid = 0.1                 # left dendro width
    ldpad = 0.003
    lduni = ldwid + ldpad

    nmwid = 0.1
    nmpad = 0.003

    nmuni = nmwid + nmpad
    lcwid = 0.01                # left colorbar width
    lcpad = 0.001               # left colorbar spacing
    lcuni = lcwid + lcpad       # left colorbar unit width

    tchei = 0.02  # top colorbar height
    if 'tcpad' in kwargs and kwargs['tcpad']:
        tcpad = kwargs['tcpad']
    else:
        tcpad = 0.001 # top colorbar spacing
    tcuni = tchei + tcpad       # top colorbar unit height

    # print left_cbars, top_cbars
    mawid = (1 - 2*mar - ldpad - ldwid - nmwid - 2*nmpad - ((len(left_cbars)-1)*lcpad if len(left_cbars)>1 else 0) - len(left_cbars)*lcwid)
    mahei = (1 - 2*mar - tdpad - tdhei - ((len(top_cbars)-1)*lcpad if len(top_cbars)>1 else 0) - len(left_cbars)*lcwid)
    # mahei /= 2
    mahei_c = mahei * contrast_hei_ratio
    mahei = mahei - mahei_c
    maheipad = 0.001

    figfile = figfile

    assert(mawid>0)
    assert(mahei>0)

    print 'left dendro contrast'
    cd_tumor_contrast.D_lft.plot(fig, [mar, mar, ldwid, mahei_c], orientation='left')

    print 'left color bars contrast'
    lclft = mar+lduni
    for i, cbar in enumerate(left_cbars_contrast.cbars()):
        cbar.plot([lclft+i*lcuni, mar, lcwid, mahei_c], fig, orientation='v')

    print 'normal heatmap contrast'
    nmlft = lclft + len(left_cbars)*lcuni
    plot_heatmap(cd_normal_contrast.df, [nmlft, mar, nmwid, mahei_c], fig, interpolation='none')

    print 'tumor heatmap contrast'
    malft = nmlft + nmuni
    plot_heatmap(cd_tumor_contrast.df, [malft, mar, mawid, mahei_c], fig, interpolation='none')
    
    print "left dendro"
    ldbtm = mar + mahei_c + maheipad
    cd_tumor.D_lft.plot(fig, [mar, ldbtm, ldwid, mahei], orientation="left")

    print "left color"
    for i, cbar in enumerate(left_cbars.cbars()):
        cbar.lineanno = 'topleft'
        cbar.annlft = lclft
        cbar.plot([lclft+i*lcuni, ldbtm, lcwid, mahei], fig, orientation='v')

    print "normal heatmap"
    plot_heatmap(cd_normal.df, [nmlft, ldbtm, nmwid, mahei], fig, interpolation='none')

    print "tumor heatmap"
    plot_heatmap(cd_tumor.df, [malft, ldbtm, mawid, mahei], fig, interpolation='none')

    print "top normal color"
    tcbtm = ldbtm + mahei + maheipad
    for i, top_normal_cbar in enumerate(top_normal_cbars.cbars()):
        top_normal_cbar.plot([nmlft, tcbtm+i*tcuni, nmwid, tchei], fig)

    print "top color"
    for i, cbar in enumerate(top_cbars.cbars()):
        cbar.plot([malft, tcbtm+i*tcuni, mawid, tchei], fig)

    print "top dendro"
    tdbtm = tcbtm + len(top_cbars)*tcuni
    cd_tumor.D_top.plot(fig, [malft, tdbtm, mawid, tdhei], orientation='top')

    print "Saving.."
    fig.savefig(figfile, bbox_inches='tight', dpi=150)

    # plot legend separately
    print 'Legend'
    legfigwid = 10
    legfighei = 8
    legbtm = 0.1
    leglft = 0.1
    legpad = 1.0 / legfigwid
    legwid = 0.6 / legfigwid
    fig = plt.figure(figsize=(10,8))
    for cbar in left_cbars.colls+top_cbars.colls:
        cbar.plot_legend([leglft, legbtm, legwid, 0.4], fig)
        leglft += legwid + legpad
    ff = os.path.splitext(figfile)
    legfigfile = ff[0]+'_legend'+ff[1]
    fig.savefig(legfigfile, bbox_inches='tight', dpi=150)
    print 'Saving legend..'
    print "Done."


def savefig(fn):
    plt.savefig(fn, bbox_inches='tight', dpi=150)

def violin(data, colors, labels=None, bw=0.4, cut=2, gridsize=1000, poses=None, dwidth=0.5):
    # res = plt.violinplot(data, showextrema=False, bw_method="silverman")
    # res = plt.violinplot(data, showextrema=False, bw_method="scott")
    # res = plt.violinplot(data, showextrema=False, bw_method=0.3)
    # for i, r in enumerate(res['bodies']):
    #     r.set_facecolor(colors[i])
    #     r.set_edgecolor('none')
    #     r.set_alpha(0.9)

    ax = plt.gca()
    from scipy.stats import gaussian_kde
    supports = []
    densities = []
    max_dens = 0
    for i, datum in enumerate(data):
        kde = gaussian_kde(datum, bw)
        bw2 = bw*np.std(datum, ddof=1)
        support_min = np.min(datum) - bw2 * cut
        support_max = np.max(datum) + bw2 * cut
        support = np.linspace(support_min, support_max, gridsize)
        density = kde.evaluate(support)
        supports.append(support)
        max_dens = max(np.max(density), max_dens)
        densities.append(density)

    for i, density in enumerate(densities):
        density /= max_dens
        grid = poses[i] if poses is not None else i+1
        ax.fill_betweenx(supports[i], grid - density*dwidth, grid + density*dwidth, facecolor=colors[i], edgecolor="none")

    if labels is not None:
        plt.xticks(range(1,len(labels)+1), labels)

    if poses is None:
        plt.xlim(0, len(densities)+0.5)
    else:
        plt.xlim(min(poses)-0.5, max(poses)+0.5)


