import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as mcolls
import matplotlib.colorbar as mcolorbar
import matplotlib.lines as mlines
import numpy as np
import wzcolors
import os

def plot_heatmap(data, dim=[0.1,0.1,0.85,0.85], fig=None, xlabels=None, ylabels=None, label_fontsize=8, cmap=None, interpolation=None):

    # data = np.ma.array(data.as_matrix(), mask=data.isnull())
    
    if cmap is None:
        cmap = cm.jet

    cmap.set_bad('#E6E6E6', 1)
    if fig is None:
        fig = plt.figure()
            
    ax = fig.add_axes(dim, frameon=False)
    ax.imshow(data, aspect='auto', origin='lower', cmap=cmap, interpolation=interpolation)

    ax.set_xticks([])
    ax.set_yticks([])
    if xlabels is not None:
        if type(xlabels) == bool:
            xlabels = data.columns.format()
        for i in xrange(data.shape[1]):
            ax.text(i+0.2, -1, xlabels[i], rotation=90, horizontalalignment='center', verticalalignment='top', fontsize=label_fontsize)

    return ax

def colorshow_legend(l2c, (left,bottom,width,height), fig, fontsize=7,
                     horizontalspace=0.05, title=None, title_fontsize=7):

    n = len(l2c)
    ax = fig.add_axes((left, bottom, width, height), frameon=False)
    ax.set_xlim(0,1)
    ax.set_ylim(0,n)
    levels, colors = zip(*l2c.items())
    for i, color in enumerate(colors):
        ax.add_patch(mpatches.Rectangle((0,i),1,1,facecolor=color,edgecolor='white'))
        ax.text(1+horizontalspace,i+0.5,str(levels[i]),verticalalignment="center", fontsize=fontsize, fontweight='light')
    # ax.set_yticks([i+0.5 for i in xrange(n)])
    # ax.set_yticklabels(levels, fontsize=fontsize, fontweight='light')
    # ax.yaxis.tick_right()
    ax.set_yticks([])
    ax.set_xticks([])
    if title:
        ax.set_title(title, fontsize=title_fontsize)

    return ax

def discrete_array_colorshow(data, dim, fig, orientation='horizontal',
                             greyscale=False, greyscale_range=(0.1,0.9), level2color=None):

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

def continuous_array_colorshow(data, dim=[0.1,0.1,0.85,0.85], fig=None, cmap='jet', orientation='horizontal'):

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

    ax = fig.add_axes(dim, frameon=False)
    # ax.pcolor(data_mat, cmap=colormap, norm=norm)#, interpolation='none')
    ax.imshow(data_mat, cmap=colormap, aspect='auto', norm=norm, interpolation='none')
    ax.set_xticks([])
    ax.set_yticks([])
    
    return ax, norm, colormap

def continuous_array_colorshow_legend(norm, colormap, dim=[0.1,0.1,0.85,0.85], fig=None, title=None, fontsize=7, title_fontsize=7):

    if fig == None:
        fig = plt.figure()

    ax = fig.add_axes(dim, frameon=False)
    cb = mcolorbar.ColorbarBase(ax, cmap=colormap, norm=norm)
    
    cb.ax.tick_params(labelsize=fontsize)
    cb.outline.set_color('white')
    cb.outline.set_linewidth(2)
    # cb.dividers.set_color('white')
    # cb.dividers.set_linewidth(2)
    
    if title:
        ax.set_title(title, fontsize=title_fontsize)

    return ax

def fig_add_line(fig, xcoor, ycoor, color='k', linewidth=.5):
    fig.lines.append(mlines.Line2D(xcoor, ycoor, transform=fig.transFigure,
                                   figure=fig, color=color, linewidth=linewidth))

class WCbar(object):

    def __init__(self, data, continuous=False, label=None, cmap='jet',
                 greyscale=False, greyscale_range=(0.1,0.9), labelside='r',
                 labelspacing=0.001, labelfontsize=5, labelfontweight='light',
                 annhei=None, annlft=None, annlen=0.03, lineanno=None, anntan=1):

        self.continuous = continuous
        self.label = label
        self.data = data
        self.ax = None

        # for discrete color map
        self.label2color = None
        self.greyscale = greyscale
        self.greyscale_range = greyscale_range

        # for continuous color map
        self.cmap = cmap
        self.norm = None
        self.colormap = None

        # labels
        self.labelside = labelside
        self.labelspacing = labelspacing
        self.labelfontsize = labelfontsize
        self.labelfontweight = labelfontweight

        # slanted annotation
        self.lineanno = lineanno
        self.annhei = annhei
        self.annlft = annlft
        self.annlen = annlen
        self.anntan = anntan

    def __len__(self):

        return 1

    def plot(self, dim, fig, orientation='horizontal'):

        left, bottom, width, height = dim
        top = bottom + height

        if self.continuous:
            self.ax, self.norm, self.colormap = continuous_array_colorshow(
                self.data, dim, fig, cmap=self.cmap)
        else:
            self.ax, self.label2color = discrete_array_colorshow(
                self.data, dim, fig, orientation=orientation, level2color=self.label2color)

        if self.lineanno == 'topleft':
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
            fig.text(self.annlft, self.annhei, self.label,
                     fontsize = self.labelfontsize, fontweight = self.labelfontweight,
                     horizontalalignment='right', verticalalignment='bottom')

        elif self.label is not None:
            if self.labelside == 'right' or self.labelside == 'r':
                fig.text(left + width + self.labelspacing,
                         bottom + height / 2., self.label,
                         fontsize = self.labelfontsize, fontweight = self.labelfontweight,
                         horizontalalignment = 'left', verticalalignment = 'center')
            elif self.labelside == 'left' or self.labelside == 'l':
                fig.text(left - self.labelspacing, bottom + height / 2., self.label,
                         fontsize = self.labelfontsize, fontweight = self.labelfontweight,
                         horizontalalignment = 'right', verticalalignment = 'center')
            else:
                raise Exception('Unacceptable labelside %s' % self.labelside)


        return

    def plot_legend(self, (left, bottom, width), fig, patchheight=0.015,
                    continuous_height=0.1, title_fontsize=7):

        if self.continuous:
            continuous_array_colorshow_legend(
                self.norm, self.colormap, (left, bottom, width, continuous_height),
                fig, title=self.label)
            height = continuous_height
        else:
            height = len(self.label2color)*patchheight
            colorshow_legend(self.label2color, (left, bottom, width, height),
                             fig, title=self.label, title_fontsize=title_fontsize)

        self.leghei = height
        return

class WCbarGroup(object):

    def __init__(self, title=None):

        """ label2color is set at plot time """
        self.cbars = []
        self.label2color = None
        self.uid2cbar = {}
        self.title = title

    def getcat(self, catid):
        g2 = WCbarGroup(title=self.title)
        g2.cbars = [c for c in self.cbars if c.cat==catid]
        g2.label2color = self.label2color
        return g2

    def __len__(self):
        
        return len(self.cbars)
    
    def add_cbar(self, data, uid=None, cat=1, **kwargs):

        cbar = WCbar(data, **kwargs)
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

    def plot_legend(self, (left, bottom, width), fig, patchheight=0.015,
                    continuous_height=0.1, title_fontsize=7):

        # if self.continuous:
        #     continuous_array_colorshow_legend(
        #         self.norm, self.colormap, (left, bottom, width, continuous_height),
        #         fig, title=self.label)
        #     height = continuous_height
        # else:
        height = len(self.label2color)*patchheight

        colorshow_legend(self.label2color, (left, bottom, width, height), fig, title=self.title)

        self.leghei = height
        return

class WCbarCollection():

    def __init__(self, colls=[]):

        self.colls = colls

    def __len__(self):

        return sum([len(coll) for coll in self.colls])

    def cbars(self):

        for cb in self.colls:
            if isinstance(cb, WCbarGroup):
                for _cb in cb.cbars:
                    yield _cb
            elif isinstance(cb, WCbar):
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

def normal_tumor_layout(cd_tumor, cd_normal=None, top_cbars=WCbarCollection(), left_cbars=WCbarCollection(),
                        top_normal_cbars=WCbarCollection(),
                        figwid=18, fighei=8, figfile="tmp.png", **kwargs):

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

    nmwid = 0.1 if cd_normal is not None else 0
    nmpad = 0.003 if cd_normal is not None else 0

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
        cbar.plot_legend([leglft, legbtm, legwid], fig)
        leglft += legwid + legpad
    ff = os.path.splitext(figfile)
    legfigfile = ff[0]+'_legend'+ff[1]
    fig.savefig(legfigfile, bbox_inches='tight', dpi=150)
    print 'Saving legend..'
    print "Done."


def normal_tumor_with_contrast_layout(cd_tumor, cd_normal=None, cd_tumor_contrast=None, cd_normal_contrast=None, top_cbars=[], left_cbars=[], left_cbars_contrast=[], top_normal_cbars=[], figwid=18, fighei=10, figfile="tmp.png", **kwargs):

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
    mahei /= 2
    maheipad = 0.001
    maheiuni = mahei + maheipad

    figfile = figfile

    assert(mawid>0)
    assert(mahei>0)

    print 'left dendro contrast'
    cd_tumor_contrast.D_lft.plot(fig, [mar, mar, ldwid, mahei], orientation='left')

    print 'left color bars contrast'
    lclft = mar+lduni
    for i, cbar in enumerate(left_cbars_contrast.cbars()):
        cbar.plot([lclft+i*lcuni, mar, lcwid, mahei], fig, orientation='v')

    print 'normal heatmap contrast'
    nmlft = lclft + len(left_cbars)*lcuni
    plot_heatmap(cd_normal_contrast.df, [nmlft, mar, nmwid, mahei], fig, interpolation='none')

    print 'tumor heatmap contrast'
    malft = nmlft + nmuni
    plot_heatmap(cd_tumor_contrast.df, [malft, mar, mawid, mahei], fig, interpolation='none')
    
    print "left dendro"
    ldbtm = mar+maheiuni
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
    tcbtm = ldbtm + maheiuni
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
        cbar.plot_legend([leglft, legbtm, legwid], fig)
        leglft += legwid + legpad
    ff = os.path.splitext(figfile)
    legfigfile = ff[0]+'_legend'+ff[1]
    fig.savefig(legfigfile, bbox_inches='tight', dpi=150)
    print 'Saving legend..'
    print "Done."
