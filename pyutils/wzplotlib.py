import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as mcolls

# def plot_heatmap(fig, dim, data, cmap=None):
#     if cmap is None:
#         cmap = cm.jet
#         cmap.set_bad('w', 1)
#     ax = fig.add_axes(dim, frameon=False)
#     ax.imshow(data, aspect='auto', origin='lower', cmap=cmap)
#     ax.set_xticks([])
#     ax.set_yticks([])

#     return ax

# def colorshow_legend(fig, (left,bottom,width), l2c, fontsize=8, patchheight=0.015, horizontalspace=0.05):

#     n = len(l2c)
#     ax = fig.add_axes((left, bottom, width, n*patchheight), frameon=False)
#     ax.set_xlim(0,1)
#     ax.set_ylim(0,n)
#     levels, colors = zip(*l2c.items())
#     for i, color in enumerate(colors):
#         ax.add_patch(mpatches.Rectangle((0,i),1,1,facecolor=color,edgecolor='white'))
#         ax.text(1+horizontalspace,i+0.5,str(levels[i]),verticalalignment="center", fontsize=fontsize, fontweight='light')
#     # ax.set_yticks([i+0.5 for i in xrange(n)])
#     # ax.set_yticklabels(levels, fontsize=fontsize, fontweight='light')
#     # ax.yaxis.tick_right()
#     ax.set_yticks([])
#     ax.set_xticks([])

#     return ax

# def discrete_array_colorshow(fig, dim, data, orientation='horizontal'):

#     n = len(data)
#     if orientation == 'horizontal':
#         beg = '(i,0)'
#         xmax = n
#         ymax = 1
#     elif orientation == 'vertical':
#         beg = '(0,i)'
#         xmax = 1
#         ymax = n
#     else:
#         raise Exception("unknown orientation")

#     import wzcolors
#     ax = fig.add_axes(dim, frameon=False)
#     colors, level2color = wzcolors.map_distinct_colors_hex(data, other2grey=True)
#     color2patches = {}
#     for i, color in enumerate(colors):
#         if color not in color2patches:
#             color2patches[color] = []
#         color2patches[color].append(mpatches.Rectangle(eval(beg),1,1))

#     for color, patches in color2patches.iteritems():
#         ax.add_collection(mcolls.PatchCollection(patches, facecolor=color, edgecolor='none'))

#     ax.set_xlim(0,xmax)
#     ax.set_ylim(0,ymax)
#     ax.set_xticks([])
#     ax.set_yticks([])

#     return ax, level2color

def plot_heatmap(data, dim=[0.1,0.1,0.85,0.85], fig=None, xlabels=None, ylabels=None, label_fontsize=8, cmap=None):
    if cmap is None:
        cmap = cm.jet
        cmap.set_bad('w', 1)

        if fig is None:
            fig = plt.figure()
            
    ax = fig.add_axes(dim, frameon=False)
    ax.imshow(data, aspect='auto', origin='lower', cmap=cmap)

    ax.set_xticks([])
    ax.set_yticks([])
    if xlabels is not None:
        if type(xlabels) == bool:
            xlabels = data.columns.format()
        for i in xrange(data.shape[1]):
            ax.text(i+0.2, -1, xlabels[i], rotation=90, horizontalalignment='center', verticalalignment='top', fontsize=label_fontsize)

    return ax

def colorshow_legend(l2c, (left,bottom,width), fig, fontsize=8, patchheight=0.015, horizontalspace=0.05):

    n = len(l2c)
    ax = fig.add_axes((left, bottom, width, n*patchheight), frameon=False)
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

    return ax

def discrete_array_colorshow(data, dim, fig, orientation='horizontal',
                             greyscale=False, greyscale_range=(0.1,0.9)):

    n = len(data)
    if orientation == 'horizontal':
        beg = '(i,0)'
        xmax = n
        ymax = 1
    elif orientation == 'vertical':
        beg = '(0,i)'
        xmax = 1
        ymax = n
    else:
        raise Exception("unknown orientation")

    import wzcolors
    ax = fig.add_axes(dim, frameon=False)
    colors, level2color = wzcolors.map_distinct_colors_hex(
        data, other2grey=True, greyscale=greyscale, greyscale_range=greyscale_range)
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
