import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.collections as mcolls

def plot_heatmap(fig, dim, data, cmap=None):
    if cmap is None:
        cmap = cm.jet
        cmap.set_bad('w', 1)
    ax = fig.add_axes(dim, frameon=False)
    ax.imshow(data, aspect='auto', origin='lower', cmap=cmap)
    ax.set_xticks([])
    ax.set_yticks([])

    return ax

def colorshow_legend(fig, (left,bottom,width), l2c, fontsize=8, patchheight=0.015, horizontalspace=0.05):

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

def discrete_array_colorshow(fig, dim, data, orientation='horizontal'):

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
    colors, level2color = wzcolors.map_distinct_colors_hex(data, other2grey=True)
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
