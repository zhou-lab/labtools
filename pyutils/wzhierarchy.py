import numpy as np
import itertools

"""
usage:
ez_cluster(df)
ez_cluster_row(df)
ez_cluster_column(df)

"""

def p_euclidean(X):
    from numpy.linalg import norm
    m, n = X.shape
    dm = np.zeros((m * (m - 1)) // 2, dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i + 1, m):
            diff = X[i] - X[j]
            dm[k] = norm(diff[np.invert(np.isnan(diff))])
            k = k + 1
    return dm

def compute_dendrogram(Z, w=None):

    """ w is of size n is the initial weight of the leaves
    internal weights are calculated and stored into weight
    """
    n = Z.shape[0] + 1

    # post-order traversal, set weight
    weight = np.zeros(2*n-1, dtype='float32')
    cladesize = np.zeros(2*n-1, dtype='float32')
    if w is not None:
        stack = [2*n-2]
        last = None
        while stack:
            i = stack[-1]
            if i < n or Z[i-n,0] == last or Z[i-n,1] == last:
                stack.pop()
                last = i
                if i < n:       # leaf
                    cladesize[i] = 1
                    weight[i] = w[i]
                else:
                    cladesize[i] = cladesize[Z[i-n,0]] + cladesize[Z[i-n,1]]
                    weight[i] = weight[Z[i-n,0]] + weight[Z[i-n,1]]
                    # weight[i] = max(weight[Z[i-n,0]], weight[Z[i-n,1]])
                    # weight[i] = (weight[Z[i-n,0]] + weight[Z[i-n,1]]) / 2.0
                    # print Z[i-n,0], Z[i-n,1], i, weight[Z[i-n,0]], weight[Z[i-n,1]], weight[i]
            else:
                stack.append(Z[i-n,1])
                stack.append(Z[i-n,0]) # visit this first

    # post-order traversal, set order
    leaforder = []
    stack = [2*n-2]
    last = None
    heights = np.zeros(2*n-1, dtype='float32')
    locs = np.zeros(2*n-1, dtype='int')
    lefts = np.zeros(2*n-1, dtype='int')
    rights = np.zeros(2*n-1, dtype='int')
    while stack:
        i = stack[-1]
        if i < n:               # leaf
            stack.pop()
            last = i
            leaforder.append(i)
            heights[i] = 0
            locs[i] = len(leaforder)
            lefts[i] = len(leaforder)
            rights[i] = len(leaforder)
        elif Z[i-n,0] == last or Z[i-n,1] == last: # nonleaf
            stack.pop()
            last = i
            locs[i] = (locs[Z[i-n,0]]+locs[Z[i-n,1]]) / 2.0
            heights[i] = (Z[i-n,2]+heights[Z[i-n,0]]+heights[Z[i-n,1]])/2.0
            if heights[i] < heights[Z[i-n,0]]:
                heights[i] = heights[Z[i-n,0]]*1.05
            if heights[i] < heights[Z[i-n,1]]:
                heights[i] = heights[Z[i-n,1]]*1.05
            lefts[i] = min(lefts[int(Z[i-n,0])], lefts[int(Z[i-n,1])])
            rights[i] = max(rights[int(Z[i-n,0])], rights[int(Z[i-n,1])])
        else:
            # bigger weight visited last
            if weight[Z[i-n,1]]/cladesize[Z[i-n,1]] < weight[Z[i-n,0]]/cladesize[Z[i-n,0]]:
                stack.append(Z[i-n,0])
                stack.append(Z[i-n,1])
            else:
                stack.append(Z[i-n,1])
                stack.append(Z[i-n,0])

    d = Dendrogram()
    d.Z = Z
    d.locs = locs
    d.lefts = lefts
    d.rights = rights
    d.heights = heights
    d.leaforder = leaforder

    return d

class Dendrogram():

    def __init__(self):

        pass

    def plot(self, fig, dim=[0.1,0.1,0.8,0.7], color_threshold=None, orientation="top"):

        Z, heights, leaforder = self.Z, self.heights, self.leaforder
        locs, lefts, rights = self.locs, self.lefts, self.rights
        n = Z.shape[0] + 1

        # identify color information
        color = ['k']*(2*n-1)
        if color_threshold is None:
            color_threshold = 0.7*max(Z[:,2])

        import wzcolors
        clade_colors = itertools.cycle(wzcolors.get_spectral_colors_rgb_dark(6))
        curr_col = next(clade_colors)
        col_changed = False
        stack = [2*n-2]
        last = None
        while stack:
            i = stack[-1]
            if color[int(i)] is None and not col_changed:
                curr_col = next(clade_colors)
                col_changed = True
            if i < n or Z[i-n,0] == last or Z[i-n,1] == last:
                stack.pop()
                last = i
                if i < n:
                    broken = False
                    color[int(i)] = curr_col
                    col_changed = False
                else:
                    if (color[int(Z[i-n, 0])] != 'k' and
                        color[int(Z[i-n, 1])] != 'k' and
                        color[int(Z[i-n, 0])] == color[int(Z[i-n, 1])]):
                        color[int(i)] = color[int(Z[i-n, 0])]
                    else:
                        color[int(i)] = 'k'
            else:
                if Z[i-n,2] > color_threshold:
                    color[int(Z[i-n, 1])] = None
                    color[int(Z[i-n, 0])] = None
                stack.append(Z[i-n,1])
                stack.append(Z[i-n,0]) # visit this first

        import matplotlib as mpl
        import matplotlib.patches as mpatches
        axis = fig.add_axes(dim, frameon=False)

        xlines = []
        ylines = []
        clines = []
        for i in xrange(2*n-1):
            if i<n: continue
            col = color[i]
            c1 = int(Z[i-n,0])
            c2 = int(Z[i-n,1])
            # horizontal
            xlines.append((locs[c1], locs[c2]))
            ylines.append((heights[i], heights[i]))
            clines.append(col)
            # vertical of c0
            xlines.append((locs[c1], locs[c1]))
            ylines.append((heights[i], heights[c1]))
            clines.append(col)
            # vertical of c1
            xlines.append((locs[c2], locs[c2]))
            ylines.append((heights[i], heights[c2]))
            clines.append(col)
                
        if orientation == 'top':
            axis.set_ylim(0, max(heights))
            axis.set_xlim(0.5, max(locs)+0.5)
        elif orientation == 'bottom':
            axis.set_ylim(max(heights), 0)
            axis.set_xlim(0.5, max(locs)+0.5)
        elif orientation == 'left':
            tmp = xlines
            xlines = ylines
            ylines = tmp
            axis.set_xlim(max(heights), 0)
            axis.set_ylim(0.5, max(locs)+0.5)
        elif orientation == 'right':
            tmp = xlines
            xlines = ylines
            ylines = tmp
            axis.set_xlim(0, max(heights))
            axis.set_ylim(0.5, max(locs)+0.5)
        else:
            raise Exception("Invalid orientation %s" % orientation)

        col2lines = {}
        for xline, yline, col in zip(xlines, ylines, clines):
            if i >= n:
                if col not in col2lines:
                    col2lines[col] = []
                col2lines[col].append(list(zip(xline, yline)))
        for col, lines in col2lines.iteritems():
            axis.add_collection(mpl.collections.LineCollection(lines,color=col))

        axis.set_xticks([])
        axis.set_yticks([])

        return axis

class ClusterData():

    def __init__(self):

        self.D_top = None
        self.D_lft = None

    def lft_order(self):

        if self.D_lft is None:
            return None
        else:
            return self.D_lft.leaforder

    def top_order(self):

        if self.D_top is None:
            return None
        else:
            return self.D_top.leaforder

def ez_good_row(df, decision='all'):

    if decision == 'all':
        return df.apply(lambda x: all(x.notnull()), axis=1)
    elif decision == 'any':
        return df.apply(lambda x: any(x.notnull()), axis=1)


def pdist_na(df, metric):

    if metric == 'euclidean':
        p = p_euclidean(df.transpose().as_matrix())
    elif metric == 'hamming':
        p = p_hamming(df.transpose().as_matrix())
    else:
        raise Exception('unknown metric %s' % metric)

    return p
    
def ez_cluster(df, fast=False, good_row=None, categorical=False, metric='euclidean'):

    """ with fast=True option, only consider probes with notnull in all samples """

    import fastcluster

    if fast and good_row is None:
        good_row = ez_good_row(df, decision='all')
    
    if good_row is not None:
        df = df.loc[good_row,]

    if categorical:
        if metric=='euclidean':
            metric='hamming'
        vals = sorted(list(set(df.values.flatten())))
        inds = range(len(vals))
        val2ind = dict(zip(vals, inds))
        _df = df.applymap(lambda x: val2ind[x])
    else:
        _df = df

    d = ClusterData()
    if fast:
        d.Z_top = fastcluster.linkage(_df.transpose().as_matrix(), metric=metric, method='ward', preserve_input=True)
        d.Z_lft = fastcluster.linkage(_df.as_matrix(), method='ward', metric=metric, preserve_input=True)
    else:
        p_top = pdist_na(_df, metric)
        d.Z_top = fastcluster.linkage(p_top, method='ward', preserve_input=True)
        p_lft = pdist_na(_df.transpose(), metric)
        d.Z_lft = fastcluster.linkage(p_lft, method='ward', preserve_input=True)

    if categorical:
        d.w_lft = None
    else:
        d.w_lft = np.apply_along_axis(np.nansum, 1, _df)
    d.D_lft = compute_dendrogram(d.Z_lft, d.w_lft)
    if categorical:
        d.w_top = None
    else:
        d.w_top = np.apply_along_axis(np.nansum, 0, _df)
    # print d.w_top.shape, _df.shape, d.Z_top.shape, d.Z_lft.shape
    d.D_top = compute_dendrogram(d.Z_top, d.w_top)
    d.df = df.iloc[d.D_lft.leaforder,d.D_top.leaforder]

    return d

def ez_cluster_column(df, fast=False, good_row=None):

    import fastcluster
    d = ClusterData()
    d.df = df
    if good_row is not None:
        d.df = d.df.loc[good_row,]
    if fast:
        d.Z_top = fastcluster.linkage(d.df.transpose().as_matrix(), method='ward', preserve_input=True)
    else:
        d.Z_top = fastcluster.linkage(p_euclidean(d.df.transpose().as_matrix()),
                                      method='ward', preserve_input=True)
    d.w_top = np.apply_along_axis(np.nansum, 0, d.df)
    d.D_top = compute_dendrogram(d.Z_top, d.w_top)
    d.df = d.df.iloc[:,d.D_top.leaforder]

    return d


def ez_cluster_row(df, fast=False):

    import fastcluster
    d = ClusterData()
    d.df = df
    if fast:
        d.Z_lft = fastcluster.linkage(d.df.as_matrix(), method='ward', preserve_input=True)
    else:
        d.Z_lft = fastcluster.linkage(p_euclidean(d.df.as_matrix()),
                                      method='ward', preserve_input=True)
    d.w_lft = np.apply_along_axis(np.nansum, 1, d.df)
    d.D_lft = compute_dendrogram(d.Z_lft, d.w_lft)
    d.df = d.df.iloc[d.D_lft.leaforder,:]

    return d

