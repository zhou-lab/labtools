
# this is a rewrite of phylogram plotting using ETE-3
# some basic functionality can be achieved by native code
# of ETE-3
from ete3 import Tree
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import matplotlib

def PLG_tree_calc_dimension(t, isroot=False):
    
    """ the distance from current to its deepest leaf """
    if isroot:
        t.depth = 0
    if t.is_leaf():
        t.height = 0
        t.nleaves = 1
        return
    t.height = 1
    t.nleaves = 0
    for c in t.children:
        c.depth = t.depth + 1
        PLG_tree_calc_dimension(c)
        if c.height + 1 > t.height:
            t.height = c.height + 1
        t.nleaves += c.nleaves

def PLG_normalize_angle(angle):
    
    if angle % 360 > 90 and angle % 360 < 270:
        angle += 180

    return angle


def PLG_normalize_angle2(angle):
    
    if angle % 360 > 0 and angle % 360 < 180:
        angle += 180

    return angle


def PLG_calc_text_aln(angle, hc=0, vc=40):
    ra = angle % 360
    if ra > 90+hc and ra < 270-hc:
        ha = 'right'
    elif ra < 90-hc or ra > 270+hc:
        ha = 'left'
    else:
        ha = 'center'
        
    if ra > 0 and ra < 180:
        va = 'bottom'
    elif ra > 180+vc and ra < 360-vc:
        va = 'top'
    else:
        va = 'baseline'

    return ha, va

def polar2cart(theta, r):
    """ polar coordinate to cartesian coordinate """
    return (r*np.cos(theta), r*np.sin(theta));

class PhyloGram():

    def __init__(self, root, angle_beg, angle_end,
                 target_clade=None, # a special clade to zoom in
                 target_clade_angle_span=np.pi, # set when target_clade is not None
                 clade_colors = pd.Series(),
                 phylo_contraction = 1,
                 leaffix = False, # whether all leaves are of the same depth
                 fig = None
             ):

        PLG_tree_calc_dimension(root, True)

        self.root = root   # ete TreeNode object
        self.root._angle_beg = angle_beg # set the _angle_beg for the root
        self.angle_beg = angle_beg
        self.angle_end = angle_end
        self.angle_span = angle_end - angle_beg
        self.clade_colors = clade_colors
        self.leaffix = leaffix

        # reset all node angle span
        for n in root.traverse():
            n._angle_span = -1

        # set target clade angle span
        if target_clade is None:
            target_clade_angle_span = 0
            n_clade_leaves = 0
        else:
            n_clade_leaves = len(target_clade.get_leaves())
            for n in target_clade.get_leaves():
                n._angle_span = target_clade_angle_span / n_clade_leaves

        # set non-target-clade angle
        nontarget_angle = ((self.angle_span - target_clade_angle_span) /
                            (len(root.get_leaves()) - n_clade_leaves))
        for n in root.get_leaves():
            if n._angle_span < 0:
                n._angle_span = nontarget_angle

        # set angle span for internal nodes
        for n in root.traverse('postorder'):
            if not n.is_leaf():
                n._angle_span = sum([c._angle_span for c in n.children])

        # aesthetics
        self.radialtext = True
        self.innerlabels = []     # internal branch label
        self.internal_drawer_f = None  # draw on internal node
        self.leaf_drawer_f = None     # draw on leaf node
        self.leafratio = 3 # leaf branch length ratio w.r.t internal branches
        self.blv = 1 * phylo_contraction # unit branch length radial
        self.lw = 1        # line width
        self.alpha = 0.7   # alpha
        self.fontsize = 12 # font size
        # self.radius = 100  # plot radius

        """ draw TreeNode t """
        if fig is None:
            self.fig = plt.figure(figsize=(13,13))
        else:
            self.fig = fig
        self.ax = self.fig.add_axes([0,0,1,1], projection='polar')
        self.ax.set_axis_off()
        self.drawNode(self.root)

    def drawCircle(self, t, r, circle_radius, **kargs):
        """ r is the radius of the center in polar coordiantes, not the radius of the circle """
        self.ax.add_artist(plt.Circle(
            polar2cart(t._angle_mid, r), circle_radius, transform=self.ax.transData._b, **kargs))

    def drawNode(self, t, ec='k'):

        """ 
                |
        |------ t -------|
        |                |
        c                c

        The parent of t sets t._angle_beg
        This sets t._angle_end, t._angle_mid, t._r
        """
        # set t._r, t._angle_end, t._angle_mid based on t._angle_beg
        if t.is_leaf() and self.leaffix:
            t._r = self.root.height*self.blv
        else:
            t._r = t.depth*self.blv

        # initialize t._r - where next track will be added
        if t.is_leaf():
            t._r2 = t._r
            
        t._angle_end = t._angle_beg + t._angle_span
        t._angle_mid = (t._angle_beg + t._angle_end) / 2.0

        if t.is_leaf():
            if self.leaf_drawer_f is not None:
                self.leaf_drawer_f(self, t)
            # else:
            #     self.drawLeafText(t)
        else:
            if self.internal_drawer_f is not None:
                self.internal_drawer_f(self, t)
        
        # plot branches to children
        tangent_angle_beg = t._angle_mid
        tangent_angle_end = t._angle_mid
        children = sorted(t.get_children(), key=lambda c: c.nleaves)
        # reverse=random.random() < 0.5)
        
        ab = t._angle_beg        # angle begin for each children
        for c in children:
            c_angle_mid = ab + c._angle_span / 2.

            # whether to change edge color?
            if c.name in self.clade_colors.index:
                ec2 = self.clade_colors[c.name]
            else:
                ec2 = ec
            # if self.clade_colors is not None:
            #     success, _ec = self.node2color_f(c)
            #     if success:
            #         ec2 = _ec
            #     else:
            #         ec2 = ec
            # else:
            #     ec2 = ec

            # update tangent angles
            tangent_angle_beg = min(c_angle_mid, tangent_angle_beg)
            tangent_angle_end = max(c_angle_mid, tangent_angle_end)

            # plot radial line
            if c.is_leaf():
                if self.leaffix:
                    c_r = self.root.height*self.blv
                else:
                    c_r = (c.depth-1+self.leafratio)*self.blv
            else:
                c_r = c.depth*self.blv
            self.ax.plot([c_angle_mid, c_angle_mid], # angle
                         [t._r, c_r],                # axial
                         color=ec2, lw=self.lw, alpha=self.alpha)

            c._angle_beg = ab
            self.drawNode(c, ec=ec2)
            ab += c._angle_span

        # plot horizontal line
        self.ax.plot(np.linspace(tangent_angle_beg, tangent_angle_end, 40),
                     [t._r]*40, color=ec, lw=self.lw, alpha=self.alpha)
        

    def add_leaf_bar_track(self, srs, inner_space_ratio = 0.1, track_color='k', baseline_color='k', track_height = None, tick_length=0.03, y_max = None):

        """
        srs is a pandas Series, indexed by t.name 
        inner_space is default to 0.5*track_height
        """

        if track_height is None:
            track_height = self.blv
            
        inner_space = inner_space_ratio * track_height

        if y_max is None:
            y_max = srs.max() # the data maximum

        # import pdb; pdb.set_trace()
        for t in self.root.get_leaves():
            t._r2 += inner_space + track_height
            y = srs[[t.name]].iloc[0] / y_max * track_height
            if (isinstance(track_color, pd.Series)):
                tc = track_color[[t.name]].iloc[0]
            else:
                tc = 'k'

            # import pdb; pdb.set_trace()
            self.ax.bar(t._angle_mid, height = y, width = 0.04, bottom = t._r2 - y, facecolor=tc, alpha=0.85)

        # the base line
        self.ax.plot(np.linspace(self.angle_beg, self.angle_end, 100), [t._r2]*100, color=baseline_color, lw=0.5)

        # y-axis
        self.ax.plot([self.angle_beg, self.angle_beg], [t._r2, t._r2-track_height], color=baseline_color, lw=0.5)

        # y-axis ticks
        self.ax.plot(np.linspace(self.angle_beg, self.angle_beg-tick_length, 3), [t._r2]*3, color=baseline_color, lw=1)
        self.ax.plot(np.linspace(
            self.angle_beg, self.angle_beg-tick_length, 3), [t._r2-track_height]*3, color=baseline_color, lw=1)
        self.ax.text(self.angle_beg, t._r2 - track_height, "%1.2g" % y_max,
                     rotation=self.angle_beg/np.pi*180+90, verticalalignment='bottom', horizontalalignment='right')

    def add_leaf_bubble_track(self, srs, inner_space = None, track_color='k', track_height = None, max_radius = 9999):

        """
        srs is a pandas Series, indexed by t.name
        inner_space is default to 0.5*track_height
        """

        if track_height is None:
            track_height = self.blv
            
        if inner_space is None:
            inner_space = 0.5 * track_height
            
        y_max = srs.max() # the data maximum
        for t in self.root.get_leaves():
            t._r2 += inner_space + track_height
            circle_radius = srs[t.name] / y_max * min(track_height/2, max_radius)
            self.ax.add_artist(plt.Circle(polar2cart(
                t._angle_mid, t._r2-track_height/2), circle_radius, transform=self.ax.transData._b))

    def add_leaf_labels(self, labels=None, r_pad=1):

        for t in self.root.get_leaves():

            if labels is None:
                label = t.common_name
            else:
                label = labels[t.name]
                
            r = t._r2 + r_pad
            # + self.blv * self.leafratio

            if self.radialtext:
                rot = PLG_normalize_angle(t._angle_mid / np.pi * 180)
                ha, va = PLG_calc_text_aln(t._angle_mid / np.pi * 180, hc=0, vc=0)
            else:
                rot = 0
                ha, va = PLG_calc_text_aln(t._angle_mid / np.pi * 180, hc=0, vc=0)

            # if label[:4] == 'Boli':
            # import pdb; pdb.set_trace()

            # currently text are all center-adjusted
            # this is a hard adjustment and may be inaccturate due to
            # un-equal space of different characters
            # Is there a better solution?
            self.ax.text(t._angle_mid, r, label, rotation=rot, ha=ha, va=va, fontsize=self.fontsize)
            # return rot, ha,va

    def add_internal_labels(self, internal_node_labels):

        """ internal_node_labels is a Series mapping node name to label """

        for t in self.root.traverse():
            if t.name in internal_node_labels.index:
                r = t._r - 0.7 * self.blv
                text_angle = PLG_normalize_angle2(t._angle_mid / np.pi * 180) + 90
                self.ax.text(t._angle_mid, r, internal_node_labels[t.name],
                             fontsize=12, rotation=text_angle, ha='center', va='center')


    def add_radial_extension(self):

        for t in self.root.get_leaves():
            self.ax.plot([t._angle_mid, t._angle_mid],[t._r, t._r2], lw=0.5, ls='dashed', color='black')
