
# this is a rewrite of phylogram plotting using ETE-3
# some basic functionality can be achieved by native code
# of ETE-3
from ete3 import Tree
import numpy as np
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
             ):

        PLG_tree_calc_dimension(root, True)

        self.root = root   # ete TreeNode object
        self.root._angle_beg = angle_beg
        self.angle_beg = angle_beg
        self.angle_end = angle_end
        self.angle_span = angle_end - angle_beg

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
        self.node2color_f = None  # map of edge color function
        self.internal_drawer_f = None  # draw on internal node
        self.leaf_drawer_f = None     # draw on leaf node
        self.leafratio = 3 # leaf branch length ratio w.r.t internal branches
        self.blv = 1       # unit branch length radial
        self.lw = 1        # line width
        self.alpha = 0.7   # alpha
        self.fontsize = 12 # font size
        # self.radius = 100  # plot radius

    def drawCircle(self, t, r, circle_radius, **kargs):
        """ r is the radius of the center in polar coordiantes, not the radius of the circle """
        self.ax.add_artist(plt.Circle(
            polar2cart(t._angle_mid, r), circle_radius, transform=self.ax.transData._b, **kargs))


    def drawInternalText(self, t, text, r=None):

        if r is None:
            r = t._r - 0.7 * self.blv
            
        text_angle = PLG_normalize_angle2(t._angle_mid / np.pi * 180) + 90
        self.ax.text(t._angle_mid, r, text,
                     fontsize=12, rotation=text_angle, ha='center', va='center')
        
    def drawLeafText(self, t, r=None):

        if r is None:
            r = t._r + self.blv * self.leafratio

        # extend_adjustment
        r += len(t.common_name)*0.2
        
        if self.radialtext:
            rot = PLG_normalize_angle(t._angle_mid / np.pi * 180)
            ha, va = PLG_calc_text_aln(t._angle_mid / np.pi * 180)
        else:
            rot = 0
            ha, va = PLG_calc_text_aln(t._angle_mid / np.pi * 180)

        # currently text are all center-adjusted
        # this is a hard adjustment and may be inaccturate due to
        # un-equal space of different characters
        # Is there a better solution?
        self.ax.text(
            t._angle_mid, r,
            t.common_name, rotation=rot, ha='center', va='center', fontsize=self.fontsize)

        return rot, ha,va
                
    def drawNode(self, t, ec='k'):

        """ 
                |
        |------ t -------|
        |                |
        c                c
        """
        # set t._r, t._angle_end, t._angle_mid based on t._angle_beg
        t._r = t.depth*self.blv
        t._angle_end = t._angle_beg + t._angle_span
        t._angle_mid = (t._angle_beg + t._angle_end) / 2.0

        if t.is_leaf():
            if self.leaf_drawer_f is not None:
                self.leaf_drawer_f(self, t)
            else:
                self.drawLeafText(t)
        elif self.internal_drawer_f is not None:
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
            if self.node2color_f is not None:
                success, _ec = self.node2color_f(c)
                if success:
                    ec2 = _ec
                else:
                    ec2 = ec
            else:
                ec2 = ec

            # update tangent angles
            tangent_angle_beg = min(c_angle_mid, tangent_angle_beg)
            tangent_angle_end = max(c_angle_mid, tangent_angle_end)

            # plot radial line
            c_r = (c.depth-1+self.leafratio)*self.blv if c.is_leaf() else c.depth*self.blv
            self.ax.plot([c_angle_mid, c_angle_mid], # angle
                         [t._r, c_r],                 # axial
                         color=ec2, lw=self.lw, alpha=self.alpha)

            c._angle_beg = ab
            self.drawNode(c, ec=ec2)
            ab += c._angle_span

        # plot horizontal line
        self.ax.plot(np.linspace(tangent_angle_beg, tangent_angle_end, 40),
                     [t._r]*40, color=ec, lw=self.lw, alpha=self.alpha)
        

    def draw(self, fig=None):
        
        """ draw TreeNode t """
        if fig is None:
            self.fig = plt.figure(figsize=(13,13))
        else:
            self.fig = fig
        self.ax = self.fig.add_axes([0,0,1,1], projection='polar')
        self.ax.set_axis_off()
        self.drawNode(self.root)

