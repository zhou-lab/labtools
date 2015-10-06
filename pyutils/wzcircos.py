from wzcore import *
import numpy as np
import re
import faidx
import matplotlib.patches as mpatches
import matplotlib.path as mpath

def wzpolygon(angle_beg, angle_end, inner_radius, outer_radius, angle_inc=0.01, fc='blue', alpha=0.8, edgecolor='none'):
    angle_num = max(3,int(float(angle_end - angle_beg) / angle_inc))
    angle_num = min(500, angle_num)
    p = mpatches.Polygon(zip([angle_beg, angle_beg]+
                             list(np.linspace(angle_beg, angle_end, angle_num))+
                             [angle_end,angle_end]+
                             list(np.linspace(angle_end, angle_beg, angle_num)), 
                             [inner_radius, outer_radius]+
                             [outer_radius]*angle_num+
                             [outer_radius, inner_radius]+
                             [inner_radius]*angle_num), alpha=alpha, fc=fc, edgecolor=edgecolor)
    return p

def getsortedchrms(genome):
    special = []
    chrms = []
    specialord = {'X':1,'Y':2,'M':3}
    for chrm in genome.faidx.keys():
        m = re.match(r'chr([0-9XY]*)$', chrm)
        if m:
            if m.group(1) in 'XYM':
                special.append((specialord[m.group(1)], chrm))
            else:
                chrms.append((int(m.group(1)), chrm))
                
    chrms = sorted(chrms)+sorted(special)
    return [(_[1], genome.faidx[_[1]][0]) for _ in chrms]

def normalize_text_angle(angle):

    while angle > 90:
        angle -= 180

    while angle < -90:
        angle += 180

    return angle

class CircosLayout:

    def __init__(self, genome, angle_beg=0, angle_end=np.pi*2,
                 angle_chrm_space=0.01, inner_radius=45, outer_radius=50, angle_inc=0.01,
                 bezier_anchor = 4.0, # higher the number the higher the peak
    ):
        
        chrms = getsortedchrms(genome)
        self.totalbases = sum([_[1] for _ in chrms])
        self.angle_per_base = (angle_end - angle_beg - angle_chrm_space*len(chrms)) / float(self.totalbases)
        self.angle_chrm_space = angle_chrm_space
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.chrms = [_[0] for _ in chrms]
        self.chrm2len = dict(chrms)
        self.chrm2angles = {}
        self.angle_inc = angle_inc
        self.bezier_anchor = bezier_anchor
        tmp = angle_beg
        self.connect_radius = inner_radius*0.95
        for chrm in chrms:
            tmp1 = tmp + chrm[1]*self.angle_per_base
            self.chrm2angles[chrm[0]] = (tmp, tmp1)
            tmp = tmp1
            tmp += angle_chrm_space

        return

    def loc2angle(self, chrm, pos):

        return self.chrm2angles[chrm][0] + self.angle_per_base*pos
    
    def plot_chromosomes(self, ax, fc='blue'):

        for chrm in self.chrms:
            angles = self.chrm2angles[chrm]
            ax.add_patch(wzpolygon(angles[0], angles[1], self.inner_radius, self.outer_radius, fc=fc))
            text_angle = (angles[0]+angles[1])/2.0
            ax.text(text_angle, self.outer_radius*1.05, chrm, ha='center', va='center', rotation=normalize_text_angle(text_angle/(np.pi*2)*360))

    def _plot_arc(self, angle1, angle2, radius):
        angle_num = max(3,int(float(angle2 - angle1) / self.angle_inc))
        angle_num = min(500, angle_num)
        for angle in np.linspace(angle1, angle2, angle_num):
            self.verts.append((angle, radius))
            self.codes.append(mpath.Path.LINETO)

    def _plot_bezier(self, angle1, angle2, radius):

        while abs(angle2 - angle1) > np.pi:
            if angle1 < angle2:
                angle1 += np.pi*2
            else:
                angle2 += np.pi*2
        anchor_radius = radius*np.cos(abs(angle2 - angle1)/2.0)/self.bezier_anchor
        anchor_angle = (angle1+angle2)/2.0
        self.verts.extend([(angle1, radius), (anchor_angle, anchor_radius), (angle2, radius)])
        if self.codes:
            a = mpath.Path.LINETO
        else:
            a = mpath.Path.MOVETO
        self.codes.extend([a, mpath.Path.CURVE3, mpath.Path.CURVE3])


    def _init_path(self):

        self.verts = []
        self.codes = []

    def _close_path(self, theta, r):

        self.verts.append((theta, r))
        self.codes.append(mpath.Path.STOP)

    def path(self):
        return mpath.Path(self.verts, self.codes)

    def plot_enclosure(
            self, ax, _angle1, _angle2, fc='grey',
            
    ):
        
        angle1 = min(_angle1, _angle2)
        angle2 = max(_angle1, _angle2)

        # connect
        self._init_path()
        self._plot_bezier(angle1, angle2, self.connect_radius)

        # arc
        self._plot_arc(angle1, angle2, self.connect_radius)
        self._close_path(angle1, self.connect_radius)
        
        ax.add_patch(mpatches.PathPatch(self.path(), alpha=0.4, fc=fc, ec='none'))

    def plot_connection(
            self, ax, _source_angle1, _source_angle2, _target_angle1, _target_angle2, fc='grey'
    ):

        source_angle1 = min(_source_angle1, _source_angle2)
        source_angle2 = max(_source_angle1, _source_angle2)
        target_angle1 = min(_target_angle1, _target_angle2)
        target_angle2 = max(_target_angle1, _target_angle2)

        # source 2 to target 1
        self._init_path()
        self._plot_bezier(source_angle2, target_angle1, self.connect_radius)
        self._plot_arc(target_angle1, target_angle2, self.connect_radius)
        self._plot_bezier(target_angle2, source_angle1, self.connect_radius)
        self._plot_arc(source_angle1, source_angle2, self.connect_radius)
        self._close_path(source_angle2, self.connect_radius)

        ax.add_patch(mpatches.PathPatch(self.path(), alpha=0.4, fc=fc, ec='none'))

    def plot_connection_chrm(
            self, ax,
            _source_chrm, _source_beg, _source_end,
            _target_chrm, _target_beg, _target_end, fc='grey'):

        self.plot_connection(ax, self.loc2angle(_source_chrm, _source_beg),
                             self.loc2angle(_source_chrm, _source_end),
                             self.loc2angle(_target_chrm, _target_beg),
                             self.loc2angle(_target_chrm, _target_end),)

    def plot_tracks(self, ax, track, track_inner_radius, track_outer_radius, direction='outer'):

        """ track is [(chrm, loc, value)] """

        vals = [v for c,l,v in track]
        minval = np.min(vals)
        maxval = np.max(vals)
        valrange = maxval - minval
        track_height = track_outer_radius - track_inner_radius
        for chrm, loc, value in tracks:
            angle = self.loc2angle(chrm, loc)
            if direction == 'outer':
                ax.plot([angle, angle],[track_inner_radius+(value-minval)/valrange*track_height])
            else:
                ax.plot([angle, angle],[track_inner_radius+(1.0-(value-minval)/valrange)*track_height])

        return
