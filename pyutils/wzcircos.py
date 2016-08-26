from wzcore import *
import wzcolors
import numpy as np
import re
import faidx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.path as mpath
import matplotlib.collections as mcollections
import pandas as pd

cytoband_color = {'gneg':'k',
                  'gpos25':'k','gpos33':'k','gpos50':'k','gpos66':'k','gpos75':'k','gpos100':'k',
                  'gvar':'k','acen':'k','stalk':'k'}

cytoband_alpha = {'gneg':0,
                  'gpos25':0.25,'gpos33':0.3,'gpos50':0.50,'gpos66':0.6,'gpos75':0.75,'gpos100':1.0,
                  'gvar':0,'acen':0.9,'stalk':1.0}

def wzpolygon(angle_beg, angle_end, inner_radius, outer_radius, angle_inc=0.01, fc='blue', alpha=0.8, ec='none'):
    angle_num = max(3,int(float(angle_end - angle_beg) / angle_inc))
    angle_num = min(500, angle_num)
    p = mpatches.Polygon(zip([angle_beg, angle_beg]+
                             list(np.linspace(angle_beg, angle_end, angle_num))+
                             [angle_end,angle_end]+
                             list(np.linspace(angle_end, angle_beg, angle_num)),
                             [inner_radius, outer_radius]+
                             [outer_radius]*angle_num+
                             [outer_radius, inner_radius]+
                             [inner_radius]*angle_num), alpha=alpha, fc=fc, edgecolor=ec, lw=0.5)
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
    return [(_[1], genome.faidx[_[1]][0]) for _ in chrms] # [(chromosome name, chromosome length)]

def normalize_text_angle(angle, tangent=False):

    while angle > 90:
        angle -= 180
            
    while angle < -90:
        angle += 180

    if tangent:
        if angle < 0:
            angle += 90
        else:
            angle -= 90

    return angle

def round_to_10(n):

    if n>10:
        while n>10:
            n/=10
    return n

class CircosTrack(object):

    def __init__(self, df, layout):
        self.df = df
        self.layout = layout
        self.track_height = None
        self.track_bottom = None
        
    def plot(self):
        pass

def polar2cart(theta, r):
    return (r*np.cos(theta/180.*np.pi), r*np.sin(theta/180.*np.pi));

class CircosTrackBar(CircosTrack):

    def __init__(self, df, layout, **plot_kwargs):

        # summarize value distribution
        CircosTrack.__init__(self, df, layout)
        self.vals = [v for k,(c,beg,end,v) in df.iterrows()]
        self.minval = np.min(self.vals)
        self.maxval = np.max(self.vals)
        self.valrange = float(self.maxval - self.minval)
        self.direction = 'inner'
        self.maxr = 0.9         # scale the actual plot area
        self.plot_kwargs = plot_kwargs
        
    def plot(self, background=False):

        kwargs1 = {}
        kwargs2 = {}
        if 'fc' in self.plot_kwargs:
            kwargs1['fc'] = self.plot_kwargs['fc']

        if 'alpha' in self.plot_kwargs:
            kwargs1['alpha'] = self.plot_kwargs['alpha']
            
        if 'label_fontsize' in self.plot_kwargs:
            kwargs2['fontsize'] = self.plot_kwargs['label_fontsize']
        else:
            kwargs2['fontsize'] = 9

        if 'label' in self.plot_kwargs:
            label = self.plot_kwargs['label']
        else:
            label = self.df.columns[3]

        ly = self.layout

        angle1s = []
        anglewids = []
        heights = []
        for k, (chrm, beg, end, value) in self.df.iterrows():
            if chrm not in ly.chrm2angles:
                continue
            angle1 = ly.loc2angle(chrm, beg)
            anglewid = ly.loc2angle(chrm, end) - angle1
            angle1s.append(angle1)
            anglewids.append(anglewid)
            heights.append(self.track_height*(value-self.minval)/self.valrange*self.maxr)

        # ha='center', va='center', rotation=normalize_text_angle(text_angle/(np.pi*2)*360,tangent=True),
        
        # plotting background
        bg_start = []
        bg_width = []
        bg_height = []
        if background:
            for chrm in ly.chrm2angles:
                angle_beg = ly.loc2angle(chrm,0)
                angle_end = ly.loc2angle(chrm,ly.chrm2len[chrm])
                bg_start.append(angle_beg)
                bg_width.append(angle_end-angle_beg)
                bg_height.append(self.track_height*self.maxr)
        ly.ax.bar(bg_start, bg_height, bg_width, bottom=self.track_bottom, ec='none', fc='grey', alpha=0.1)
        
        # plot data
        ly.ax.bar(angle1s, heights, anglewids, bottom=self.track_bottom, ec='none', **kwargs1)
        if 'labelside_text_angle' in self.plot_kwargs:
            text_angle = self.plot_kwargs['labelside_text_angle'] / 180.0 * np.pi
        else:
            text_angle = ly.angle_beg
        ly.ax.text(text_angle, self.track_bottom + self.track_height/2.0, label, fontname=ly.fontname, va='center', ha='left', **kwargs2)

        if 'labelside_circle' in self.plot_kwargs:
            fc = self.plot_kwargs['fc'] if 'fc' in self.plot_kwargs else 'r'
            if 'labelside_circ_angle' in self.plot_kwargs:
                circ_angle = self.plot_kwargs['labelside_circ_angle']
            else:
                circ_angle = ly.angle_beg
            ly.ax.add_artist(mpatches.Circle(
                polar2cart(circ_angle, self.track_bottom),
                self.plot_kwargs['labelside_circle'],
                edgecolor=fc, color=fc, alpha=0.4, lw=0.1, transform=ly.ax.transData._b))

# usage:
# cl = CircosLayout("/Users/wandingzhou/reference/hg19/hg19.fa")
# cl = CircosLayout(faidx.RefGenome("/Users/wandingzhou/reference/hg19/hg19.fa"))
# cl = CircosLayout([('Alu',10000), ('L1', 50000), ..])
# cl.plot_chromosomes()
# cl.plot_interaction_list([(('Alu','L1'),20000), (('Alu','Alu'), 10000), ...])
class CircosLayout:

    def __init__(self,
                 genome,        # faidx.RefGenome
                 cytoband=None,
                 angle_beg=0, angle_end=360,
                 angle_chrm_space=0.02, inner_radius=98, outer_radius=101,
                 angle_inc=0.01, # when plotting an arc, this is the increment
                 bezier_anchor = 15.0, # higher the number the higher the peak
                 track_height = 10,
                 track_space = 1,
                 fontname = 'Arial Narrow',
    ):

        genome_path = genome
        if isinstance(genome, str):
            genome = faidx.RefGenome(genome)
        if isinstance(genome, list):
            chrms = genome
        else:
            chrms = getsortedchrms(genome)

        if cytoband:
            if re.search('mm10', genome_path) is not None:
                self.cytoband_table = pd.read_table('/Users/wandingzhou/references/mm10/cytoband.tsv',header=None, names=['_','chrm','band','band2','beg','end','bandtype'])[['chrm','beg','end','bandtype']]
            elif re.search('hg38', genome_path) is not None:
                self.cytoband_table = pd.read_table('/Users/wandingzhou/references/hg38/cytoBand.txt', header=None, names=['chrm','beg','end','band','bandtype'])[['chrm','beg','end','bandtype']]
        else:
            self.cytoband_table = None

        self.chrms = [_[0] for _ in chrms]
        self.chrms_plot = self.chrms[:] # plot all by default
        self.chrm2len = dict(chrms)
        self.chrm2angles = {}

        ## plot parameters
        self.angle_beg = angle_beg / 180.0 * np.pi
        self.angle_end = angle_end / 180.0 * np.pi
        self.angle_chrm_space = angle_chrm_space # space between chromosomes
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.angle_inc = angle_inc
        self.bezier_anchor = bezier_anchor
        self.ax = None
        self.fontname = fontname
        self.track_height = track_height
        self.track_space = track_space
        
        ## derived parameters
        self.connect_radius = inner_radius*0.95
        # self.ax_radius = self.outer_radius*1.05
        self.inner_track_outer_radius = self.inner_radius - self.track_space
        self.inner_track_inner_radius = self.inner_track_outer_radius - self.track_height
        self.outer_track_inner_radius = self.outer_radius + self.track_space
        self.outer_track_outer_radius = self.outer_track_inner_radius + self.track_height

        self.tracks = []

        return

    def show(self):

        plt.show()

    def save(self, fn, dpi=300):

        plt.savefig(fn, bbox_inches='tight', dpi=dpi)

    def loc2angle(self, chrm, pos):

        return self.chrm2angles[chrm][0] + self.angle_per_base*pos
    
    def plot_chromosomes(self, ax=None, figsize=(10,10), fc='blue',
                         chrms_plot=None, # which chromosome to plot
                         shade=False, fontsize=10,
                         ticks=False, tickspace=10000000, chrmtext_padding=5):

        if ax is not None:
            self.ax = ax

        if self.ax is None:
            fig = plt.figure(figsize=figsize)
            self.ax = fig.add_axes([0.1,0.1,0.9,0.9],projection='polar')
            self.ax.set_axis_off()

        if chrms_plot is not None:
            self.chrms_plot = chrms_plot

        self.totalbases = sum([self.chrm2len[_] for _ in self.chrms_plot])

        ## setup chrm2angles
        self.angle_per_base = (self.angle_end - self.angle_beg -
                               self.angle_chrm_space*len(self.chrms_plot)) / float(self.totalbases)

        tmp = self.angle_beg
        for chrm in self.chrms_plot:
            tmp1 = tmp + self.chrm2len[chrm]*self.angle_per_base
            self.chrm2angles[chrm] = (tmp, tmp1) # (start angle, end angle)
            tmp = tmp1
            tmp += self.angle_chrm_space

        self.chrm2color = dict(zip(self.chrms_plot, wzcolors.getncolors(len(self.chrms_plot),cm='Paired')))
        for chrm in self.chrms_plot:
            angles = self.chrm2angles[chrm]
            self.ax.add_patch(wzpolygon(angles[0], angles[1], self.inner_radius, self.outer_radius, ec='k', fc=self.chrm2color[chrm]))
            if shade:
                self.ax.bar(angles[0],self.inner_radius,angles[1]-angles[0],bottom=0.0, ec='none', fc='grey', alpha=0.1)
            # plot text
            text_angle = (angles[0]+angles[1])/2.0
            self.ax.text(text_angle, self.outer_radius+chrmtext_padding, chrm, ha='center', va='center', rotation=normalize_text_angle(text_angle/(np.pi*2)*360,tangent=True), fontsize=fontsize, fontname=self.fontname)

        # cytoband
        if self.cytoband_table is not None:
            for k, (chrm,beg,end,bandtype) in self.cytoband_table.iterrows():
                if chrm not in self.chrm2angles:
                    continue
                angle_beg = self.loc2angle(chrm,beg)
                angle_end = self.loc2angle(chrm,end)
                self.ax.bar(angle_beg, self.outer_radius-self.inner_radius,
                            width=angle_end-angle_beg, bottom=self.inner_radius,
                            ec='none', fc=cytoband_color[bandtype], alpha=cytoband_alpha[bandtype])

        if ticks:
            tickheight_major = 2
            tickheight_minor = 1
            ticklines = []
            for chrm in self.chrms_plot:
                chrmlen = self.chrm2len[chrm]
                chrmbeg = self.chrm2angles[chrm][0]
                for pos1M in xrange(chrmlen/tickspace+1):
                    tickangle = chrmbeg+pos1M*tickspace*self.angle_per_base
                    if pos1M % 5 == 0:
                        tick_height = tickheight_major
                        self.ax.text(tickangle, self.outer_radius+tick_height+1.5, str(pos1M*round_to_10(tickspace)), ha='center', va='center', rotation=normalize_text_angle(tickangle/(np.pi*2)*360), fontsize=8, fontname=self.fontname)
                    else:
                        tick_height = tickheight_minor
                    ticklines.append([(tickangle, self.outer_radius),
                                      (tickangle, self.outer_radius+tick_height)])

            # print ticklines
            self.ax.add_collection(mcollections.LineCollection(ticklines, colors='k', linewidth=0.5))

        # self.ax.set_xlim(-self.ax_radius,self.ax_radius)
        # self.ax.set_ylim(0,self.ax_radius)

    def _plot_arc(self, angle1, angle2, radius):
        angle_num = max(3,int(float(angle2 - angle1) / self.angle_inc))
        angle_num = min(500, angle_num)
        for angle in np.linspace(angle1, angle2, angle_num):
            self.verts.append((angle, radius))
            self.codes.append(mpath.Path.LINETO)


    def _plot_bezier(self, angle1, angle2, radius, anchor_angle=None, anchor_radius=None):

        while abs(angle2 - angle1) > np.pi:
            if angle1 < angle2:
                angle1 += np.pi*2
            else:
                angle2 += np.pi*2

        if anchor_radius is None:
            anchor_radius = radius*np.cos(abs(angle2 - angle1)/2)/self.bezier_anchor

        ### through trial and error, special treatment
        if anchor_radius < radius/20.:
            anchor_radius = -45 # the more negative the closer the cross the center connection to the center

        if anchor_angle is None:
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
            self, _angle1, _angle2, fc='grey',
            alpha=0.4,
    ):
        
        angle1 = min(_angle1, _angle2)
        angle2 = max(_angle1, _angle2)

        # connect
        self._init_path()
        self._plot_bezier(angle1, angle2, self.connect_radius)

        # arc
        self._plot_arc(angle1, angle2, self.connect_radius)
        self._close_path(angle1, self.connect_radius)
        
        self.ax.add_patch(mpatches.PathPatch(self.path(), alpha=alpha, fc=fc, ec='none'))

    def plot_ribbon(
            self, _source_angle1, _source_angle2, _target_angle1, _target_angle2, fc='grey'
    ):

        source_angle1 = min(_source_angle1, _source_angle2)
        source_angle2 = max(_source_angle1, _source_angle2)
        target_angle1 = min(_target_angle1, _target_angle2)
        target_angle2 = max(_target_angle1, _target_angle2)

        # source 2 to target 1
        self._init_path()

        source_angle_mid = (source_angle1+source_angle2)/2.
        target_angle_mid = (target_angle1+target_angle2)/2.
        while abs(target_angle_mid - source_angle_mid) > np.pi:
            if source_angle_mid < target_angle_mid:
                source_angle_mid += np.pi*2
            else:
                target_angle_mid += np.pi*2

        anchor_angle = (source_angle_mid+target_angle_mid)/2.
        anchor_radius = self.connect_radius*np.cos(abs(source_angle_mid - target_angle_mid)/2)/self.bezier_anchor
        
        self._plot_bezier(source_angle2, target_angle1, self.connect_radius,
                          anchor_angle=anchor_angle, anchor_radius=anchor_radius)
        self._plot_arc(target_angle1, target_angle2, self.connect_radius)
        self._plot_bezier(target_angle2, source_angle1, self.connect_radius,
                          anchor_angle=anchor_angle, anchor_radius=anchor_radius)
        self._plot_arc(source_angle1, source_angle2, self.connect_radius)
        self._close_path(source_angle2, self.connect_radius)

        self.ax.add_patch(mpatches.PathPatch(self.path(), alpha=0.4, fc=fc, ec='none'))

    def chrm2angle_span(self, chrm):

        beg, end = self.chrm2angles[chrm]
        return end-beg
        
    def plot_interaction_list(self, interactions):

        ## squeeze all counts into the angle span
        # calculate supp2angle
        supp2angle = {}
        chrm2totalcnt = {}
        for inter, supp in interactions:
            for chrm in set(inter):
                if chrm in self.chrm2angles:
                    if chrm in chrm2totalcnt:
                        chrm2totalcnt[chrm] += supp
                    else:
                        chrm2totalcnt[chrm] = supp

        # plot interaction
        chrm2cumulative_angle = {}
        for inter, supp in interactions:
            if len(inter) == 1: # no colocalization
                chrm = inter[0]
                if chrm not in self.chrm2angles:
                    continue
                if chrm not in chrm2cumulative_angle:
                    chrm2cumulative_angle[chrm] = self.chrm2angles[chrm][0]
                cum_angle_beg = chrm2cumulative_angle[chrm]
                cum_angle_end = cum_angle_beg+supp/float(chrm2totalcnt[chrm])*self.chrm2angle_span(chrm)
                self.plot_enclosure(cum_angle_beg, cum_angle_end, fc=self.chrm2color[chrm])
                chrm2cumulative_angle[chrm] = cum_angle_end
            elif len(inter) == 2:
                if inter[0] != inter[1]: # between chromosome
                    chrm1 = inter[0]
                    chrm2 = inter[1]
                    if chrm1 not in self.chrm2angles or chrm2 not in self.chrm2angles:
                        continue
                    if chrm1 not in chrm2cumulative_angle:
                        chrm2cumulative_angle[chrm1] = self.chrm2angles[chrm1][0]
                    if chrm2 not in chrm2cumulative_angle:
                        chrm2cumulative_angle[chrm2] = self.chrm2angles[chrm2][0]
                    cum_angle_beg1 = chrm2cumulative_angle[chrm1]
                    cum_angle_end1 = cum_angle_beg1+supp/float(chrm2totalcnt[chrm1])*self.chrm2angle_span(chrm1)
                    cum_angle_beg2 = chrm2cumulative_angle[chrm2]
                    cum_angle_end2 = cum_angle_beg2+supp/float(chrm2totalcnt[chrm2])*self.chrm2angle_span(chrm2)
                    self.plot_ribbon(cum_angle_beg1, cum_angle_end1, cum_angle_beg2, cum_angle_end2)
                    chrm2cumulative_angle[chrm1] = cum_angle_end1
                    chrm2cumulative_angle[chrm2] = cum_angle_end2
                else:           # within chromosome bivalent
                    chrm = inter[0]
                    if chrm not in self.chrm2angles:
                        continue
                    if chrm not in chrm2cumulative_angle:
                        chrm2cumulative_angle[chrm] = self.chrm2angles[chrm][0]
                    cum_angle_beg = chrm2cumulative_angle[chrm]
                    cum_angle_end = cum_angle_beg+supp/float(chrm2totalcnt[chrm])*self.chrm2angle_span(chrm)
                    self.plot_enclosure(cum_angle_beg, cum_angle_end, fc=self.chrm2color[chrm], alpha=0.7)
                    chrm2cumulative_angle[chrm] = cum_angle_end
            
                
    def plot_ribbon_chrm(
            self, ax,
            _source_chrm, _source_beg, _source_end,
            _target_chrm, _target_beg, _target_end, fc='grey'):

        self.plot_ribbon(ax, self.loc2angle(_source_chrm, _source_beg),
                         self.loc2angle(_source_chrm, _source_end),
                         self.loc2angle(_target_chrm, _target_beg),
                         self.loc2angle(_target_chrm, _target_end),)

    def add_bar_track(self, df, direction='inner', track_height=None, **plot_kwargs):

        """ track is [(chrm, beg, end, value)]
        norm_fac is the normalization factor for the height, sometimes we want
        different track to have different height

        @ params: fc, alpha, label_fontsize (9), label
        """

        t = CircosTrackBar(df, self, **plot_kwargs)
        self.tracks.append(t)
        t.direction = direction

        # set track_height
        if track_height is None:
            t.track_height = self.track_height
        else:
            t.track_height = track_height

        return

    def plot(self, track_height_proportion=False, track_background=True, inner_plot_frac=0.8):

        """ plot everything """

        # renormalize height of all the inner tracks
        inner_tracks = [t for t in self.tracks if t.direction == 'inner']
        total_track_heights = 0
        for t in inner_tracks:
            if track_height_proportion:
                t.track_height = t.df.iloc[:,3].max()
            total_track_heights += t.track_height
        inner_track_height_range = self.inner_radius*inner_plot_frac
        _inner_track_inner_radius = self.inner_radius
        for t in inner_tracks:
            t.track_height = inner_track_height_range/float(total_track_heights)*t.track_height
            _inner_track_inner_radius -= t.track_height
            t.track_bottom = _inner_track_inner_radius

        # set height and bottom of outer tracks
        outer_tracks = [t for t in self.tracks if t.direction == 'outer']
        _outer_track_inner_radius = self.outer_radius
        for t in outer_tracks:
            t.track_height = self.track_height
            _outer_track_inner_radius += self.track_height
            t.track_bottom = _outer_track_inner_radius

        for t in self.tracks:
            t.plot(background=track_background)
