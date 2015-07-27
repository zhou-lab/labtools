import colorsys
import itertools
from fractions import Fraction

def zenos_dichotomy():
    """
    http://en.wikipedia.org/wiki/1/2_%2B_1/4_%2B_1/8_%2B_1/16_%2B_%C2%B7_%C2%B7_%C2%B7
    """
    for k in itertools.count():
        yield Fraction(1,2**k)

def getfracs():
    """
    [Fraction(0, 1), Fraction(1, 2), Fraction(1, 4), Fraction(3, 4), Fraction(1, 8), Fraction(3, 8), Fraction(5, 8), Fraction(7, 8), Fraction(1, 16), Fraction(3, 16), ...]
    [0.0, 0.5, 0.25, 0.75, 0.125, 0.375, 0.625, 0.875, 0.0625, 0.1875, ...]
    """
    yield 0
    for k in zenos_dichotomy():
        i = k.denominator # [1,2,4,8,16,...]
        for j in range(1,i,2):
            yield Fraction(j,i)

bias = lambda x: (math.sqrt(x/3)/Fraction(2,3)+Fraction(1,3))/Fraction(6,5) # can be used for the v in hsv to map linear values 0..1 to something that looks equidistant

def genhsv(h):
    for s in [Fraction(6,10)]: # optionally use range
        for v in [Fraction(8,10),Fraction(5,10)]: # could use range too
            yield (h, s, v) # use bias for v here if you use range

genrgb = lambda x: colorsys.hsv_to_rgb(*x)

flatten = itertools.chain.from_iterable

gethsvs = lambda: flatten(itertools.imap(genhsv,getfracs()))

getrgbs = lambda: itertools.imap(genrgb, gethsvs())

import numpy as np

def get_spectral_colors_rgb(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (30 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append(colorsys.hls_to_rgb(hue, lightness, saturation))
        
    return colors

""" adjacent objects have distinct colors """

def get_distinct_colors_rgb(num_colors):
    return list(itertools.islice(itertools.imap(
        lambda x: colorsys.hsv_to_rgb(*x), gethsvs()), num_colors))

def get_distinct_colors_html(num_colors):
    return list(itertools.islice(itertools.imap(
        lambda x: "rgb({},{},{})".format(*[int(xx*255) for xx in x]),
        getrgbs()), num_colors))

def _round_hex(x):
    rx = int(x*256)
    if rx == 256: rx -= 1
    return rx

def hex2rgb(x):

    if x.startswith('#'):
        x = x[1:]

    return tuple(map(lambda _: ord(_)/256.0, x.decode('hex')))

def get_spectral_colors_rgb_dark(n, lightness=0.4, saturation=0.25):
    return [colorsys.hls_to_rgb(1.0/n*i,lightness,saturation) for i in xrange(n)]

def map2l_spectral_colors(p2c, lightness=0.4, saturation=0.25):
    """ map two levels of spectral colors """
    
    nump = len(p2c)
    pcolors_hl = [(1.0/nump*i, lightness) for i in xrange(nump)]
    p2colors = {}
    c2colors = {}
    for i, (p, cl) in enumerate(p2c.iteritems()):
        h,l = pcolors_hl[i]
        p2colors[p] = colorsys.hls_to_rgb(h,l,saturation)
        for c, s in zip(cl, np.linspace(saturation, 0.8, len(cl))):
            c2colors[c] = colorsys.hls_to_rgb(h,l,s)

    return p2colors, c2colors

def map2sub_saturation(p2c, p2color):

    c2color = {}
    for p, cl in p2c.iteritems():
        cl = [_ for _ in cl if not isnan(_)]
        pcolor = p2color[p]
        if isinstance(pcolor, str):
            pcolor = hex2rgb(pcolor)
        elif max(pcolor) > 1.0:
            pcolor = tuple([_ / 256.0 for _ in pcolor])

        h, l, saturation = colorsys.rgb_to_hls(*pcolor)
        for c, s in zip(cl, np.linspace(saturation, 0.8, len(cl))):
            c2color[c] = colorsys.hls_to_rgb(h,l,s)

    return c2color

def map2sub_alpha(p2c, p2color):

    c2color = {}
    for p, cl in p2c.iteritems():
        cl = [_ for _ in cl if not isnan(_)]
        pcolor = p2color[p]
        if isinstance(pcolor, str):
            r,g,b = hex2rgb(pcolor)
        elif max(pcolor) > 1.0:
            r,g,b = tuple([_ / 256.0 for _ in pcolor])

        for c, a in zip(cl, np.linspace(1.0, 0.3, len(cl))):
            c2color[c] = (r,g,b,a)

    return c2color


def map2sub_light(p2c, p2color):

    c2color = {}
    for p, cl in p2c.iteritems():
        cl = [_ for _ in cl if not isnan(_)]
        pcolor = p2color[p]
        if isinstance(pcolor, str):
            pcolor = hex2rgb(pcolor)
        elif max(pcolor) > 1.0:
            pcolor = tuple([_ / 256.0 for _ in pcolor])

        h, darkest, s = colorsys.rgb_to_hls(*pcolor)
        for c, l in zip(cl, np.linspace(darkest, 0.6, len(cl))):
            c2color[c] = colorsys.hls_to_rgb(h,l,s)

    return c2color

        
def get_distinct_colors_hex(num_colors):
    return ['#%02x%02x%02x' % (int(r*256), int(g*256), int(b*256))
            for r,g,b in get_distinct_colors_rgb(num_colors)]

def get_grey_scale_rgb(num_colors, greyscale_range=None):
    if greyscale_range is None:
        greyscale_range = (0.1,0.9)
    darkest, lightest = greyscale_range
    if num_colors == 1:
        return [(lightest,lightest,lightest,1)]
    # print [(c,c,c,1.) for c in np.linspace(darkest,lightest,num_colors)]
    return [(c,c,c,1.) for c in np.linspace(darkest,lightest,num_colors)]

def kwargs_or_none(kwargs, op):
    if op in kwargs:
        return kwargs[op]
    else:
        return None

def kwargs_or_false(kwargs, op):
    if op in kwargs:
        return kwargs[op]
    else:
        return False

def map_level2color(data, **kwargs):

    """
    Note: nan is not included in level2color

    kwargs: other2grey=False, greyscale=False, greyscale_range=(0.1,0.9)
    """
    levels = set([_ for _ in data if not isnan(_)])

    if levels == set([True, False]):
        levels = [True, False]
        colors = ['r', '#A9F5BC']
    elif 'greyscale' in kwargs and kwargs['greyscale']:
        colors = get_grey_scale_rgb(len(levels), greyscale_range=kwargs_or_none(kwargs, 'greyscale_range'))
    else:
        colors = get_distinct_colors_hex(len(levels))
        
    level2color = dict(zip(levels, colors))
    if kwargs_or_false(kwargs, 'other2grey') and "other" in level2color:
        level2color["other"] = "#E6E6E6"

    if kwargs_or_false(kwargs, 'other2grey') and 'NA' in level2color:
        level2color['NA'] = "#E6E6E6"

    return level2color

def isnan(x):

    try:
        return np.isnan(x)
    except TypeError:
        return False

def map_distinct_colors_hex(data, **kwargs):

    """ kwargs: other2grey=False, greyscale=False, greyscale_range=(0.1,0.9)
    """

    level2color = map_level2color(data, **kwargs)

    return (['#E6E6E6' if isnan(datum) else level2color[datum] for datum in data], level2color)


color_curation = {

    'IndianRed'



}

if __name__ == "__main__":
    # print(list(itertools.islice(gethtmlcolors(), 1,100)))
    # print get_distinct_colors_html(100)
    # print get_distinct_colors_rgb(100)
    # print _get_colors(8)
    # print get_distinct_colors_hex(10)
    # print map_distinct_colors_hex([1,2,3,1,2,3])
    # import matplotlib.pyplot as plt
    # ncolors = 30
    # plt.pie([1]*ncolors, colors=get_distinct_colors_rgb(ncolors))
    # plt.pie([1]*ncolors, colors=get_spectral_colors_rgb(ncolors))
    # plt.pie([1]*ncolors, colors=get_distinct_colors_hex(ncolors))
    # plt.pie([1]*ncolors, colors=get_distinct_colors_rgb(ncolors))
    # plt.colorbar()
    # plt.show()

    pass
