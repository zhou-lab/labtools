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

def get_spectral_colors_rgb_dark(n, lightness=0.4, saturation=0.25):
    return [colorsys.hls_to_rgb(1.0/n*i,lightness,saturation) for i in xrange(n)]

def get_distinct_colors_hex(num_colors):
    return ['#%02x%02x%02x' % (int(r*256), int(g*256), int(b*256))
            for r,g,b in get_distinct_colors_rgb(num_colors)]

def get_grey_scale_rgb(num_colors, greyscale_range=(0.1,0.9)):
    darkest, lightest = greyscale_range
    if num_colors == 1:
        return [(lightest,lightest,lightest,1)]
    # print [(c,c,c,1.) for c in np.linspace(darkest,lightest,num_colors)]
    return [(c,c,c,1.) for c in np.linspace(darkest,lightest,num_colors)]

def map_distinct_colors_hex(data, other2grey=False, greyscale=False, greyscale_range=(0.1,0.9)):

    levels = set(data)
    if greyscale:
        colors = get_grey_scale_rgb(len(levels), greyscale_range=greyscale_range)
    else:
        colors = get_distinct_colors_hex(len(levels))
    level2color = dict(zip(levels, colors))
    if other2grey and "other" in level2color:
        level2color["other"] = "#E6E6E6"
    return ([level2color[datum] for datum in data], level2color)

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
