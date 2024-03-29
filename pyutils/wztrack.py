#!/usr/bin/env python
import argparse
import sys
# TCGA_HM450/hg38/output.bb
def main(args):

    if args.shortLabel is None:
        args.shortLabel = args.name
        
    if args.longLabel is None:
        args.longLabel = args.name
    
    if args.type == 'bigBed9':

        if args.indent: sys.stdout.write('    ')
        print('track %s' % args.name)
        if args.indent: sys.stdout.write('    ')
        print('type bigBed 9')
        if args.indent: sys.stdout.write('    ')
        print('shortLabel %s' % args.shortLabel)
        if args.indent: sys.stdout.write('    ')
        print('longLabel %s' % args.longLabel)
        if args.indent: sys.stdout.write('    ')
        print('bigDataUrl %s' % args.url)
        if args.indent: sys.stdout.write('    ')
        print('itemRgb on')
        if args.indent: sys.stdout.write('    ')
        print('parent %s on' % args.parent)
        print('')

    elif args.type == 'bigWigMeth':

        if args.indent: sys.stdout.write('    ')
        print('track %s' % args.name)
        if args.indent: sys.stdout.write('    ')
        print('type bigWig')
        if args.indent: sys.stdout.write('    ')
        print('shortLabel %s' % args.shortLabel)
        if args.indent: sys.stdout.write('    ')
        print('longLabel %s' % args.longLabel)
        if args.indent: sys.stdout.write('    ')
        print('bigDataUrl %s' % args.url)
        if args.indent: sys.stdout.write('    ')
        print('parent %s on' % args.parent)
        if args.indent: sys.stdout.write('    ')
        print('color 0,102,255')
        if args.indent: sys.stdout.write('    ')
        print('maxHeightPixels 128:25:10')
        if args.indent: sys.stdout.write('    ')
        print('viewLimits 0.0:1.0')
        print('')

    elif args.type == 'bigWigMethParent':

        print('track %s' % args.name)
        print('type bigWig')
        print('shortLabel %s' % args.shortLabel)
        print('longLabel %s' % args.longLabel)
        print('compositeTrack off')
        print('viewLimits 0.0:1.0')
        print('maxHeightPixels 100:15:5')
        print('visibility hide')
        print('allButtonPair on')
        print('dragAndDrop on')
        print('')

    elif args.type == 'bigBed9Parent':

        print('track %s' % args.name)
        print('type bigBed 9')
        print('shortLabel %s' % args.shortLabel)
        print('longLabel %s' % args.longLabel)
        print('superTrack on show')
        # print('autoScale on')
        print('maxHeightPixels 100:15:5')
        print('visibility hide')
        print('')

    elif args.type == 'bigBed9ParentComp':

        print('track %s' % args.name)
        print('type bigBed 9')
        print('shortLabel %s' % args.shortLabel)
        print('longLabel %s' % args.longLabel)
        print('compositeTrack off')
        # print('autoScale on')
        print('maxHeightPixels 100:15:5')
        print('visibility hide')
        print('')


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""create formatted track
Ex:
> wztrack.py --type bigBed9 --name xxx --url https://zwdzwd.s3.amazonaws.com/trackHubs/Mouse_MM285/mm10/204875570008_R06C01.bb --parent TCGA_ACC
produces
==
track xxx
type bigBed 9
shortLabel xxx
longLabel xxx
bigDataUrl https://zwdzwd.s3.amazonaws.com/trackHubs/Mouse_MM285/mm10/204875570008_R06C01.bb
itemRgb on
parent TCGA_ACC on
==""", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--name', type=str)
    parser.add_argument('--url')
    parser.add_argument('--parent')
    parser.add_argument('--shortLabel', default=None, type=str)
    parser.add_argument('--longLabel', default=None, type=str)
    parser.add_argument('--indent', action='store_true')
    parser.add_argument('--type', help="can be: bigBed9, bigWigMeth, bigWigMethParent, bigBed9Parent, bigBed9ParentComp", default='bigBed9')
    parser.set_defaults(func=main)

    args = parser.parse_args()
    try:
        args.func(args)
    except IOError as e:
        sys.exit()
