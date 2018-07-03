"""
Generate intervals for meta-gene plot.
Usage:
python ~/wzlib/pyutils/wzmetagene.py ~/projects/xiaotian/Canyons.conservation.bed | les
"""

import argparse
import numpy as np

def main(args):


    if args.flankByBase > 0:
        args.flank = args.flankByBase * args.numflank

    for line in args.table:
        fields = line.strip('\n').split('\t')
        chrm = fields[0]
        beg = int(fields[1])
        end = int(fields[2])

        if args.strand is None:
            strand = '+'
        else:
            strand = fields[args.strand-1]

        if strand == '+':
            # upstream intervals
            if not args.noflank:
                sentinels = np.linspace(beg - args.flank, beg, args.numflank+1)
                for i in xrange(len(sentinels)-1):
                    window_beg = int(sentinels[i])
                    window_end = int(sentinels[i+1])
                    if window_end > window_beg:
                        if args.middle:
                            window_mid = int((window_beg + window_end)/2)
                            window_beg = window_mid-1
                            window_end = window_mid

                        index = i - args.numflank

                        # -1 for upstream
                        if window_beg > 0:
                            print('%s\t%d\t%d\t%d\t%d\t-1\t%s' % (chrm, window_beg, window_end, index, index*args.flank/args.numflank, '\t'.join(fields)))

            # internal intervals
            sentinels = np.linspace(beg, end, args.numinternal+1)
            for i in xrange(len(sentinels)-1):
                window_beg = int(sentinels[i])
                window_end = int(sentinels[i+1])
                if window_end > window_beg:
                    if args.middle:
                        window_mid = int((window_beg + window_end)/2)
                        window_beg = window_mid-1
                        window_end = window_mid

                    if args.fold and args.strand is None:
                        index = min(args.numinternal - i, i + 1)
                    else:
                        index = i+1
                      
                    # 0 for internal
                    if window_beg > 0:
                        print('%s\t%d\t%d\t%d\t%d%%\t0\t%s' % (chrm, window_beg, window_end, index, float(index-1)/args.numinternal*100,'\t'.join(fields)))

            # downstream intervals
            if not args.noflank:
                sentinels = np.linspace(end, end + args.flank, args.numflank+1)
                for i in xrange(len(sentinels)-1):
                    window_beg = int(sentinels[i])
                    window_end = int(sentinels[i+1])
                    if window_end > window_beg:
                        if args.middle:
                            window_mid = int((window_beg + window_end)/2)
                            window_beg = window_mid-1
                            window_end = window_mid

                        if args.fold and args.strand is None: 
                            index = -i-1
                        else:
                            index = 1+i+args.numinternal

                        # +1 for upstream
                        if window_beg > 0:
                            print('%s\t%d\t%d\t%d\t%d\t1\t%s' % (chrm, window_beg, window_end, index, index*args.flank/args.numflank, '\t'.join(fields)))

        else:
            # upstream
            if not args.noflank:
                sentinels = np.linspace(end, end + args.flank, args.numflank+1)
                for i in xrange(len(sentinels)-2,0,-1):
                    window_end = int(sentinels[i+1])
                    window_beg = int(sentinels[i])
                    if window_end > window_beg:
                        if args.middle:
                            window_mid = int((window_beg + window_end)/2)
                            window_beg = window_mid-1
                            window_end = window_mid

                        index = -i

                        # +1 for upstream
                        if window_beg > 0:
                            print('%s\t%d\t%d\t%d\t%d\t1\t%s' % (chrm, window_beg, window_end, index, index*args.flank/args.numflank, '\t'.join(fields)))

            # internal
            sentinels = np.linspace(beg, end, args.numinternal+1)
            for i in xrange(len(sentinels)-2,0,-1):
                window_end = int(sentinels[i+1])
                window_beg = int(sentinels[i])
                if window_end > window_beg:
                    if args.middle:
                        window_mid = int((window_beg + window_end)/2)
                        window_beg = window_mid-1
                        window_end = window_mid

                    if args.fold:
                        index = min(args.numinternal - i, i + 1)
                    else:
                        index = args.numinternal - i

                    # 0 for internal
                    if window_beg > 0:
                        print('%s\t%d\t%d\t%d\t%d\t0\t%s' % (chrm, window_beg, window_end, index, float(index-1)/args.numinternal*100, '\t'.join(fields)))

            # downstream
            if not args.noflank:
                sentinels = np.linspace(beg - args.flank, beg, args.numflank+1)
                for i in xrange(len(sentinels)-2,0,-1):
                    window_end = int(sentinels[i+1])
                    window_beg = int(sentinels[i])
                    if window_end > window_beg:
                        if args.middle:
                            window_mid = int((window_beg + window_end)/2)
                            window_beg = window_mid-1
                            window_end = window_mid

                        index = i + args.numinternal

                        # -1 for upstream
                        if window_beg > 0:
                            print('%s\t%d\t%d\t%d\t%d\t-1\t%s' % (chrm, window_beg, window_end, index, index*args.flank/args.numflank, '\t'.join(fields)))


    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate meta gene')
    parser.add_argument('table', help="Input bed file", type = argparse.FileType('r'), default='-')
    # report middle point for each sampled interval
    parser.add_argument('--middle', action = 'store_true', help = 'use middle point of each interval as the sentinel')
    # controls flanking length
    parser.add_argument('-f', '--flank', default=10000, type=int, help = 'length of flanking to plot, default 10kb')
    parser.add_argument('-F', '--flankByBase', type = int, default=-1, help = 'plot each X bases for flanking sequences, by default false (-1), this overrides -f')
    parser.add_argument('--noflank', action = 'store_true', help = 'suppress flanking region output')
    parser.add_argument('-m', '--numflank', type = int, default=30, help = 'number of points to sample in the flanking region')
    # controls internal sampling
    parser.add_argument('-n', '--numinternal', type = int, default=30, help = 'number of points to sample in the genic/internal region, --middle ignores this')
    # others
    parser.add_argument('--fold', action = 'store_true', help = 'use the same index for intervals from two sides the target, usually used when strand is irrelevant')
    parser.add_argument('-s', '--strand', type=int, default=None, help = 'the field which contains strand information, if None then ignore strand')

    parser.set_defaults(func=main)
    args = parser.parse_args()
    try:
      args.func(args)
    except IOError:
      exit
