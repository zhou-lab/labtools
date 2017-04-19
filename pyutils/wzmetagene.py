"""
Generate intervals for meta-gene plot.
Usage:
python ~/wzlib/pyutils/wzmetagene.py ~/projects/xiaotian/Canyons.conservation.bed | les
"""

import argparse
import numpy as np

def main(args):

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
            sentinels = np.linspace(beg - args.flank, beg, args.numflank+1)
            for i in xrange(len(sentinels)-1):
                window_beg = int(sentinels[i])
                window_end = int(sentinels[i+1])
                if window_end > window_beg:
                    # -1 for upstream
                    print('%s\t%d\t%d\t%d\t-1\t%s' % (chrm, window_beg, window_end, i-args.numflank, '\t'.join(fields)))

            # internal intervals
            sentinels = np.linspace(beg, end, args.numinternal+1)
            for i in xrange(len(sentinels)-1):
                window_beg = int(sentinels[i])
                window_end = int(sentinels[i+1])
                if window_end > window_beg:
                    # 0 for internal
                    print('%s\t%d\t%d\t%d\t0\t%s' % (chrm, window_beg, window_end, i+1, '\t'.join(fields)))

            # downstream intervals
            sentinels = np.linspace(end, end + args.flank, args.numflank+1)
            for i in xrange(len(sentinels)-1):
                window_beg = int(sentinels[i])
                window_end = int(sentinels[i+1])
                if window_end > window_beg:
                    # +1 for upstream
                    print('%s\t%d\t%d\t%d\t1\t%s' % (chrm, window_beg, window_end, 1+i+args.numinternal, '\t'.join(fields)))
        else:
            # upstream
            sentinels = np.linspace(end, end + args.flank, args.numflank+1)
            for i in xrange(len(sentinels),0,-1):
                window_beg = int(sentinels[i+1])
                window_end = int(sentinels[i])
                if window_end > window_beg:
                    # +1 for upstream
                    print('%s\t%d\t%d\t%d\t1\t%s' % (chrm, window_beg, window_end, -i, '\t'.join(fields)))

            # internal
            sentinels = np.linspace(beg, end, args.numinternal+1)
            for i in xrange(len(sentinels),0,-1):
                window_beg = int(sentinels[i+1])
                window_end = int(sentinels[i])
                if window_end > window_beg:
                    # 0 for internal
                    print('%s\t%d\t%d\t%d\t0\t%s' % (chrm, window_beg, window_end, args.numinternal-i, '\t'.join(fields)))

            # downstream
            sentinels = np.linspace(beg - args.flank, beg, args.numflank+1)
            for i in xrange(len(sentinels),0,-1):
                window_beg = int(sentinels[i+1])
                window_end = int(sentinels[i])
                if window_end > window_beg:
                    # -1 for upstream
                    print('%s\t%d\t%d\t%d\t-1\t%s' % (chrm, window_beg, window_end, i+args.numinternal, '\t'.join(fields)))


    return

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Generate meta gene')
    parser.add_argument('table', help="Input bed file", type = argparse.FileType('r'), default='-')
    parser.add_argument('-f', '--flank', default=10000, help = 'length of flanking to plot, default 10kb')
    parser.add_argument('-m', '--numflank', type = int, default=30, help = 'number of points to sample in the flanking region')
    parser.add_argument('-n', '--numinternal', type = int, default=30, help = 'number of points to sample in the genic/internal region')
    parser.add_argument('-s', '--strand', default=None, help = 'the strand information, if None then ignore strand')

    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)
