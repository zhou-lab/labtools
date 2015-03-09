#!/usr/bin/env python
import argparse
import faidx
""" convert bed to bed graph """

def main(args):

    ref = faidx.RefGenome(args.ref)
    prev_chrm = None
    prev_beg = None
    prev_end = None
    for line in args.bed:
        fields = line.strip().split('\t')
        chrm = fields[0]
        beg = int(fields[1])
        end = int(fields[2])
        val = float(fields[3])
        if chrm != prev_chrm or beg != prev_beg or end != prev_end:
            if not prev_end is None:
                print '\t'.join(map(str, [prev_chrm, prev_beg-1, prev_end, tval]))
            tval = 0.0
            if chrm != prev_chrm:
                if not prev_end is None:
                    print '\t'.join(map(str, [prev_chrm, prev_end, ref.chrm2len(prev_chrm), 0]))
                print '\t'.join(map(str, [chrm, 0, beg-1, 0]))
            elif prev_end != beg-1:
                print '\t'.join(map(str, [chrm, prev_end, beg-1, 0]))
                
        if args.op == 'sum':
            tval += val

        prev_chrm = chrm
        prev_beg = beg
        prev_end = end

    if prev_chrm is not None:
        print '\t'.join(map(str, [prev_chrm, prev_beg-1, prev_end, tval]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.set_defaults(func=main)
    parser.add_argument("-bed", help="data table", type=argparse.FileType('r'), default='-')
    parser.add_argument("-ref", help="reference file", required=True)
    parser.add_argument("-op", help="operations {sum, }", default="sum")

    args = parser.parse_args()

    try:
        args.func(args)
    except IOError:
        pass
