#!/usr/bin/env python

import sys
import argparse

class Indices:

    def __init__(self):
        self.spans = []

    def extend(self, start, end):
        self.spans.append((start, end))

    def extract(self, lst):
        result = []
        for start, end in self.spans:
            if not end:
                end = len(lst)
            result.extend([lst[_] for _ in xrange(start, end)])

        return result

def parse_indices(indstr):
    rgs = indstr.split(',')
    indices = Indices()
    for rg in rgs:
        if rg.find('-') >= 0:
            pair = rg.split('-')
            if not pair[0]:
                pair[0] = 0
            if not pair[1]:
                pair[1] = None
            indices.extend(int(pair[0])-1 if pair[0] else 0,
                           int(pair[1]) if pair[1] else None)
        else:
            indices.extend(int(rg)-1, int(rg))

    return indices

def main_reorder(args):

    indices = parse_indices(args.c)

    for i, line in enumerate(args.table):
        fields = line.strip().split(args.delim)
        print '\t'.join(indices.extract(fields))

    return

def main_transpose(args):

    data = []
    for line in args.table:
        fields = line.strip('\n').split(args.delim)
        data.append(fields)

    for fields in zip(*data):
        print args.delim.join(fields)

    return


def _build_map(args, index):

    args_c = getattr(args, "c%d" % index)
    args_p = getattr(args, "p%d" % index)
    args_fc = getattr(args, "fc%d" % index)
    args_fp = getattr(args, "fp%d" % index)

    c = parse_indices(args_c) if args_c else None
    p = parse_indices(args_p) if args_p else None
    lmap = {}
    for line in getattr(args, "t%d" % index):
        f = line.strip().split("\t")
        key = '\t'.join(c.extract(f)) if c else args_fc.format(f=f)
        if p:
            val = '\t'.join(p.extract(f))
        elif args_fp:
            val = args_fp.format(f=f)
        else:
            val = line.strip()
            
        lmap[key] = val

    return lmap


def main_compare3(args):

    map1 = _build_map(args, 1)
    set1 = set(map1.keys())

    map2 = _build_map(args, 2)
    set2 = set(map2.keys())

    map3 = _build_map(args, 3)
    set3 = set(map3.keys())

    if not args.n:
        sys.stderr.write("123: %d\n" % len(set1 & set2 & set3))
        sys.stderr.write("12not3: %d\n" % len((set1 & set2) - set3))
        sys.stderr.write("23not1: %d\n" % len((set2 & set3) - set1))
        sys.stderr.write("13not2: %d\n" % len((set1 & set3) - set2))
        sys.stderr.write("1not23: %d\n" % len(set1 - set2 - set3))
        sys.stderr.write("2not13: %d\n" % len(set2 - set1 - set3))
        sys.stderr.write("3not12: %d\n" % len(set3 - set1 - set2))

    if args.p == '123':
        for e in set1 & set2 & set3:
            print e, map1[e], map2[e], map3[e]
    if args.p == '12not3':
        for e in (set1 & set2) - set3:
            print e, map1[e], map2[e]
    if args.p == '13not2':
        for e in (set1 & set3) - set2:
            print e, map1[e], map3[e]
    if args.p == '23not1':
        for e in (set2 & set3) - set1:
            print e, map2[e], map3[e]
    if args.p == '1not23':
        for e in set1 - set2 - set3:
            print map1[e]
    if args.p == '2not13':
        for e in set2 - set1 - set3:
            print map2[e]
    if args.p == '3not12':
        for e in set3 - set1 - set2:
            print map3[e]

    return
    
def main_compare(args):

    # if 3-way comparison
    if args.t3:
        return main_compare3(args)

    map1 = _build_map(args, 1)
    set1 = set(map1.keys())

    map2 = _build_map(args, 2)
    set2 = set(map2.keys())

    if not args.n:
        sys.stderr.write("12: %d\n" % len(set1 & set2))
        sys.stderr.write("1not2: %d\n" % len(set1 - set2))
        sys.stderr.write("2not1: %d\n" % len(set2 - set1))

    if args.p == '1and2':
        for e in set1 & set2:
            print map1[e]+'\t'+map2[e]
    if args.p == '1or2':
        for e in set1 | set2:
            print e
    if args.p == '1not2':
        for e in set1 - set2:
            print map1[e]
    if args.p == '2not1':
        for e in set2 - set1:
            print map2[e]
    if args.p == '1':
        for e in set1:
            print map1[e]
    if args.p == '2':
        for e in set2:
            print map2[e]

def main_tabulate(args):

    """ tabulate a list into table """

    rc_map = {}
    colnames = []
    if args.cnt:
        for line in args.list:
            fields = line.strip().split(args.delim)
            row = fields[args.r-1]
            col = fields[args.c-1]
            if row not in rc_map:
                rc_map[row] = {}
            if col in rc_map[row]:
                rc_map[row][col] += 1
            else:
                rc_map[row][col] = 1
            if col not in colnames:
                colnames.append(col)
    else:
        for line in args.list:
            fields = line.strip().split(args.delim)
            row = fields[args.r-1]
            col = fields[args.c-1]
            data = fields[args.d-1]
            if row not in rc_map:
                rc_map[row] = {}
            if col not in rc_map[row]:
                rc_map[row][col] = data
            if col not in colnames:
                colnames.append(col)

    rownames = rc_map.keys()
    if args.s:
        colnames.sort()
        rownames.sort()

    print 'header\t%s' % '\t'.join(colnames)
    for row in rownames:
        line = row
        for col in colnames:
            if col in rc_map[row]:
                line += '\t%s' % str(rc_map[row][col])
            else:
                line += '\t0' if args.cnt else '\tNA'
        print line

def main_untabulate(args):

    """ untabulate a table into a list """
    
    rownames = args.table.readline().strip('\n').split('\t')
    for line in args.table:
        fields = line.strip().split('\t')
        for i, d in enumerate(fields):
            if i>0:
                print '%s\t%s\t%s' % (fields[0], rownames[i], d)
    
    return


def main_match(args):

    """ match two table by common column """

    c1 = parse_indices(args.c1) if args.c1 else None
    c2 = parse_indices(args.c2) if args.c2 else None
    p1 = parse_indices(args.p1) if args.p1 else None
    p2 = parse_indices(args.p2) if args.p2 else None

    key2prints = {}
    for line in args.t1:
        f = line.strip().split(args.delim)
        key = '\t'.join(c1.extract(f)) if c1 else args.fc1.format(f=f)
        if p1:
            val = '\t'.join(p1.extract(f))
        elif args.fp1:
            val = args.fp1.format(f=f)
        else:
            val = None
        key2prints[key] = val

    for line in args.t2:
        f = line.strip().split(args.delim)
        key = '\t'.join(c2.extract(f)) if c2 else args.fc2.format(f=f)

        if (key in key2prints):

            if p2:
                val = '\t'.join(p2.extract(f))
            elif args.fp2:
                val = args.fp2.format(f=f)
            else:
                val = None

            prnstr = key
            if key2prints[key]:
                prnstr += '\t'+key2prints[key]
            if val:
                prnstr += '\t'+val

            print prnstr

        elif args.pu:

            if p2:
                val = '\t'.join(p2.extract(f))
            elif args.fp2:
                val = args.fp2.format(f=f)
            else:
                val = None

            prnstr = key+'\tUNMATCHED'
            if val:
                prnstr += '\t'+val
            print prnstr

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Manipulate tables')
    subparsers = parser.add_subparsers()

    parser_reorder = subparsers.add_parser("reorder", help="reorder columns in a table")
    parser_reorder.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    parser_reorder.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_reorder.add_argument('-c', default=None, help="columns to be printed in the output, 1-based. E.g., -c 1,3-4 [None]")
    parser_reorder.set_defaults(func=main_reorder)

    parser_transpose = subparsers.add_parser("transpose", help="transpose table")
    parser_transpose.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    parser_transpose.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_transpose.set_defaults(func=main_transpose)

    parser_compare = subparsers.add_parser("compare", help="compare two columns of two table")
    parser_compare.add_argument("-t1", type=argparse.FileType('r'), default=None, help="data table 1")
    parser_compare.add_argument("-t2", type=argparse.FileType('r'), default=None, help="data table 2")
    parser_compare.add_argument("-t3", type=argparse.FileType('r'), default=None, help="data table 3 (optional)")
    parser_compare.add_argument("-c1", default=None, help="column(s) to be compared in table 1, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-c2", default=None, help="column(s) to be compared in table 2, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-c3", default=None, help="column(s) to be compared in table 3, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-fc1", default=None, help="format key in table 1, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-fc2", default=None, help="format key in table 2, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-fc3", default=None, help="format key in table 3, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-p", choices=[None, '1not2', '2not1', '1and2', '1or2', '1', '2', '3', '12not3', '13not2', '23not1', '1not23', '2not13', '3not12', '123'], default=None, help="optional print")
    parser_compare.add_argument('-n', action="store_true", help="suppress the statistics output")
    parser_compare.add_argument("-p1", default=None, help="column(s) to be compared in table 1, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-p2", default=None, help="column(s) to be compared in table 2, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-p3", default=None, help="column(s) to be compared in table 3, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-fp1", default=None, help="format output in table 1, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-fp2", default=None, help="format output in table 2, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-fp3", default=None, help="format output in table 3, e.g., {f[0]}:{f[1]}")
    parser_compare.set_defaults(func=main_compare)


    parser_tabulate = subparsers.add_parser("tabulate", help="tabulate list")
    parser_tabulate.add_argument("-r", type=int, required=True, 
                                 help="column for row names (1-based)")
    parser_tabulate.add_argument("-c", type=int, required=True,
                                 help="column for column names (1-based)")
    parser_tabulate.add_argument("-d", type=int, help="column for data (1-based, optional)")
    parser_tabulate.add_argument("list", type=argparse.FileType('r'),
                                 default='-')
    parser_tabulate.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_tabulate.add_argument('-s', action="store_true", help="have the rows and columns sorted by their names")
    parser_tabulate.add_argument('--cnt', action="store_true", help="count occurrence as data")
    parser_tabulate.set_defaults(func=main_tabulate)

    parser_untabulate = subparsers.add_parser("untabulate", help="untabulate a table into list")
    parser_untabulate.add_argument("table", type=argparse.FileType('r'), default='-')
    parser_untabulate.set_defaults(func=main_untabulate)

    # this function is similar to compare, but not exactly
    parser_match = subparsers.add_parser("match", help="help match two table by common column")
    parser_match.add_argument("-t1", type=argparse.FileType('r'), required=True, help='table 1 (used as key-value map)')
    parser_match.add_argument("-t2", type=argparse.FileType('r'), required=True, help='table 2 (one to be translated)')
    parser_match.add_argument('-c1', default=None, help="column(s) to match, e.g., 1,3,5 (1-based)")
    parser_match.add_argument('-c2', default=None, help="column(s) to match, e.g., 1,3,5 (1-based)")
    parser_match.add_argument('-p1', default=None, help="column(s) to print in table 1, e.g., 4,5,6")
    parser_match.add_argument('-p2', default=None, help="column(s) to print in table 2, e.g., 3,4")
    parser_match.add_argument("-fc1", default=None, help="format match in table 1, e.g., {f[0]},{f[2]}:{f[4]} (fields are 0-based)")
    parser_match.add_argument("-fc2", default=None, help="format match in table 2, e.g., {f[0]},{f[2]}:{f[4]} (fields are 0-based)")
    parser_match.add_argument('-fp1', default=None, help="format print in table 1, e.g., {f[0]},{f[2]}:{f[4]}")
    parser_match.add_argument('-fp2', default=None, help="format print in table 2, e.g., {f[0]},{f[2]}")
    parser_match.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_match.add_argument('-pu', action='store_true', help='print unmatched entry in table 2')
    parser_match.set_defaults(func=main_match)
    
    args = parser.parse_args()
    args.func(args)
