#!/usr/bin/env python

import sys, re
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
    indices = Indices()
    if not indstr: return indices
    rgs = indstr.split(',')
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

    if args.n:
        for i in xrange(args.nskip):
            args.table.readline()
        header = args.table.readline().strip()
        names = args.n.split(',')
        indices = [i for i, name in enumerate(header.split(args.delim)) if name in names]
        for line in args.table:
            fields = line.strip().split(args.delim)
            print '\t'.join([fields[i] for i in indices])
    else:
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
    args_np = getattr(args, "np%d" % index)

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
        elif args_np:
            val = ''
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
            prncol = []
            if not args.rk: prncol.append(e)
            if map1[e]: prncol.append(map1[e])
            if map2[e]: prncol.append(map2[e])
            if map3[e]: prncol.append(map3[e])
            print '\t'.join(prncol)

    if args.p == '12not3':
        for e in (set1 & set2) - set3:
            prncol = []
            if not args.rk: prncol.append(e)
            if map1[e]: prncol.append(map1[e])
            if map2[e]: prncol.append(map2[e])
            print '\t'.join(prncol)

    if args.p == '13not2':
        for e in (set1 & set3) - set2:
            prncol = []
            if not args.rk: prncol.append(e)
            if map1[e]: prncol.append(map1[e])
            if map3[e]: prncol.append(map3[e])
            print '\t'.join(prncol)

    if args.p == '23not1':
        for e in (set2 & set3) - set1:
            prncol = []
            if not args.rk: prncol.append(e)
            if map2[e]: prncol.append(map2[e])
            if map3[e]: prncol.append(map3[e])
            print '\t'.join(prncol)

    if args.p == '1not23':
        for e in set1 - set2 - set3:
            prncol = []
            if not args.rk: prncol.append(e)
            if map1[e]: prncol.append(map1[e])
            print '\t'.join(prncol)

    if args.p == '2not13':
        for e in set2 - set1 - set3:
            prncol = []
            if not args.rk: prncol.append(e)
            if map2[e]: prncol.append(map2[e])
            print '\t'.join(prncol)

    if args.p == '3not12':
        for e in set3 - set1 - set2:
            prncol = []
            if not args.rk: prncol.append(e)
            if map3[e]: prncol.append(map3[e])
            print '\t'.join(prncol)

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
            prncol = []
            if not args.rk: prncol.append(e)
            if map1[e]: prncol.append(map1[e])
            if map2[e]: prncol.append(map2[e])
            print '\t'.join(prncol)

    if args.p == '1or2':
        for e in set1 | set2:
            prncol = []
            if not args.rk: prncol.append(e)
            if e in map1 and map1[e]: prncol.append(map1[e])
            if e in map2 and map2[e]: prncol.append(map2[e])
            print '\t'.join(prncol)

    if args.p == '1not2':
        for e in set1 - set2:
            prncol = []
            if not args.rk: prncol.append(e)
            if map1[e]: prncol.append(map1[e])
            print '\t'.join(prncol)

    if args.p == '2not1':
        for e in set2 - set1:
            prncol = []
            if not args.rk: prncol.append(e)
            if map2[e]: prncol.append(map2[e])
            print '\t'.join(prncol)

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
        if key in key2prints:
            key2prints[key].append(val)
        else:
            key2prints[key] = [val]

    keyprinted = set()
    for line in args.t2:
        f = line.strip().split(args.delim)
        key = '\t'.join(c2.extract(f)) if c2 else args.fc2.format(f=f)

        if p2:
            val2 = '\t'.join(p2.extract(f))
        elif args.fp2:
            val2 = args.fp2.format(f=f)
        else:
            val2 = None

        if (key in key2prints):

            if not args.sm:
                for val1 in key2prints[key]:

                    prncols = []
                    if not args.rk:
                        prncols.append(key)
                    if val1:
                        prncols.append(val1)
                    if val2:
                        prncols.append(val2)

                    print '\t'.join(prncols)

            keyprinted.add(key)

        elif args.um2:

            prncols = []
            if not args.rk:
                prncols.append(key)
            if val:
                prncols.append(val2)
                
            print '\t'.join(prncols)

    if args.um1:
        for key in key2prints:
            if key not in keyprinted:
                for val in key2prints[key]:
                    prncols = []
                    if not args.rk:
                        prncols.append(key)
                    if val:
                        prncols.append(val)
                    print '\t'.join(prncols)


def main_colindex(args):

    p = re.compile(args.r)
    for line in args.table:
        fields = line.strip().split(args.delim)
        print [(i+1, _) for i, _ in enumerate(fields) if p.search(_)]

def main_classify(args):

    k2v = {}
    for line in args.t:
        fields = line.strip().split(args.delim)
        k = fields[args.k-1]
        v = fields[args.v-1]
        if k in k2v:
            k2v[k].append(v)
        else:
            k2v[k] = [v]

    for k, vs in k2v.iteritems():
        prncol = []
        if not args.rk:
            prncol.append(k)
        if args.nv:
            prncol.append(str(len(vs)))
        prncol.append('\t'.join(vs))
        print '\t'.join(prncol)

def main_dedupmax(args):

    k2v = {}
    for line in args.t:
        fields = line.strip().split(args.delim)
        k = fields[args.k-1]
        v = (float(fields[args.v-1]), fields)
        if k in k2v:
            if k2v[k][0] < v[0]:
                k2v[k] = v
        else:
            k2v[k] = v

    for k, v in k2v.iteritems():
        print '\t'.join(v[1])

def main_dupcompress(args):

    ind = parse_indices(args.k)
    k2v = {}
    for line in args.t:
        fields = line.strip().split(args.delim)
        k = tuple(ind.extract(fields))
        if args.v <= len(fields):
            v = fields[args.v-1]
            if k in k2v:
                if v not in k2v[k]: k2v[k].append(v)
            else:
                k2v[k] = [v]

    for k, v in k2v.iteritems():
        print '%s\t%d\t%s' % ('\t'.join(k), len(v), '\t'.join(v))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Manipulate tables')
    subparsers = parser.add_subparsers()

    parser_reorder = subparsers.add_parser("reorder", help="reorder columns in a table")
    parser_reorder.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    parser_reorder.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_reorder.add_argument('-c', default=None, help="columns to be printed in the output, 1-based. E.g., -c 1,3-4 [None]")
    parser_reorder.add_argument('-n', default=None, help='header names, e.g., col1,col2 ...')
    parser_reorder.add_argument('-nskip', default=0, type=int, help='number of lines to skip before header')
    parser_reorder.set_defaults(func=main_reorder)

    parser_colindex = subparsers.add_parser("colindex", help="find the index of a particular column")
    parser_colindex.add_argument('table', help="data table", type = argparse.FileType('r'), default='-')
    parser_colindex.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_colindex.add_argument('-r', default=None, help='regular expression')
    parser_colindex.set_defaults(func=main_colindex)

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
    parser_compare.add_argument("-p", choices=[None, '1not2', '2not1', '1and2', '1or2', '12not3', '13not2', '23not1', '1not23', '2not13', '3not12', '123'], default=None, help="optional print")
    parser_compare.add_argument('-n', action="store_true", help="suppress the statistics output")
    parser_compare.add_argument("-p1", default=None, help="column(s) to be compared in table 1, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-p2", default=None, help="column(s) to be compared in table 2, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-p3", default=None, help="column(s) to be compared in table 3, e.g., '1,2-5' (1-based)")
    parser_compare.add_argument("-fp1", default=None, help="format output in table 1, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-fp2", default=None, help="format output in table 2, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument("-fp3", default=None, help="format output in table 3, e.g., {f[0]}:{f[1]}")
    parser_compare.add_argument('-rk', action='store_true', help='suppress key output')
    parser_compare.add_argument('-np1', action='store_true', help='no print of table 1')
    parser_compare.add_argument('-np2', action='store_true', help='no print of table 2')
    parser_compare.add_argument('-np3', action='store_true', help='no print of table 3')
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
    parser_match.add_argument('-um1', action='store_true', help='print unmatched entry in table 1')
    parser_match.add_argument('-um2', action='store_true', help='print unmatched entry in table 2')
    parser_match.add_argument('-sm', action='store_true', help='suppress match print')
    parser_match.add_argument('-rk', action='store_true', help='repress key output')
    parser_match.set_defaults(func=main_match)

    parser_classify = subparsers.add_parser("classify", help="classify a column by another column")
    parser_classify.add_argument('-t',type=argparse.FileType('r'), default='-')
    parser_classify.add_argument('-k', type=int, required=True, help="column for classification key (1-based)")
    parser_classify.add_argument('-v', type=int, required=True, help="column for classification value (1-based)")
    parser_classify.add_argument('-nv', action='store_true', help="output number of classification values")
    parser_classify.add_argument('-rk', action='store_true', help="no output of key")
    parser_classify.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_classify.set_defaults(func=main_classify)

    parser_dupmax = subparsers.add_parser("dedupmax", help="remove dup in one column by maxing a column")
    parser_dupmax.add_argument('-t',type=argparse.FileType('r'), default='-')
    parser_dupmax.add_argument('-k', type=int, required=True, help="column to dedup (1-based)")
    parser_dupmax.add_argument('-v', type=int, required=True, help="column to maximize (1-based)")
    parser_dupmax.add_argument('-rk', action='store_true', help="no output of key")
    parser_dupmax.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_dupmax.set_defaults(func=main_dedupmax)

    parser_dupcompress = subparsers.add_parser('dupcompress', help='remove dup in one column and list all value in another')
    parser_dupcompress.add_argument('-t',type=argparse.FileType('r'), default='-')
    parser_dupcompress.add_argument('-k', required=True, help="column to dedup (1-based)")
    parser_dupcompress.add_argument('-v', type=int, required=True, help="column to list (1-based)")
    parser_dupcompress.add_argument('--delim', default="\t", help="table delimiter [\\t]")
    parser_dupcompress.set_defaults(func=main_dupcompress)

    args = parser.parse_args()
    args.func(args)
