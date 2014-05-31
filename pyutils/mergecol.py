#!/RIS/HPC_apps/AMD/python/python-2.7.2/bin/python

import argparse
import re
import sys

def main(args):

    fhs = []
    colinds = []
    for fncol in args.tables:
        pair = re.split(":|,", fncol)
        fn = pair[0]
        colinds.append([int(_)-1 for _ in pair[1:]])
        fhs.append(open(fn))

    for line in fhs[0]:
        pair = line.strip().split(args.d)
        s = '\t'.join([pair[j] for j in colinds[0]])
        for i in xrange(1, len(fhs)):
            line = fhs[i].readline()
            pair = line.strip().split(args.d)
            if colinds[i]:
                s += '\t'
                s += '\t'.join([pair[j] for j in colinds[i]])
        try:
            sys.stdout.write("%s\n" % s)
            sys.stdout.flush()
        except IOError:
            break

    for fh in fhs:
        fh.close()
        

if __name__ == '__main__':

    psr = argparse.ArgumentParser(description='merge columns from multiple files', epilog="""
Example:
mergecol.py file1:3 file2:5
 """)
    psr.add_argument('tables', nargs="+", help="table file name with the columns to extract. Format: fn:col1,col2,col3")
    psr.add_argument('-d', default="\t", help="table delimiter [\\t]")
    psr.add_argument('--skipheader', action='store_true', help='skip header')
    psr.set_defaults(func=main)


    args = psr.parse_args()
    args.func(args)
