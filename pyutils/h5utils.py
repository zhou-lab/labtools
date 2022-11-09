import sys
import os
import gzip
import h5py
import numpy as np
import argparse

def pipeprint(s):
    try:
        print(s)
    except IOError:
        sys.exit(1)
        
def opengz(fn,m='r'):

    if fn == "-":
        fh = sys.stdin
    else:
        if not os.path.isfile(fn):
            print("File path {} does not exist. Exiting...".format(fn))
            sys.exit()
        if fn.endswith('.gz'):
            import gzip
            fh = gzip.open(fn, m)
        else:
            fh = open(fn,m)

    return fh

def main_importMU(args):
    
    MU = []
    with opengz(args.input) as fh:
        for i, line in enumerate(fh):
            MU.append(int(line.strip()))
            # fields = line.strip().split("\t")
            # MU.append((int(fields[0]), int(fields[1])))

    # note that this is different from compound datatype
    # dt = np.dtype([("M", np.uint32), ("U", np.uint32)])
    pairs = args.out.split(".h5")
    if len(pairs) != 2:
        sys.exit(1)
    fpath_out = pairs[0]+".h5"
    fpath_dat = pairs[1]
    with h5py.File(fpath_out, "w") as f:
        f.create_dataset(
            fpath_dat, data=MU, dtype=np.uint32,
            compression="gzip", compression_opts=9)

def main_importBinary(args):

    BI = []
    with opengz(args.input) as fh:
        for i, line in enumerate(fh):
            fields = line.strip().split("\t")
            BI.append(int(fields[0]))

    pairs = args.out.split(".h5")
    if len(pairs) != 2:
        sys.exit(1)
    fpath_out = pairs[0]+".h5"
    fpath_dat = pairs[1]
    with h5py.File(fpath_out, "w") as f:
        f.create_dataset(
            fpath_dat, data=BI, dtype=np.uint8,
            compression="gzip", compression_opts=9)

def main_importGenome(args):
    chrmset = []
    chrm = []
    pos = []
    with opengz(args.input) as fh:
        for i, line in enumerate(fh):
            fields = line.decode().strip().split("\t")
            if fields[0] not in chrmset:
                chrmset.append(fields[0])
            chrm.append(chrmset.index(fields[0]))
            pos.append(int(fields[1])+1)

    with h5py.File(args.out,"r+") as f:
        f.create_dataset("rowData/chrmset", data=chrmset, compression="gzip")
        f.create_dataset("rowData/chrm", data=chrm, dtype=np.uint8, compression="gzip")
        f.create_dataset("rowData/pos", data=pos, dtype=np.uint64, compression="gzip")

def main_print(args):

    with h5py.File(args.input, "r") as f:
        chrmset = [x.decode() for x in f["rowData/chrmset"][...]]
        chrm = f["rowData/chrm"][...]
        pos = f["rowData/pos"][...]
        dset = f[args.samples]
        if isinstance(dset, h5py._hl.group.Group):
            dset = dset[list(dset.keys())[0]]

        dset = dset[...]
        for i,d in enumerate(dset):
            pipeprint("%s\t%d\t%d\t%d\t%d" % (chrmset[chrm[i]], pos[i]-1, pos[i]+1, d[0], d[1]))


def main_cbind(args):

    ## /bin/ls tmp/step2b/*.h5 | awk '{split($1,a,"/"); split(a[3],aa,"."); print $0"/MU/"aa[1];}' | python2022 ~/zhoulab/labtools/pyutils/h5utils.py cbind - myout.h5/MU
    ## this takes longer for big datasets. takes about 3 hours for a 2k datasets, and has to be done serially
    
    fnames = []
    fpaths = []
    with opengz(args.input) as fh:
        fnames = [line.strip() for line in fh]
        
        for fname in fnames:
            pairs = fname.split(".h5/")
            if len(pairs) != 2:
                sys.exit(1)
            fpaths.append((pairs[0]+".h5", pairs[1]))

    dset0 = h5py.File(fpaths[0][0], "r")[fpaths[0][1]]
    pairs = args.out.split(".h5")
    if len(pairs) != 2:
        sys.exit(1)
    fpath_out = pairs[0]+".h5"
    fpath_dat = pairs[1]
    fout = h5py.File(fpath_out,"w")
    nrow = dset0.shape[0]
    ncol = len(fnames)
    dset = fout.create_dataset(fpath_dat, shape=(nrow, ncol), dtype = dset0.dtype, chunks = True, compression = "gzip")
    print("Creating dataset %s/%s: %d x %d" % (fpath_out, fpath_dat, nrow, ncol))
    for i, (fpath, dat) in enumerate(fpaths):
        print("Adding %d: %s/%s" % (i, fpath, dat))
        with h5py.File(fpath,"r") as f:
            dset[:,i] = f[dat]

    fout.close()
        
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(prog='h5utils')
    subs = parser.add_subparsers()
    
    sub = subs.add_parser('importGenome', help='Import genome coordinates')
    sub.add_argument('input', help='input file')
    sub.add_argument('out', help='out.h5')
    sub.set_defaults(func=main_importGenome)

    sub = subs.add_parser('importMU', help='Import MU')
    sub.add_argument('input', help='input file')
    sub.add_argument('out', help='out.h5/path')
    sub.set_defaults(func=main_importMU)

    sub = subs.add_parser('importBinary', help='Import Binary')
    sub.add_argument('input', help='input file')
    sub.add_argument('out', help='out.h5')
    sub.set_defaults(func=main_importBinary)

    sub = subs.add_parser('print', help='Print sample data')
    sub.add_argument('input', help='in.h5 file')
    sub.add_argument("samples", help="sample names (if None, print the first sample)")
    sub.set_defaults(func=main_print)

    sub = subs.add_parser("cbind", help="column-bind, while chunking, should use this in parallel")
    sub.add_argument('input', help='a list of in.h5/path files')
    sub.add_argument("out", help="out.h5/path file")
    sub.set_defaults(func=main_cbind)
    
    args = parser.parse_args()
    try:
        func = args.func
    except AttributeError:
        parser.error("too few arguments")

    args.func(args)
