import pysam
import argparse
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def main(args):

    samfile=pysam.Samfile(args.bam)
    for read in samfile.fetch(args.reg):
	if read.qname != args.reg:
            continue
	for op, len in read.cigar:
            if op == 
		


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r')
    parser.add_argument('-reg')
    parser.add_argument('-ref')
    parser.add_argument('-bam')
    parser.set_defaults(func=main)
    args = parser.parse_args()
    args.func(args)

