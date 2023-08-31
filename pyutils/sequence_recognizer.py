#!/usr/bin/env python3
import argparse
from fuzzysearch import find_near_matches

def check_one_read_standalone(read, seq, max_len):
    return find_near_matches(seq, read, max_l_dist=max_len)

# python3 -c 'import SHAREseq_tools; SHAREseq_tools.true_false_table("SEPA", "seqs", "nums", "1.csv")'
# python3 -c 'import SHAREseq_tools; SHAREseq_tools.table_analysis("1.csv",[1], 0, 100000000000)'
# python3 -c 'import SHAREseq_tools; SHAREseq_tools.check_one_read_standalone("ATGCATGCATGC", "ATGC", 0)'

parser = argparse.ArgumentParser(description="This is a command line tool to analyze your SHARE-seq data or other data that includes barcodes.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-c", "--contains", action="store_true", help="Very simple tool, returns T/F for whether a read contains a sequence.")
group.add_argument("-t", "--table", action="store_true", help="Prints a table with information of the matches.")
group.add_argument("-a", "--annotate", action="store_true", help="Annotates the finding on the sequence.")
parser.add_argument("-r", "--read", type=str, help="The enter read.")
parser.add_argument("-s", "--sequences", type=str, nargs='+', help="The sequences that should be contained.")
parser.add_argument("-m", "--max_distances", type=int, default=2, nargs='+', help="The maximum edit distance.")
args = parser.parse_args()

if len(args.sequences) != len(args.max_distances):
    print("The sequences and the maximum distances do not match!")
    exit()

for (seq, max_len) in zip(args.sequences, args.max_distances):
    if args.contains:
        print(str(len(check_one_read_standalone(args.read, seq, max_len)) > 0), end="\t")
    elif args.table:
        print("Sequence: " + seq)
        print(check_one_read_standalone(args.read, seq, max_len))
    elif args.annotate:
        print("Sequence: " + seq)
        for match in check_one_read_standalone(args.read, seq, max_len):
            CRED = '\033[91m'
            CEND = '\033[0m'
            print("Match:", end=" ")
            print(args.read[:match.start] + CRED + args.read[match.start:match.end] + CEND + args.read[match.end:])


# Example runs:
# python3 sequence_recognizer.py -c -r ACTTCCGTTCAGTTACGTATTGCTAGATGTGTATAAGAGACAGAGGATAAATCATGCTAAGGCGAGGATGAAGCCGATATCGCCGATACGGTTGTATAGGATTGCTTGAATGGCTGCTGTGTTGGCATCTGCTCGGGCGTATCATCAACTGATGAGCAAAAAGAAGGATGTAATTCCGTT -s CTGTCTCTTATACAC GTGTATAAGAGACAG -m 2 2
# False   True    
# python3 sequence_recognizer.py -t -r ACTTCCGTTCAGTTACGTATTGCTAGATGTGTATAAGAGACAGAGGATAAATCATGCTAAGGCGAGGATGAAGCCGATATCGCCGATACGGTTGTATAGGATTGCTTGAATGGCTGCTGTGTTGGCATCTGCTCGGGCGTATCATCAACTGATGAGCAAAAAGAAGGATGTAATTCCGTT -s CTGTCTCTTATACAC GTGTATAAGAGACAG -m 2 2
# Sequence: 
# []
# Sequence: 
# [Match(start=28, end=43, dist=0, matched='GTGTATAAGAGACAG')]
# python3 sequence_recognizer.py -a -r ACTTCCGTTCAGTTACGTATTGCTAGATGTGTATAAGAGACAGAGGATAAATCATGCTAAGGCGAGGATGAAGCCGATATCGCCGATACGGTTGTATAGGATTGCTTGAATGGCTGCTGTGTTGGCATCTGCTCGGGCGTATCATCAACTGATGAGCAAAAAGAAGGATGTAATTCCGTT -s CTGTCTCTTATACAC GTGTATAAGAGACAG -m 2 2
# Sequence: CTGTCTCTTATACAC
# Sequence: GTGTATAAGAGACAG
# ACTTCCGTTCAGTTACGTATTGCTAGATGTGTATAAGAGACAGAGGATAAATCATGCTAAGGCGAGGATGAAGCCGATATCGCCGATACGGTTGTATAGGATTGCTTGAATGGCTGCTGTGTTGGCATCTGCTCGGGCGTATCATCAACTGATGAGCAAAAAGAAGGATGTAATTCCGTT (highlighted)




