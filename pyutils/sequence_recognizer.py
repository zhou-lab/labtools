#!/usr/bin/env python3
import argparse
from fuzzysearch import find_near_matches
from sequence_operations_common import reverse_complement

def check_one_read_standalone(read, seq, max_len):
    return find_near_matches(seq, read, max_l_dist=max_len)

def sort_key_dist(obj):
    return obj.dist

parser = argparse.ArgumentParser(description="This is a command line tool to analyze your SHARE-seq data or other data that includes barcodes.")
group = parser.add_mutually_exclusive_group()
group.add_argument("-U", "--ultimate_checking", action="store_true", help="Returns true if it is a good read for a split seq based sequencing")
group.add_argument("-c", "--contains", action="store_true", help="Very simple tool, returns T/F for whether a read contains a sequence.")
group.add_argument("-t", "--table", action="store_true", help="Prints a table with information of the matches. The columns are: start, end, edit distance, matched sequence in the read.")
group.add_argument("-a", "--annotate", action="store_true", help="Annotates the finding on the sequence.")
parser.add_argument("-f", "--flip", action="store_true", help="If you include this, all the reverse complementary sequences of the sequences you enter will also be computed.")
parser.add_argument("-r", "--read", type=str, help="The enter read.")
parser.add_argument("-s", "--sequences", type=str, nargs='+', help="The sequences that should be contained.")
parser.add_argument("-m", "--max_distances", type=int, default=2, nargs='+', help="The maximum edit distance.")
args = parser.parse_args()

if len(args.sequences) != len(args.max_distances):
    print("The sequences and the maximum distances do not match!")
    exit()

sequences = [element for element in args.sequences]
max_distances = [element for element in args.max_distances]

if args.ultimate_checking:
    contain_all = True
    contain_all_flip = True
    for (seq, max_len) in zip(sequences, max_distances):
        if not len(check_one_read_standalone(args.read, seq, max_len)) > 0:
            contain_all = False
        if not len(check_one_read_standalone(args.read, reverse_complement(seq), max_len)) > 0:
            contain_all_flip = False
    if contain_all:
        forward = []
        forward_check = True
        for (seq, max_len) in zip(sequences, max_distances):
            forward.append(sorted(check_one_read_standalone(args.read, seq, max_len), key=sort_key_dist)[0])
        for i in range(len(forward)):
            if i < len(forward)-1:
                if forward[i].end > forward[i+1].start:
                    forward_check = False
                    break
        if forward_check:
            print("F", end="\t")
            for i in range(len(forward)):
                # if i < len(forward)-1:
                #     print(args.read[forward[i].end:forward[i+1].start], end="\t")
                print(str(forward[i].start) + "\t" + str(forward[i].end) + "\t" + str(forward[i].dist) + "\t" + forward[i].matched, end="\t")
        else:
            print("False", end="\t")
    elif contain_all_flip:
        backward = []
        backward_check = True
        for (seq, max_len) in zip(sequences, max_distances):
            backward.append(sorted(check_one_read_standalone(args.read, reverse_complement(seq), max_len), key=sort_key_dist)[0])
        for i in range(len(backward)):
            if i < len(backward)-1:
                if backward[i].start < backward[i+1].end:
                    backward_check = False
                    break
        if backward_check:
            print("B", end="\t")
            for i in range(len(backward)):
                # if i < len(backward)-1:
                #     print(args.read[backward[i+1].end:backward[i].start], end="\t")
                print(str(backward[i].start) + "\t" + str(backward[i].end) + "\t" + str(backward[i].dist) + "\t" + backward[i].matched, end="\t")
        else:
            print("False", end="\t")
    else:
        print("False", end="\t")
    print()
    exit()

if args.flip:
    for (seq, max_len) in zip(args.sequences, args.max_distances):
        sequences.append(reverse_complement(seq))
        max_distances.append(max_len)


for (seq, max_len) in zip(sequences, max_distances):
    if args.contains:
        print(str(len(check_one_read_standalone(args.read, seq, max_len)) > 0), end="\t")
    elif args.table:
        #print("Sequence: " + seq)
        if len(check_one_read_standalone(args.read, seq, max_len)) > 0:
            matched_dic = check_one_read_standalone(args.read, seq, max_len)[0]
            print("Seq:" + seq + "\t"  + str(matched_dic.start) + "\t" + str(matched_dic.end) + "\t" + str(matched_dic.dist) + "\t" + str(matched_dic.matched), end="\t")
        #else:
        #    print("-1\t-1\t-1\tX")
    elif args.annotate:
        print("Sequence: " + seq)
        for match in check_one_read_standalone(args.read, seq, max_len):
            CRED = '\033[91m'
            CEND = '\033[0m'
            print("Match:", end=" ")
            print(args.read[:match.start] + CRED + args.read[match.start:match.end] + CEND + args.read[match.end:])

print() ## This is to give a new line.

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




