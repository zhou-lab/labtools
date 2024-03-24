#!/usr/bin/env python3
from fuzzysearch import find_near_matches
import fastq as fq
import os
import argparse
from sequence_operations_common import reverse_complement

def check_one_read_standalone(read, seq, max_l_dist):
    return find_near_matches(seq, read, max_l_dist=max_l_dist)

def sort_key_start(obj):
    return obj.start

def read_bam(file, check_sq=False):
    reads = []
    bamfile = pysam.AlignmentFile(file, "rb", check_sq=check_sq)
    header = bamfile.header
    for read in bamfile:
        reads.append(read)
    bamfile.close()
    return reads, header

def write_bam(reads, file, header):
    output_bamfile = pysam.AlignmentFile(file, "wb", header=header)
    for read in reads:
        output_bamfile.write(read)
    output_bamfile.close()

parser = argparse.ArgumentParser(description="This is a program to demultiplex the nanopore output data. So far, it only supports dual indices.")
parser.add_argument("-f", "--file", type=str, help="The bam file that includes all reads.")
parser.add_argument("-i", "--id", type=str, help="The id of this dataset.")
parser.add_argument("-p", "--possible_barcodes_file", type=str, help="The sequences that should be contained.")
parser.add_argument("-m", "--max_distance", type=int, default=3, help="The maximum edit distance.")
parser.add_argument("-e", "--edge_distance", type=int, default=50, help="The edge distance the barcodes should be.")
args = parser.parse_args()

possible_barcodes = []
demultiplexed_bam_objs = [[], []]
with open(args.possible_barcodes_file, 'r') as file:
    possible_barcodes = file.read().splitlines()
    for _ in possible_barcodes: 
        demultiplexed_bam_objs.append([])

bam_reads, header = read_bam(args.file)
for bam_read in bam_reads:
    seq = bam_read.query_sequence
    if seq is None:
        continue
    all_matches = []
    for index, barcode in enumerate(possible_barcodes):
        matches = []
        founds = check_one_read_standalone(seq, barcode, args.max_distance)
        founds.extend(check_one_read_standalone(seq, reverse_complement(barcode), args.max_distance))
        for found in founds:
            if not (found.start > args.edge_distance and found.end < len(seq) - args.edge_distance):
                matches.append(found)
        if matches:
            all_matches.append(index)
    if len(all_matches) == 0:
        demultiplexed_bam_objs[len(demultiplexed_bam_objs)-2].append(bam_read) # no match: last one
    elif len(all_matches) == 1:
        demultiplexed_bam_objs[all_matches[0]].append(bam_read)
    else:
        demultiplexed_bam_objs[len(demultiplexed_bam_objs)-1].append(bam_read) # multiple matches: one before last one

for index, demultiplexed_per_barcode in enumerate(demultiplexed_bam_objs):
    if index < len(possible_barcodes):
        file = args.file.replace(args.id, "").replace("/pass/", "")
        output_directory = args.id + "_barcode" + str(index+1) + "/pass"
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        output_bam_path = os.path.join(output_directory, file)
    elif index == len(possible_barcodes):
        output_directory = args.id + "_no_match/pass"
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        output_bam_path = os.path.join(output_directory, file)
    elif index == len(possible_barcodes)+1:
        output_directory = args.id + "_multiple_matches/pass"
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        output_bam_path = os.path.join(output_directory, file)
    write_bam(demultiplexed_per_barcode, output_bam_path, header)
