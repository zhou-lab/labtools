#!/usr/bin/env python3
from sequence_operations_common import check_one_read_standalone, reverse_complement, read_fastq, write_fastq
import fastq as fq
import argparse

def sort_key_dist(obj):
    return obj.dist

def parse_barcodes(direction, barcodes):
    barcodes_read = []
    for barcode in barcodes:
        if isinstance(barcode, float):
            barcodes_read.append("")
        else:
            barcodes_read.append(barcode)
    return barcodes_read

def sequence_extraction(read, sequences, max_distances):
    contain_all = True
    contain_all_flip = True
    extracted = []
    for (seq, max_len) in zip(sequences, max_distances):
        if not len(check_one_read_standalone(read, seq, max_len)) > 0:
            contain_all = False
        if not len(check_one_read_standalone(read, reverse_complement(seq), max_len)) > 0:
            contain_all_flip = False
    if contain_all:
        forward = []
        forward_check = True
        for (seq, max_len) in zip(sequences, max_distances):
            forward.append(sorted(check_one_read_standalone(read, seq, max_len), key=sort_key_dist)[0])
        for i in range(len(forward)):
            if i < len(forward)-1:
                if forward[i].end > forward[i+1].start:
                    forward_check = False
                    break
        if not forward_check:
            extracted = [False]
        else:
            forward.insert(0, "F")
            extracted = forward
    elif contain_all_flip:
        backward = []
        backward_check = True
        for (seq, max_len) in zip(sequences, max_distances):
            backward.append(sorted(check_one_read_standalone(read, reverse_complement(seq), max_len), key=sort_key_dist)[0])
        for i in range(len(backward)):
            if i < len(backward)-1:
                if backward[i].start < backward[i+1].end:
                    backward_check = False
                    break
        if not backward_check:
            extracted = [False]
        else:
            backward.insert(0, "B")
            extracted = backward
    else:
        extracted = [False]
    return extracted

def extract_barcodes(read, extracted_sequences, lengths):
    direction = extracted_sequences[0]
    linkers = extracted_sequences[1:]
    num_linkers = len(linkers)
    barcodes = []
    if lengths[1] != 0:
        if direction == "B":
            if linkers[num_linkers-1].start <= lengths[0]:
                return [False]
            else: 
                barcodes.append(read[linkers[num_linkers-1].start-lengths[1]:linkers[num_linkers-1].start])
    if lengths[0] != 0:
        if direction == "F":
            if linkers[0].start <= lengths[1]:
                return [False]
            else:
                barcodes.append(read[linkers[0].start-lengths[0]:linkers[0].start])
    for i in range(num_linkers):
        if direction == "F":
            if i < num_linkers-1:
                barcodes.append(read[linkers[i].end:linkers[i+1].start])
        elif direction == "B":
            if i < num_linkers-1:
                barcodes.append(read[linkers[num_linkers-1-i].end:linkers[num_linkers-2-i].start])
    if lengths[1] != 0:
        if direction == "F":
            if len(read)-linkers[num_linkers-1].end < lengths[0]:
                return [False]
            else: 
                barcodes.append(read[linkers[num_linkers-1].end:linkers[num_linkers-1].end+lengths[1]])
    if lengths[0] != 0:
        if direction == "B":
            if len(read)-linkers[0].end <= lengths[1]:
                return [False]
            else:
                barcodes.append(read[linkers[0].end:linkers[0].end+lengths[0]])
    extract_barcodes = [direction]
    if direction == "F":
        for barcode in barcodes:
            extract_barcodes.append(barcode)
    elif direction == "B":
        reverse = barcodes[::-1]
        for barcode in reverse:
            extract_barcodes.append(reverse_complement(barcode))
    return extract_barcodes

def barcodes_validation(extracted_barcodes, possible_barcodes):

    for i in range(len(extracted_barcodes)):
        barcode_found = False
        for j in range(len(possible_barcodes[i])):
            if extracted_barcodes[i] == possible_barcodes[i][j]:
                barcode_found = True
                break
        if not barcode_found:
            extracted_barcodes[i] = "ERROR"

    result_string = '_'.join(map(str, extracted_barcodes))
    if result_string.find("ERROR") == -1:
        return result_string
    return "ERROR"

parser = argparse.ArgumentParser(description="This is a command line tool designed to extract barcodes for SPLIT seq based single cell sequencing.")
parser.add_argument("-f", "--file", type=str, help="The fastq file that contains raw data (basecalled).")
parser.add_argument("-s", "--sequences", type=str, nargs='+', help="The sequences that should be contained.")
parser.add_argument("-m", "--max_distances", type=int, default=2, nargs='+', help="The maximum edit distance.")
parser.add_argument("-l", "--lengths", type=int, nargs='+', help="Lengths of the barcodes not bewteen the linkers.")
parser.add_argument("-t", "--trim_barcodes", type=int, nargs='+', help="If any barcode needs to be trimmed, enter here. If not, enter 0.")
parser.add_argument("-p", "--possible_barcodes_files", type=str, nargs='+', help="The text files that contain barcodes. How many files depends on how many rounds of barcodes.")
args = parser.parse_args()

possible_barcodes = [[] for _ in range(len(args.possible_barcodes_files))]
for index, file_name in enumerate(args.possible_barcodes_files):
    with open(file_name, 'r') as file:
        possible_barcodes[index] = file.read().splitlines()

fqs = read_fastq(args.file)
fqs_barcodes = []
for fastq_obj in fqs:
    extracted_sequences = sequence_extraction(fastq_obj.getSeq(), args.sequences, args.max_distances)
    if not extracted_sequences[0]:
        continue
    extracted_barcodes = extract_barcodes(fastq_obj.getSeq(), extracted_sequences, args.lengths)
    if not extracted_barcodes[0]:
        continue
    for i in range(len(extracted_barcodes)-1):
        if args.trim_barcodes[i] != 0:
            extracted_barcodes[i+1] = extracted_barcodes[i+1][:args.trim_barcodes[i]]
    barcode_string = barcodes_validation(extracted_barcodes[1:], possible_barcodes)
    if barcode_string == "ERROR":
        continue
    # else:
    #     barcode_string = str(extracted_barcodes[0]) + "_" + barcode_string
    new_header = "@[" + barcode_string + "]_" + str(fastq_obj.getHead()).replace("@", "")
    fqs_barcodes.append(fq.fastq_object(new_header, fastq_obj.getSeq(), fastq_obj.getQual()))


write_fastq(fqs_barcodes, args.file.replace("portion", "barcodes").replace("tmp", "tmp/tmp_fastq"))
    

