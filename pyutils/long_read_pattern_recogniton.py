#!/usr/bin/env python3
from sequence_operations_common import reverse_complement
import fastq as fq
import argparse
import os
import multiprocessing
import numpy as np
from collections import Counter
import statistics


def sort_key_dist(obj):
    return obj.dist

# This function read the patten file and record the information. 
# Some illegal cases will be eliminated.
def parse_pattern(file_name):
    patterns = []
    known_sequences = []
    with open(file_name, 'r') as file:
        lines = file.readlines()
    for index, line in enumerate(lines):
        parsed_line = line.replace("\n", "").split("\t")
        if len(parsed_line) < 2:
            print("Format error in line " + str(index+1))
            exit()
        else:
            if not any(char.isdigit() for char in parsed_line[0]):
                known_sequences.append(parsed_line)
                patterns.append(parsed_line)
            else:
                parsed_nums = parsed_line[0].split("-")
                if len(parsed_nums) == 2 and parsed_nums[0].isnumeric() and parsed_nums[1].isnumeric():
                    patterns.append([int(parsed_nums[0]), int(parsed_nums[1]), parsed_line[1]])
                else: 
                    print("Format error in line " + str(index+1))
                    exit()
    return patterns, known_sequences

# This function tries to pinpoint the known sequences (mostly are linkers, adapters, etc.), 
# and return their match objects from fuzzy search.
def extract_known_sequences(read, sequences, max_distances):
    contain_all = True
    extracted = []
    for (seq, max_len) in zip(sequences, max_distances):
        if not len(check_one_read_standalone(read, seq, minimum_score=max_len)) > 0:
            contain_all = False
    if contain_all:
        forward = []
        forward_check = True
        for (seq, max_len) in zip(sequences, max_distances):
            forward.append(sorted(check_one_read_standalone(read, seq, minimum_score=max_len), key=sort_key_dist)[0])
        for i in range(len(forward)):
            if i < len(forward)-1:
                if forward[i].end > forward[i+1].start:
                    forward_check = False
                    break
        if not forward_check:
            extracted = [False]
        else:
            extracted = forward
    else:
        extracted = [False]
    return extracted

# This function uses the known sequences coordinates to locate and extract the sequences with ranges and lengths. 
# It will return True if the read is good, False otherwise. It will also call the save function to save the good reads.
def extract_unknown_sequences(read, patterns, extracted_known_sequence, ro):
    num_linkers = len(extracted_known_sequence)
    sequences = []
    quals = []
    coordinates = []
    sequences.append(read[:extracted_known_sequence[0].start])
    quals.append(ro.getQual()[:extracted_known_sequence[0].start])
    coordinates.append([0,extracted_known_sequence[0].start])
    for i in range(num_linkers):
        if i < num_linkers-1:
            sequences.append(read[extracted_known_sequence[i].end:extracted_known_sequence[i+1].start])
            quals.append(ro.getQual()[extracted_known_sequence[i].end:extracted_known_sequence[i+1].start])
            coordinates.append([extracted_known_sequence[i].end,extracted_known_sequence[i+1].start])
    sequences.append(read[extracted_known_sequence[num_linkers-1].end:])
    quals.append(ro.getQual()[extracted_known_sequence[num_linkers-1].end:])
    coordinates.append([extracted_known_sequence[num_linkers-1].end,len(read)])
    patterns_correlated = []
    patterns_tmp = []
    for pattern in patterns:
        if len(pattern) == 2 and patterns_tmp == []:
            continue
        if len(pattern) == 2:
            patterns_correlated.append(patterns_tmp[:])
            patterns_tmp = []
        else:
            patterns_tmp.append(pattern)
    if patterns_tmp != []:
        patterns_correlated.append(patterns_tmp[:])
    for index, patterns_correlated_single in enumerate(patterns_correlated):
        if len(patterns_correlated_single) > 1:
            difference = 0
            for pattern_correlated_single in patterns_correlated_single:
                if pattern_correlated_single[0] != pattern_correlated_single[1] or (pattern_correlated_single[0] == 0 and pattern_correlated_single[1] == 0):
                    difference = difference + 1
                elif pattern_correlated_single[1] < pattern_correlated_single[0]:
                    print("Range error.")
                    exit()
            if difference > 1:
                print("Format error: there cannot have 2 ranges in between 2 fuzzy searches (come up with better message here).")
                exit()
        ranges = [0,0]
        jump = False
        for patterns_correlated_single in patterns_correlated_single:
            if patterns_correlated_single[0] == 0 and patterns_correlated_single[1] == 0:
                jump = True
            ranges[0] = ranges[0] + patterns_correlated_single[0]
            ranges[1] = ranges[1] + patterns_correlated_single[1]
        if (len(sequences[index]) > ranges[1] or len(sequences[index]) < ranges[0]) and not jump:
            return None, False
    contain_genome = False
    contain_barcode = False
    for pattern in patterns:
        if len(pattern) > 2:
            if "barcode" in pattern[2].lower():
                contain_barcode = True
            if "genome" in pattern[2].lower():
                contain_genome = True
    if contain_barcode and contain_genome:
        ro = config_trimmed_good_read(ro.getHead(), sequences, quals, coordinates, patterns_correlated, ro)
        if ro is None:
            return None, False
    return ro, True

# This function trim the good read, extract the barcodes, and config it for saving. 
def config_trimmed_good_read(header, sequences, quals, coordinates, patterns_correlated, ro):
    barcodes = []
    genome_seq = ""
    genome_qual = ""
    genome_coor = []

    for sequence, qual, patterns_correlated_single, coordinate_pair in zip(sequences, quals, patterns_correlated, coordinates):

        front_fixed_patterns = []
        back_fixed_patterns = []

        front_fixed_sequences = []
        front_fixed_quals = []
        front_fixed_coor = []

        back_fixed_sequences = []
        back_fixed_quals = []
        back_fixed_coor = []

        seqs_cor_patterns = []
        quals_cor_patterns = []
        coor_cor_patterns = []

        got_non_fixed = False
        for pattern_correlated_single in patterns_correlated_single:
            if pattern_correlated_single[0] != pattern_correlated_single[1] or (pattern_correlated_single[0] == 0 and pattern_correlated_single[1] == 0):
                got_non_fixed = True
            else:
                if got_non_fixed:
                    back_fixed_patterns.insert(0, pattern_correlated_single)
                else:
                    front_fixed_patterns.append(pattern_correlated_single)

        remaining_seq = sequence
        remaining_qual = qual
        remaining_coor = coordinate_pair[:]
        for front_fixed_pattern in front_fixed_patterns:
            front_fixed_sequences.append(remaining_seq[:front_fixed_pattern[0]])
            front_fixed_quals.append(remaining_qual[:front_fixed_pattern[0]])
            front_fixed_coor.append([remaining_coor[0], remaining_coor[0]+front_fixed_pattern[0]-1])
            remaining_seq = remaining_seq[front_fixed_pattern[0]:]
            remaining_qual = remaining_qual[front_fixed_pattern[0]:]
            remaining_coor = [remaining_coor[0]+front_fixed_pattern[0], remaining_coor[1]]
        
        for back_fixed_pattern in back_fixed_patterns:
            back_fixed_sequences.insert(0,remaining_seq[-back_fixed_pattern[0]:])
            back_fixed_quals.insert(0,remaining_qual[-back_fixed_pattern[0]:])
            back_fixed_coor.insert(0,[remaining_coor[1]-back_fixed_pattern[0]+1, remaining_coor[1]])
            remaining_seq = remaining_seq[:-back_fixed_pattern[0]]
            remaining_qual = remaining_qual[:-back_fixed_pattern[0]]
            remaining_coor = [remaining_coor[0], remaining_coor[1]-back_fixed_pattern[0]]

        seqs_cor_patterns.extend(front_fixed_sequences)
        if got_non_fixed:
            seqs_cor_patterns.append(remaining_seq)
        seqs_cor_patterns.extend(back_fixed_sequences)
        quals_cor_patterns.extend(front_fixed_quals)
        if got_non_fixed:
            quals_cor_patterns.append(remaining_qual)
        quals_cor_patterns.extend(back_fixed_quals)
        coor_cor_patterns.extend(front_fixed_coor)
        if got_non_fixed:
            coor_cor_patterns.append(remaining_coor)
        coor_cor_patterns.extend(back_fixed_coor)

        for seq_cor_pattern, qual_cor_pattern, coor_cor_pattern, pattern_correlated_single in zip(seqs_cor_patterns, quals_cor_patterns, coor_cor_patterns, patterns_correlated_single):
            if "barcode" in pattern_correlated_single[2].lower():
                barcodes.append(seq_cor_pattern)
            elif "genome" in pattern_correlated_single[2].lower():
                genome_seq = seq_cor_pattern
                genome_qual = qual_cor_pattern
                genome_coor = coor_cor_pattern
    
    if genome_seq != "":
        if args.format == 'fastq':
            combined_barcodes = '_'.join(map(str, barcodes))
            new_header = "@[" + combined_barcodes + "]_" + str(header).replace("@", "")
            return fq.fastq_object(new_header, genome_seq, genome_qual)
        else:
            combined_barcodes = '-'.join(map(str, barcodes))
            ro.sam_obj.set_tag("CR", combined_barcodes)
            ro.trimming(genome_coor[0], genome_coor[1])
            ro.alignment_eraser()
            return ro


# This function is the workflow for one entire file. 
def process_one_file(file):
    read_sequences = []
    if args.format == 'fastq':
        ros = read_fastq(file)
    elif args.format == 'sam':
        ros, header = read_sam_edit(file, mode="r")
        ros = sam_filter(ros)
    elif args.format == 'bam':
        ros, header = read_sam_edit(file, mode="rb")
        ros = sam_filter(ros)

    for index, ro in enumerate(ros):
        read_sequences.append(ro.getSeq())
    
    results = []
    ros_for_saves = []
    for pattern_files in args.patterns:
        pattern, known_sequences = parse_pattern(pattern_files)
        result = []
        ros_for_save = []
        for index, read_seq in enumerate(read_sequences):
            sequences = []
            max_dist = []
            for combo in known_sequences:
                sequences.append(combo[0])
                try: 
                    max_dist.append(float(combo[1]))
                except:
                    max_dist.append(0.8)
            extracted = extract_known_sequences(read_seq, sequences, max_dist)
            if not extracted[0]:
                result.append(False)
                ros_for_save.append(None)
                continue
            ro, whether = extract_unknown_sequences(read_seq, pattern, extracted, ros[index])
            ros_for_save.append(ro)
            result.append(whether)
        ros_for_saves.append(ros_for_save[:])
        results.append(result[:])
    num_patterns = [0] * len(args.patterns)
    save_patterns = [[] for _ in range(len(args.patterns))]
    save_multiple = []
    multiple = 0
    none = 0
    for i in range(len(results[0])):
        matches = 0
        for j in range(len(args.patterns)):
            if results[j][i]:
                matches = matches + 1
        if matches == 0:
            none = none + 1
        elif matches > 1:
            multiple = multiple + 1
            save_multiple.append(ros[i])
        else: 
            for j in range(len(args.patterns)):
                if results[j][i]:
                    num_patterns[j] = num_patterns[j] + 1
                    save_patterns[j].append(ros_for_saves[j][i])
                    break

    good = []
    for index, pattern_name in enumerate(pattern_names):
        if "forward" in pattern_name.lower():
            good.extend(save_patterns[index])
        elif "reverse" in pattern_name.lower():
            save_pattern_reversing = []
            for jndex in range(len(save_patterns[index])):
                ro = save_patterns[index][jndex]
                if args.format == 'fastq':
                    old_header = ro.getHead()
                    start_index = old_header.find("[")
                    end_index = old_header.find("]")
                    barcodes = old_header[start_index + 1 : end_index]
                    barcodes_list = barcodes.split("_")
                    barcodes_list.reverse()
                    for i in range(len(barcodes_list)):
                        barcodes_list[i] = reverse_complement(barcodes_list[i])
                    combined_barcodes = '_'.join(map(str, barcodes_list))
                    new_header = old_header.replace(barcodes, combined_barcodes)
                    save_patterns[index][jndex] = fq.fastq_object(new_header, ro.getSeq(), ro.getQual())
                else:
                    old_barcodes = ro.sam_obj.get_tag("CR")
                    barcodes_list = old_barcodes.split("-")
                    barcodes_list.reverse()
                    for i in range(len(barcodes_list)):
                        barcodes_list[i] = reverse_complement(barcodes_list[i])
                    combined_barcodes = '-'.join(map(str, barcodes_list))
                    ro.sam_obj.set_tag("CR", combined_barcodes)
                    save_patterns[index][jndex] = ro
            good.extend(save_patterns[index])
    if args.format == 'fastq':
        write_fastq(good, args.output)
    elif args.format == 'sam':
        new_header = sam_header_dealignment(header.to_dict())
        write_sam_edit(good, args.output, new_header, mode="w")
    elif args.format == 'bam':
        new_header = sam_header_dealignment(header.to_dict())
        write_sam_edit(good, args.output, new_header, mode="wb")
    return num_patterns, multiple, none

def QC(file):
    read_sequences = []
    QC_data = dict()
    if args.format == 'fastq':
        ros = read_fastq(file)
    elif args.format == 'sam':
        ros, header = read_sam_edit(file, mode="r")
    elif args.format == 'bam':
        ros, header = read_sam_edit(file, mode="rb")

    bases = 0
    GC = 0
    longest = 0
    lengths = []
    for ro in ros:
        seq = ro.getSeq().upper()
        length = len(seq)
        read_sequences.append(seq)
        bases = bases + length
        lengths.append(length)
        for char in seq:
            if char == "G" or char == "C":
                GC = GC + 1

    half_length = bases / 2
    lengths.sort(reverse=True)
    cumulative_length = 0
    n50 = 0
    for contig_length in lengths:
        cumulative_length += contig_length
        if cumulative_length >= half_length:
            n50 = contig_length
            break

    QC_data["Total reads"] = len(ros)
    QC_data["Total bases"] = bases
    QC_data["N50"] = n50
    QC_data["Mean read length"] = bases / QC_data["Total reads"]
    QC_data["Median read length"] = statistics.median(lengths)
    QC_data["GC content"] = GC / bases

    edges = [0]
    max_length = max(lengths)
    while edges[-1] < max_length:
        edges.append(edges[-1]+1)
    hist, edges = np.histogram(lengths, bins=edges)
    QC_data["Trimmed size hist"] = {}
    for index, hist_s in enumerate(hist):
        QC_data["Trimmed size hist"][str(edges[index])] = hist_s

    if args.format == 'fastq':
        header = ros[0].getHead().split(" ")[0]
        start_index = header.find("[")
        end_index = header.find("]")
        if start_index != -1 and end_index != -1:
            contain_barcodes = True
        else:
            contain_barcodes = False
    else:
        contain_barcodes = True
        try:
            ros[0].sam_obj.get_tag("CR")
        except:
            contain_barcodes = False
    
    if contain_barcodes:
        barcodes_total = []
        reads_decipher = []
        for ro in ros:
            if args.format == 'fastq':
                header = ro.getHead().split(" ")[0]
                start_index = header.find("[")
                end_index = header.find("]")
                barcodes = header[start_index + 1 : end_index]
                barcodes_list = barcodes.split("_")
            else:
                barcodes = ro.sam_obj.get_tag("CR")
                barcodes_list = barcodes.split("-")

            reads_decipher.append(barcodes)

            barcodes_total.append(barcodes_list)

        string_counts = Counter(reads_decipher)
        string_counts_new = {}
        for key, value in string_counts.items():
            if value == 1:
                continue
            string_counts_new[key] = value
        sorted_strings = sorted(string_counts_new.items(), key=lambda x: x[1], reverse=True)
        counts = dict(sorted_strings)
        key = "Unique barcode combination counts"
        QC_data[key] = {}
        for sorted_string, count in counts.items():
            QC_data[key][sorted_string] = count
    
        num_barcodes = len(barcodes_total[0])
        list_of_barcodes = [[] for _ in range(num_barcodes)]
        for barcodes in barcodes_total:
            for i in range(num_barcodes):
                list_of_barcodes[i].append(barcodes[i])

        for i in range(num_barcodes):
            string_counts = Counter(list_of_barcodes[i])
            sorted_strings = sorted(string_counts, key=lambda x: string_counts[x], reverse=True)
            counts = []
            for sorted_string in sorted_strings:
                counts.append(string_counts[sorted_string])
            key = "barcode " + str(i+1) + " counts"
            QC_data[key] = {}
            j = 0
            for sorted_string, count in zip(sorted_strings, counts):
                QC_data[key][sorted_string] = count
                j = j + 1
                if j >= 200:
                    break

        string_counts = Counter(reads_decipher)
        sorted_strings = sorted(string_counts, key=lambda x: string_counts[x], reverse=True)
        elements = []
        for sorted_string in sorted_strings:
            elements.append(string_counts[sorted_string])
        edges = [0]
        max_ele = max(elements)
        while edges[-1] < max_ele:
            edges.append(edges[-1]+1)
        hist, edges = np.histogram(elements, bins=edges)
        QC_data["Barcode density hist"] = {}
        for index, hist_s in enumerate(hist):
            QC_data["Barcode density hist"][str(edges[index])] = hist_s
    return QC_data


parser = argparse.ArgumentParser(description="This is a command line tool designed to parse long read sequencing data based on the pattern you provide.")
parser.add_argument('--format', choices=['fastq', 'sam', 'bam'], help='Choose one of the options.')
parser.add_argument('--mode', choices=['alignment', 'fuzzy'], help='Choose one of the options.')
group = parser.add_mutually_exclusive_group()
group.add_argument("-i", "--input", type=str, help="The input file that contains raw data (basecalled).")
group.add_argument("-d", "--directory", type=str, help="The input directory that contains raw data (basecalled).")
parser.add_argument("-p", "--patterns", type=str, nargs='+', help="The pattern files.")
parser.add_argument("-o", "--output", type=str, help="The output directory for the trimmed data plus the barcodes.")
args = parser.parse_args()

if args.mode == 'alignment':
    from sequence_operations_common import check_one_read_standalone_alignment as check_one_read_standalone
elif args.mode == 'fuzzy':
    from sequence_operations_common import check_one_read_standalone_fuzzy as check_one_read_standalone

if args.format == 'fastq':
    from sequence_operations_common import read_fastq, write_fastq
elif args.format == 'sam' or args.format == 'bam':
    from sequence_operations_common import read_sam_edit, write_sam_edit, sam_filter, sam_header_dealignment

pattern_names = []
for pattern in args.patterns:
    pattern_file_name = pattern.split("/")[-1]
    pattern_name = pattern_file_name.split(".")[0]
    pattern_names.append(pattern_name)

if args.input is not None:
    num_patterns, multiple, none = process_one_file(args.input)
elif args.directory is not None:
    files = []
    for filename in os.listdir(args.directory):
        if args.format == 'fastq':
            if not filename.startswith('.') and (filename.endswith('fastq') or filename.endswith('fq') or filename.endswith('fastq.gz') or filename.endswith('fq.gz')):
                f = os.path.join(args.directory, filename)
                if os.path.isfile(f):
                    files.append(f)
        else:
            if not filename.startswith('.') and filename.endswith(args.format):
                f = os.path.join(args.directory, filename)
                if os.path.isfile(f):
                    files.append(f)
    with multiprocessing.Pool() as pool:
        results = pool.imap(process_one_file, files)
        pool.close()
        pool.join()
    num_patterns = [0] * len(args.patterns)
    multiple = 0
    none = 0
    for num_patterns_s, multiple_s, none_s in results:
        for i in range(len(num_patterns)):
            num_patterns[i] = num_patterns[i] + num_patterns_s[i]
        multiple = multiple + multiple_s
        none = none + none_s

fastq_file = args.output
log_content = QC(fastq_file)
if args.format == 'fastq':
    log_file = args.output.replace(".fastq", "") + "_patterns_log.txt"
elif args.format == 'sam':
    log_file = args.output.replace(".sam", "") + "_patterns_log.txt"
elif args.format == 'bam':
    log_file = args.output.replace(".bam", "") + "_patterns_log.txt"
with open(log_file, 'a') as file:
    for key, value in log_content.items():
        if type(value) is dict:
            file.write('\n')
            file.write(f'{key}:\n')
            for key_s, value_s in value.items():
                file.write(f'{key_s}: {value_s}\n')
        else:
            file.write(f'{key}: {value}\n')

    file.write('\n')
    for index, pattern in enumerate(args.patterns):
        file.write(pattern.replace('.txt', '') + ": " + str(num_patterns[index]) + "\n")
    file.write("Multiple matches: " + str(multiple) + "\n")
    file.write("No matches: " + str(none) + "\n")


    