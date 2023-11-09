from fuzzysearch import find_near_matches
import fastq as fq
import os
import argparse

def read_fastq(file):
    fos = fq.read(file)
    fos = list(fos)
    return fos

def read_fastq_file(folder, filename):
    foss = []
    if not filename.startswith('.') and filename.endswith('fastq'):
        f = os.path.join(folder, filename)
        if os.path.isfile(f):
            fos = read_fastq(f)
            for fo in fos:
                foss.append(fo)
    return foss

def reverse_complement(sequence):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reverse_comp_sequence = ''.join(complement[base] for base in reversed(sequence))
    return reverse_comp_sequence

def demultiplex(fastq_folder, seqs, output_files):

    max_dis = len(seqs[0]) // 7
    print(max_dis)

    for filename in os.listdir(fastq_folder):

        fos = read_fastq_file(fastq_folder, filename)
        for fastq_obj in fos:
            file_match_num = 0
            file_match = None
            for index, output_file in enumerate(output_files):
                if (len(find_near_matches(seqs[2*index-1], fastq_obj.getSeq(), max_l_dist=max_dis)) > 0 and len(find_near_matches(seqs[2*index], fastq_obj.getSeq(), max_l_dist=max_dis)) > 0) or (len(find_near_matches(reverse_complement(seqs[2*index-1]), fastq_obj.getSeq(), max_l_dist=max_dis)) > 0 and len(find_near_matches(reverse_complement(seqs[2*index]), fastq_obj.getSeq(), max_l_dist=max_dis)) > 0):
                    file_match_num = file_match_num + 1
                    file_match = output_file
            if file_match_num == 1 and file_match != None:
                fq.write(fastq_obj, file_match, mode="a")
            elif file_match_num > 1:
                fq.write(fastq_obj, 'multiple_matches.fastq', mode="a")
            elif file_match_num == 0:
                fq.write(fastq_obj, 'no_matches.fastq', mode="a")

parser = argparse.ArgumentParser(description="This is a program to demultiplex the nanopore output data. So far, it only supports dual indices.")
parser.add_argument("-f", "--folder", type=str, help="The folder that includes all fastq files.")
parser.add_argument("-s", "--sequences", type=str, nargs='+', help="The sequences that should be contained.")
parser.add_argument("-o", "--output_files", type=str, nargs='+', help="The files for storing trimmed reads.")
args = parser.parse_args()

if len(args.sequences) != 2*len(args.output_files):
    print("The sequences and the output files do not match. For every output files, there should be two sequences.")
    exit()

demultiplex(args.folder, args.sequences, args.output_files)

