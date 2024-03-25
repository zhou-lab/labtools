#!/usr/bin/env python3
import sys
import fastq as fq
from fuzzysearch import find_near_matches
import pysam
from attr import attrs, attrib
from Levenshtein import distance
from Bio import Align
import math

@attrs(frozen=True, slots=True)
class Match:
    start = attrib(type=int, eq=True, hash=True)
    end = attrib(type=int, eq=True, hash=True)
    dist = attrib(type=int, eq=True, hash=True)
    matched = attrib(eq=False, hash=False)

def check_one_read_standalone_alignment(text, target_substring, minimum_score=0.8, match_score = 1, mismatch_score = -1, gap_score = -1):
    matches = []

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = match_score
    aligner.mismatch_score = mismatch_score
    aligner.gap_score = gap_score

    alignments = aligner.align(text, target_substring)
    for alignment in sorted(alignments):
        percentage = alignment.score / (len(target_substring) * match_score)
        start = alignment.coordinates[0][0]
        end = alignment.coordinates[0][-1]
        if percentage > minimum_score:
            matches.append(Match(start, end, percentage, text[start:end]))
    return matches

def write_fastq(foss, file):
    for fos in foss:
        fq.write(fos, file, mode="a")

def read_fastq(file):
    fos = fq.read(file)
    fos = list(fos)
    return fos

def sam_filter(reads):
    good_reads = []
    for read in reads:
        if read.getSeq() is not None:
            good_reads.append(read)
    return good_reads

def read_sam_edit(file, check_sq=False, mode="rb"):
    from sam_editor import sam_editor
    reads = []
    bamfile = pysam.AlignmentFile(file, mode, check_sq=check_sq)
    header = bamfile.header
    for read in bamfile:
        reads.append(sam_editor(read))
    bamfile.close()
    return reads, header

def write_sam_edit(reads, file, header, mode="wb"):
    output_bamfile = pysam.AlignmentFile(file, mode, header=header)
    for read in reads:
        output_bamfile.write(read.sam_obj)
    output_bamfile.close()

def sam_header_dealignment(header):
    header['SQ'] = []
    return header

def check_one_read_standalone_fuzzy(read, seq, minimum_score=0.8):
    max_l_dist = int(math.ceil(len(seq) * (1 - minimum_score)))
    return find_near_matches(seq, read, max_l_dist=max_l_dist)

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = dna_sequence[::-1]
    reverse_complementary_sequence = [complement_dict[base] for base in reversed_sequence]
    reverse_complementary_sequence = ''.join(reverse_complementary_sequence)
    return reverse_complementary_sequence



if __name__ == "__main__":
    print(reverse_complement(sys.argv[1]))