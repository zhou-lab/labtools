#!/usr/bin/env python3
import sys

def reverse_complement(dna_sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    reversed_sequence = dna_sequence[::-1]
    reverse_complementary_sequence = [complement_dict[base] for base in reversed_sequence]
    reverse_complementary_sequence = ''.join(reverse_complementary_sequence)
    return reverse_complementary_sequence

if __name__ == "__main__":
    print(reverse_complement(sys.argv[1]))