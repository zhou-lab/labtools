#!/usr/bin/env python3
from sequence_operations_common import reverse_complement
import argparse

def parse_barcodes(direction, barcodes):
    barcodes_read = []
    for barcode in barcodes:
        if isinstance(barcode, float):
            barcodes_read.append("")
        else:
            barcodes_read.append(barcode)
    return barcodes_read

parser = argparse.ArgumentParser(description="This is a command line tool to analyze the extracted barcodes from SPLIT sequencing based single cell datasets sequenced using nanopore.")
parser.add_argument("-f", "--files", type=str, nargs='+', help="The text files that contain barcodes. How many files depends on how many rounds of barcodes.")
parser.add_argument("-b", "--barcodes", type=str, help="The text files that contain barcodes. How many files depends on how many rounds of barcodes.")
args = parser.parse_args()

if args.barcodes == "FAILED":
    print("ERROR")
    exit()

input_array = args.barcodes.strip().split(':')
barcodes = parse_barcodes(input_array[0], input_array[1:])

valid_barcodes = [[] for _ in range(len(args.files))]
for index, file_name in enumerate(args.files):
    with open(file_name, 'r') as file:
        valid_barcodes[index] = file.read().splitlines()

for i in range(len(barcodes)):
    barcode_found = False
    for j in range(len(valid_barcodes[i])):
        if barcodes[i] == valid_barcodes[i][j]:
            barcode_found = True
            break
    if not barcode_found:
        barcodes[i] = "ERROR"

result_string = '_'.join(map(str, barcodes))
if result_string.find("ERROR") == -1:
    print(result_string)
    exit()
print("ERROR")
exit()





