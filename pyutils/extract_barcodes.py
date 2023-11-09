#!/usr/bin/env python3
import argparse
from sequence_operations_common import reverse_complement

read_direction = ["B", "F"]

parser = argparse.ArgumentParser(description="Extract your barcodes.")
parser.add_argument("-r", "--read", type=str, help="The enter read.")
parser.add_argument("-s", "--starts", type=int, nargs='+', help="Start positions of linkers.")
parser.add_argument("-e", "--ends", type=int, nargs='+', help="Ending positions of linkers.")
parser.add_argument("-l", "--lengths", type=int, nargs='+', help="Lengths of the barcodes not bewteen the linkers.")
parser.add_argument("--direction", choices=read_direction, help="Select one of the barcode directions")
args = parser.parse_args()

if len(args.starts) != len(args.ends):
    print("ERROR: The number of linker starting coordinates does not match the number of linker ending coordinates!")
    exit()

num_linkers = len(args.starts)
barcodes = []
if args.lengths[1] != 0:
    if args.direction == "B":
        if args.starts[num_linkers-1] <= args.lengths[0]:
            print("FAILED")
            exit()
        else: 
            barcodes.append(args.read[args.starts[num_linkers-1]-args.lengths[1]:args.starts[num_linkers-1]])
if args.lengths[0] != 0:
    if args.direction == "F":
        if args.starts[0] <= args.lengths[1]:
            print("FAILED")
            exit()
        else:
            barcodes.append(args.read[args.starts[0]-args.lengths[0]:args.starts[0]])
for i in range(num_linkers):
    if args.direction == "F":
        if i < num_linkers-1:
            barcodes.append(args.read[args.ends[i]:args.starts[i+1]])
    elif args.direction == "B":
        if i < num_linkers-1:
            barcodes.append(args.read[args.ends[num_linkers-1-i]:args.starts[num_linkers-2-i]])
if args.lengths[1] != 0:
    if args.direction == "F":
        if len(args.read)-args.ends[num_linkers-1] < args.lengths[0]:
            print("FAILED")
            exit()
        else: 
            barcodes.append(args.read[args.ends[num_linkers-1]:args.ends[num_linkers-1]+args.lengths[1]])
if args.lengths[0] != 0:
    if args.direction == "B":
        if len(args.read)-args.ends[0] <= args.lengths[1]:
            print("FAILED")
            exit()
        else:
            barcodes.append(args.read[args.ends[0]:args.ends[0]+args.lengths[0]])
print(args.direction, end='\t')
if args.direction == "F":
    for barcode in barcodes:
        print(barcode, end='\t')
elif args.direction == "B":
    reverse = barcodes[::-1]
    for barcode in reverse:
        print(reverse_complement(barcode), end='\t')
print()


