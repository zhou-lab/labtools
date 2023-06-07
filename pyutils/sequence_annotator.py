#!/usr/bin/env python3
import argparse
from fuzzysearch import find_near_matches

parser = argparse.ArgumentParser(description="Annotate a read sequence with provided subsequences and labels.")
parser.add_argument("-r", "--read", type=str, help="The read sequence.")
parser.add_argument("-s", "--sequences", type=str, nargs='+', help="The sequences to find in the read. Format: sequence:label:max_dist")
args = parser.parse_args()

seq_dict = {}

for s in args.sequences:
    seq, label, max_dist = s.split(':')
    seq_dict[seq] = {'label': label, 'max_dist': int(max_dist)}

# Make a copy of the original read for modification
CRED = '\033[91m'
CEND = '\033[0m'

CRED = '\033[91m'
CEND = '\033[0m'

# Collect all matches from all sequences
all_matches = []
for seq, params in seq_dict.items():
    matches = find_near_matches(seq, args.read, max_l_dist=params['max_dist'])
    for match in matches:
        # Add the label to the match
        all_matches.append((match, params["label"]))

# Sort all matches by their start position
all_matches.sort(key=lambda m: (m[0].start, m[0].dist))

read_fragments = []

prev_end = 0
for match, label in all_matches:
    # Add the fragment of the read before the match
    read_fragments.append(args.read[prev_end:match.start])

    # Add the annotated match
    anno = CRED + " {" + label + ":"  + args.read[match.start:match.end] + "|" + str(match.dist) + "} " + CEND
    read_fragments.append(anno)

    prev_end = match.end

# Add the remainder of the read after the last match
read_fragments.append(args.read[prev_end:])

# Combine the fragments back into a single string
read_copy = ''.join(read_fragments)

try:
    print(read_copy)
except BrokenPipeError:
    # Handle broken pipes caused by piping the output (e.g., to 'head')
    sys.stderr.close()  # Close the error pipe
    exit(0)  # Exit with a clean status