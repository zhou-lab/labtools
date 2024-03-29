from fuzzysearch import find_near_matches
import argparse
import math
import fastq as fq
import os

def write_fastq(foss, file):
    for fos in foss:
        fq.write(fos, file, mode="a")

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

def trimming(fos, seqs, max_lens, output):

    for fastq_obj in fos:
        fos_trimmed = []
        length = len(fastq_obj.getSeq())
        start_coors = -1
        end_coors = math.inf

        for (seq, max_len) in zip(seqs, max_lens):
            ad_matches = find_near_matches(seq, fastq_obj.getSeq(), max_l_dist=max_len)
            if len(ad_matches) > 0:
                if ad_matches[0].end - 0 > length - ad_matches[0].start:
                    if ad_matches[0].start < end_coors:
                        end_coors = ad_matches[0].start
                else:
                    if ad_matches[0].end > start_coors:
                        start_coors = ad_matches[0].end
    
        if start_coors == -1 and end_coors == math.inf:
            fos_trimmed.append(fastq_obj)
        elif start_coors == -1:
            head_trimmed = fastq_obj.getHead()[:37]
            seq_trimmed = fastq_obj.getSeq()[:end_coors]
            qual_trimmed = fastq_obj.getQual()[:end_coors]
            fos_trimmed.append(fq.fastq_object(head_trimmed, seq_trimmed, qual_trimmed))
        elif end_coors == math.inf:
            head_trimmed = fastq_obj.getHead()[:37]
            seq_trimmed = fastq_obj.getSeq()[start_coors:]
            qual_trimmed = fastq_obj.getQual()[start_coors:]
            fos_trimmed.append(fq.fastq_object(head_trimmed, seq_trimmed, qual_trimmed))
        elif end_coors - start_coors > 10:
            head_trimmed = fastq_obj.getHead()[:37]
            seq_trimmed = fastq_obj.getSeq()[start_coors:end_coors]
            qual_trimmed = fastq_obj.getQual()[start_coors:end_coors]
            fos_trimmed.append(fq.fastq_object(head_trimmed, seq_trimmed, qual_trimmed))
        write_fastq(fos_trimmed, output)

parser = argparse.ArgumentParser(description="This is a command line tool trim off any adapter/barcode/sequence you don't want. An example of running this tool: python3 targeted_trimadapters.py -d /path/to/directory -s TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG CTGTCTCTTATACACATCTGACGCTGCCGACGA GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -m 3 3 3 3 -o /path/to/directory/name.fastq")
group = parser.add_mutually_exclusive_group()
group.add_argument("-d", "--directory", type=str, help="The folder that includes all fastq files.")
group.add_argument("-f", "--file", type=str, help="The file that is a fastq.")
parser.add_argument("-s", "--sequences", type=str, nargs='+', help="The sequences that should be contained.")
parser.add_argument("-m", "--max_distances", type=int, nargs='+', help="The maximum edit distance.")
parser.add_argument("-o", "--output", type=str, help="The file for storing trimmed reads.")
args = parser.parse_args()

if len(args.sequences) != len(args.max_distances):
    print("The sequences and the maximum distances do not match!")
    print(str(len(args.sequences)) + " does not equal to " + str(len(args.max_distances)))
    exit()

if args.directory:
    for filename in os.listdir(fastq_folder):
        fos = read_fastq_file(fastq_folder, filename)
        trimming(fos, args.sequences, args.max_distances, args.output)
else:
    fos = read_fastq(args.file)
    trimming(fos, args.sequences, args.max_distances, args.output)

    

        

