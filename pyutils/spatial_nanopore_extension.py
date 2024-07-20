import argparse
import fastq as fq
from sequence_operations_common import reverse_complement

def main(fastq, output):
    fos = fq.read(fastq) # Iterator of fastq entries.
    fos = list(fos)

    fos_result = []

    for fo in fos:
        input_string = fo.getHead()
        start_index = input_string.index('[') + 1
        end_index = input_string.index(']')
        barcodes_string = input_string[start_index:end_index]
        header = input_string[end_index + 1:]
        if not barcodes_string:
            print(f"Error: {fo.getHead()} doesn't have barcodes.")
            exit()
        barcodes = barcodes_string.split("_")
        if len(barcodes) < 2:
            print(f"Error: {fo.getHead()} has the wrong amount of barcodes.")
            exit()
        if len(barcodes[0]) == 8 and len(barcodes[1]) == 8:
            barcodes_ready = reverse_complement(barcodes[0] + barcodes[1])
            new_header = "@" + barcodes_ready + header
            fos_result.append(fq.fastq_object(new_header, fo.getSeq(), fo.getQual()))

    fq.write(fos_result, output)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--fastq', help='The fastq file that are trimmed and barcode extracted using pattern recognition.', dest='fastq')
    parser.add_argument('--output', help='The output file.', dest='output')

    args = parser.parse_args()

    main(args.fastq, args.output)