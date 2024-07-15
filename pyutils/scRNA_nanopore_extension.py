import argparse

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--sam', help='The bam file that are trimmed and barcode extracted using pattern recognition.', dest='sam')
    parser.add_argument('--output', help='The output file.', dest='output')

    args = parser.parse_args()

    main(args.fastq, args.output)