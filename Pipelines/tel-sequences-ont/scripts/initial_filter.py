import pysam
import argparse
import count_tel_repeats as ct

__author__ = "David Wilson"

parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', 
    '--input',
    required=True
    )
parser.add_argument(
    '-o', 
    '--output',
    required=True
    )
parser.add_argument(
    '-w', 
    '--window_size',
    required=True,
    type=int,
    default='300'
    )
parser.add_argument(
    '-t',
    '--threshold',
    required=True,
    type=float,
    default='0.9'
    )
parser.add_argument(
    '-m',
    '--mapq_threshold',
    required=True,
    type=int
    )
args = parser.parse_args()

repeat_threshold = int(args.window_size*args.threshold/6)

# Run an initial filter for reads that meet the mapq minimum and have enough repeats
# to have a chance of having at least one telomere window

with pysam.AlignmentFile(args.input, "rb") as in_file:
    with pysam.AlignmentFile(args.output, "wb", template=in_file) as out_file:
        for read in in_file.fetch(until_eof=True):
            if read.mapping_quality >= args.mapq_threshold:
                seq = read.query_sequence
                if ct.seq_is_tel(seq, ct.g_repeats, ct.c_repeats, ct.g_shared, ct.c_shared, repeat_threshold):
                    out_file.write(read)

