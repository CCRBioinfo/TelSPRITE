import argparse
from collections import OrderedDict
import assembly_methods as am
from tel_table import Tel_table

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
    '-b', 
    '--bam_path',
    required=True
    )
parser.add_argument(
    '-a', 
    '--assembly',
    required=True
    )
parser.add_argument(
    '-r', 
    '--ref_path',
    required=True
    )
parser.add_argument(
    '-c', 
    '--chr_end_prox_threshold',
    type=int,
    required=True
    )
parser.add_argument(
    '-t', 
    '--locus_qual_threshold',
    type=int,
    required=True
    )
parser.add_argument(
    '-f', 
    '--locus_qual_frac',
    type=float,
    required=True
    )
parser.add_argument(
    '-s', 
    '--start_similarity_threshold',
    type=int,
    required=True
    )
args = parser.parse_args()

# Initialize reference assembly
assembly = am.get_assembly(args.assembly)
assembly.load_ref_seq(args.ref_path)

# Create organized dictionary of only valid neotelomere sites
condensed = OrderedDict([(chr, []) for chr in assembly.chr_names])

with open(args.input) as tel_table:
    for line in tel_table:
        row = line.split()
        chr = row[0]
        start = int(row[1])
        if assembly.valid_locus((chr, start), args.chr_end_prox_threshold, args.bam_path, args.locus_qual_threshold, args.locus_qual_frac):
            condensed[chr].append(row[:4])

for chr in condensed:
    condensed[chr] = sorted(condensed[chr], key=lambda x: int(x[1]))

# Combine loci that are the same or soft-clipped just a few base pairs apart
for chr in condensed:
    loci = condensed[chr]
    if len(loci) == 0:
        continue
    combined = []
    
    for locus in loci:#[1:]:
        chr, start, dir, orient = locus
        
        # Find most recent item in list that has the same direction and orientation
        combined_matching = [x for x in combined if (x[2] == dir and x[3] == orient)]
        if len(combined_matching) == 0:
            combined.append([chr, [int(start)], dir, orient])
        else:
            i = combined.index(combined_matching[-1])
            start = int(start)
            combine = False
            
            # Determine if start similarity threshold is met for combining two loci
            for start_pt in combined[i][1]:
                if start-start_pt <= args.start_similarity_threshold:
                    combine = True
            if combine:
                combined[-1][1].append(start)
            else:
                combined.append([chr, [int(start)], dir, orient])
    
    # Add column for number of reads and convert start column back to a single number
    i = 0
    while i < len(combined):
        locus = combined[i]
        locus.append(len(locus[1]))
        
        # Determine most common start point
        start_pts = list(set(locus[1]))
        most_common = start_pts[0]
        for start_pt in start_pts:
            if locus[1].count(start_pt) > locus[1].count(most_common):
                most_common = start_pt
        
        locus[1] = most_common
        combined[i] = [str(x) for x in locus]
        i += 1
    
    condensed[chr] = combined

tels = Tel_table()
tels.col_names = ['Chromosome', 'Start position', 'Direction', 'Orientation', 'Read count']
for chr in condensed:
    for new_tel in condensed[chr]:
        tels.add_tel(new_tel)

tels.write_tel_table(args.output)
