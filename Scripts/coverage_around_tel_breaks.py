import argparse
import random
import numpy as np
import pysam
from collections import OrderedDict
import matplotlib.pyplot as plt
import assembly_methods_mod as am

parser = argparse.ArgumentParser()
parser.add_argument(
    '-b', 
    '--bam_file', 
    required=True
    )
parser.add_argument(
    '-f', 
    '--ref_path', 
    required=True
    )
parser.add_argument(
    '-t', 
    '--tel_breaks', 
    required=True
    )
parser.add_argument(
    '-a', 
    '--assembly', 
    required=True
    )
parser.add_argument(
    '-m', 
    '--mapq_min', 
    default=30,
    type=int
    )
parser.add_argument(
    '-s', 
    '--size', 
    required=True,
    type=int
    )
parser.add_argument(
    '-r', 
    '--resolution', 
    required=True,
    type=int
    )
parser.add_argument(
    '-p', 
    '--output_path',
    required=True
    )
parser.add_argument(
    '-o', 
    '--output_type',
    default='data'
    )
args = parser.parse_args()

def coverage_near_locus(bam_file, locus, size, resolution, mapq_min):
    samfile = pysam.AlignmentFile(bam_file, "rb")

    chr, bp, dir = locus
    bin = bp // resolution
    min_bin_i = bin - size
    max_bin_i = bin + size
    coverage = np.zeros(size*2+1)
    #coverage = [0 for _ in range(min_bin_i, max_bin_i+1)]
    
    for read in samfile.fetch(chr, min_bin_i*resolution, (max_bin_i+1)*resolution):
        if read.is_unmapped:
            continue
        
        if read.mapping_quality < mapq_min:
            continue
            
        read_start = read.reference_start
        read_end = read.reference_end
        
        #if read_start < min_bin_i*resolution or read_end >= (max_bin_i+1)*resolution:
        #    continue
        
        start_bin = read_start // resolution
        end_bin = (read_end - 1) // resolution

        # For each bin the read overlaps
        for bin_i in range(start_bin, end_bin + 1):
            bin_start = bin_i * resolution
            bin_end = bin_start + resolution

            # Compute overlap between read and bin
            overlap_start = max(read_start, bin_start)
            overlap_end = min(read_end, bin_end)
            overlap_len = overlap_end - overlap_start
                
            # Add to coverage data
            # Need to subtract min_bin_i to find position in output list
            bin_with_offset = bin_i - min_bin_i
            if bin_with_offset < 0 or bin_with_offset > size*2: 
                # Can occur if a read overlaps but is not entirely contained within the region
                continue
            
            coverage[bin_with_offset] += overlap_len

    coverage /= resolution
    if dir == "left":
        coverage = coverage[::-1]
    return coverage

def coverage_near_multiple_loci(bam_file, loci, size, resolution, mapq_min):
    coverage_all = []
    for locus in loci:
        coverage = coverage_near_locus(bam_file, locus, size, resolution, mapq_min)
        coverage_all.append(coverage)

    coverage_all = np.array(coverage_all)
    column_averages = np.mean(coverage_all, axis=0)
    column_std = np.std(coverage_all, axis=0)
    column_stderr = column_std / np.sqrt(len(loci))
    return (column_averages, column_stderr)

def plot_and_save(breaks, controls, out_path, size, resolution):
    x_values = [(x-size)*resolution for x in range(size*2+1)]
        
    plt.figure(figsize=(8, 6))
    plt.plot(x_values, breaks, label="Actual rearrangements", color='b', marker='o')
    plt.plot(x_values, controls, label="Randomly selected controls", color='r', marker='x')
        
    plt.xlabel('Distance relative to "rearrangement"')
    plt.ylabel('Average sequencing coverage')
    plt.title('Read coverage around telomere rearrangments')
    plt.legend()
    
    plt.savefig(out_path)

def output_data(breaks_data, breaks_stderr, controls_data, controls_stderr, size, resolution, out_path):
    out_str = ''
    
    for n, c, e in zip(['breaks', 'controls'], [breaks_data, controls_data], [breaks_stderr, controls_stderr]):
        x = -size*resolution
        for i in range(len(c)):
            new_line = [str(x), str(c[i]), str(e[i]), n]
            out_str += '\t'.join(new_line) + '\n'
            x += resolution
    
    with open(out_path, "w") as out_file:
        out_file.write(out_str)

class Assembly:
    
    def initialize(self, assembly, ref_path):
        self.am_methods = am.get_assembly(assembly)
        self.am_methods.load_ref_seq(ref_path)
        
        # Parameters
        self.chr_end_prox_threshold = 100000
        self.locus_qual_threshold = 10
        self.locus_qual_frac = 0.4

    
    def valid_locus(self, loc, size, resolution, bam_path):
        # Tests if locus is valid in that the window defined by size is still 
        # within the length of the choromsome
        chr, bp = loc[:2]
        bin = bp // resolution
        if bin - size < 0 or bin + size > self.chr_sizes[chr]//resolution:
            return False
        else:
            if self.am_methods.valid_locus((chr, bp), self.chr_end_prox_threshold, bam_path, self.locus_qual_threshold, self.locus_qual_frac):
                return True
            else:
                return False
    
    def select_random_loci(self, num_loci, include_dir=True):
    # Compute the total length of the genome
        total_genome_length = sum(self.chr_sizes.values())
        
        loci = []  # List to store selected loci
        
        for _ in range(num_loci):
            # Randomly select a locus in the entire genome
            random_position = random.randint(1, total_genome_length)
            
            # Traverse through chromosomes to find where the random position falls
            current_position = 0
            for chromosome, length in self.chr_sizes.items():
                if current_position + length >= random_position:
                    # The random position is within this chromosome
                    base_position = random_position - current_position
                    if include_dir:
                        dir = random.choice(["left", "right"])
                        loci.append((chromosome, base_position, dir))
                    else:
                        loci.append((chromosome, base_position))
                    break
                current_position += length
        
        return loci
    
    def n_random_loci(self, n_target, size, resolution, bam_file):
        loci = []
        n_current = 0
        
        while n_target - n_current > 0:
            new = self.select_random_loci(n_target - n_current)
            new = [loc for loc in new if self.valid_locus(loc, size, resolution, bam_file)]
            loci += new
            n_current = len(loci)
        
        return loci

class Hg19(Assembly):
    
    assembly_name = 'hg19'

    def __init__(self):
        self.chr_sizes = OrderedDict([
            ('chr1', 249250621),
            ('chr2', 243199373),
            ('chr3', 198022430),
            ('chr4', 191154276),
            ('chr5', 180915260),
            ('chr6', 171115067),
            ('chr7', 159138663),
            ('chr8', 146364022),
            ('chr9', 141213431),
            ('chr10', 135534747),
            ('chr11', 135006516),
            ('chr12', 133851895),
            ('chr13', 115169878),
            ('chr14', 107349540),
            ('chr15', 102531392),
            ('chr16', 90354753),
            ('chr17', 81195210),
            ('chr18', 78077248),
            ('chr19', 59128983),
            ('chr20', 63025520),
            ('chr21', 48129895),
            ('chr22', 51304566),
            ('chrX', 155270560),
            ('chrY', 59373566),
            ('chrM', 16571)])
        super(Hg19, self).__init__()

class Hg38(Assembly):
    
    assembly_name = 'hg38'

    def __init__(self):
        self.chr_sizes = OrderedDict([
            ('chr1', 248956422),
            ('chr2', 242193529),
            ('chr3', 198295559),
            ('chr4', 190214555),
            ('chr5', 181538259),
            ('chr6', 170805979),
            ('chr7', 159345973),
            ('chr8', 145138636),
            ('chr9', 138394717),
            ('chr10', 133797422),
            ('chr11', 135086622),
            ('chr12', 133275309),
            ('chr13', 114364328),
            ('chr14', 107043718),
            ('chr15', 101991189),
            ('chr16', 90338345),
            ('chr17', 83257441),
            ('chr18', 80373285),
            ('chr19', 58617616),
            ('chr20', 64444167),
            ('chr21', 46709983),
            ('chr22', 50818468),
            ('chrX', 156040895),
            ('chrY', 57227415),
            ('chrM', 16569)])
        super(Hg38, self).__init__()

class T2T_CHM13v2_0(Assembly):
    
    assembly_name = 't2t_chm13'

    def __init__(self):
        self.chr_sizes = OrderedDict([
            ('NC_060925.1', 248387328),
            ('NC_060926.1', 242696752),
            ('NC_060927.1', 201105948),
            ('NC_060928.1', 193574945),
            ('NC_060929.1', 182045439),
            ('NC_060930.1', 172126628),
            ('NC_060931.1', 160567428),
            ('NC_060932.1', 146259331),
            ('NC_060933.1', 150617247),
            ('NC_060934.1', 134758134),
            ('NC_060935.1', 135127769),
            ('NC_060936.1', 133324548),
            ('NC_060937.1', 113566686),
            ('NC_060938.1', 101161492),
            ('NC_060939.1', 99753195),
            ('NC_060940.1', 96330374),
            ('NC_060941.1', 84276897),
            ('NC_060942.1', 80542538),
            ('NC_060943.1', 61707364),
            ('NC_060944.1', 66210255),
            ('NC_060945.1', 45090682),
            ('NC_060946.1', 51324926),
            ('NC_060947.1', 154259566),
            ('NC_060948.1', 62460029)])
        super(T2T_CHM13v2_0, self).__init__()

def get_assembly(assembly_name):
    for cls in Assembly.__subclasses__():
        if cls.assembly_name == assembly_name:
            return cls()

def main(bam_file, ref_path, tel_breaks, assembly_str, mapq_min, size, resolution, output_path, output_type):
    assembly = get_assembly(assembly_str)
    assembly.initialize(assembly_str, ref_path)
    
    # Actual telomere breaks
    breaks = []
    header = True
    with open(tel_breaks) as breaks_file:
        for breakpoint in breaks_file:
            if header:
                header = False
                continue
            chr, bp, dir = breakpoint.split()[:3]
            breaks.append((chr, int(bp), dir))
    
    print(f"Initial breaks: {len(breaks)}")
    valid_breaks = []
    for b in breaks:
        if assembly.valid_locus(b, size, resolution, bam_file):
            valid_breaks.append(b)
    print(f"Valid breaks: {len(valid_breaks)}")
    
    breaks_coverage, breaks_stderr = coverage_near_multiple_loci(bam_file, valid_breaks, size, resolution, mapq_min)
    
    # Control breaks
    control_breaks = assembly.n_random_loci(len(valid_breaks), size, resolution, bam_file)
    print(f"Number of control breaks (should be the same as valid breaks): {len(control_breaks)}")
    
    controls_coverage, controls_stderr = coverage_near_multiple_loci(bam_file, control_breaks, size, resolution, mapq_min)
    
    if output_type == 'data':
        output_data(breaks_coverage, breaks_stderr, controls_coverage, controls_stderr, size, resolution, output_path)
    elif output_type == 'plot':
        plot_and_save(breaks_coverage, controls_coverage, output_path, size, resolution)

main(args.bam_file, args.ref_path, args.tel_breaks, args.assembly, args.mapq_min, args.size, args.resolution, args.output_path, args.output_type)
