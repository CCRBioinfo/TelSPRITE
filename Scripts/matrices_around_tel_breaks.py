import argparse
import random
import numpy as np
from collections import OrderedDict
import assembly_methods_mod as am

__author__ = "David Wilson"

parser = argparse.ArgumentParser()
parser.add_argument(
    '-m', 
    '--in_matrix', 
    required=True
    )
parser.add_argument(
    '-b', 
    '--breaks', 
    required=True
    )
parser.add_argument(
    '-p', 
    '--bam_path', 
    required=True
    )
parser.add_argument(
    '-a', 
    '--assembly', 
    required=True
    )
parser.add_argument(
    '-f', 
    '--ref_path', 
    required=True
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
    '-o', 
    '--output',
    required=True
    )
args = parser.parse_args()

class Assembly:
    
    def __init__(self, size, resolution, assembly, ref_path, exclude_Y=True, exclude_M=True):
        if exclude_Y:
            self.chr_sizes = OrderedDict([(chr, size) for (chr, size) in self.chr_sizes.items() if chr not in ['chrY', 'CP086569.2']])
        if exclude_M:
            self.chr_sizes = OrderedDict([(chr, size) for (chr, size) in self.chr_sizes.items() if chr not in ['chrM']])
        self.chr_names = [chr for (chr, size) in self.chr_sizes.items()]
        self.chr_weights = [size for (chr, size) in self.chr_sizes.items()]
        self.size = size
        self.resolution = resolution
        
        self.am_methods = am.get_assembly(assembly)
        self.am_methods.load_ref_seq(ref_path)
        
        # Parameters
        self.chr_end_prox_threshold = 100000
        self.locus_qual_threshold = 10
        self.locus_qual_frac = 0.4
    
#    def multiple_random_loci(self, n, include_dir=True):
#        # Returns a large number of random loci
#        # Will run much faster than individually calling .random_locus() many times
#        random_chrs = random.choices(self.chr_names, weights=self.chr_weights, k=n)
#        
#        random_loci = []
#        
#        for chr in self.chr_names:
#            n_chr = random_chrs.count(chr)
#            random_bps = random.choices(list(range(self.chr_sizes[chr])), k=n_chr)
#            if include_dir:
#                random_loci += [(chr, bp, random.choice(['right', 'left'])) for bp in random_bps]
#            else:
#                random_loci += [(chr, bp) for bp in random_bps]
#        
#        return random_loci
    
#    def random_locus(self, include_dir=True):
#        # Returns a randomly selected genomic locus
#        random_chr = random.choices(self.chr_names, weights=self.chr_weights, k=1)[0]
#        random_locus = random.choice(range(self.chr_sizes[random_chr]))
#        
#        if include_dir:
#            return (random_chr, random_locus, random.choice(['right', 'left']))
#        else:
#            return (random_chr, random_locus)
    
    def load_matrix(self, matrix_path):
        matrix = []
        for line in open(matrix_path):
            matrix.append([float(x) for x in line.split()])
        matrix = np.array(matrix)
        
        return matrix
    
    def valid_locus(self, loc, size, resolution, bam_path):
        # Tests if locus is valid
        chr, bp = loc[:2]
        bin = bp // resolution
        if bin - size < 0 or bin + size > self.chr_sizes[chr]//resolution:
            return False
        else:
            if self.am_methods.valid_locus((chr, bp), self.chr_end_prox_threshold, bam_path, self.locus_qual_threshold, self.locus_qual_frac):
                return True
            else:
                return False
    
    def contacts_around_pos(self, pos, matrix):
        chr, bp, dir = pos
        bin = bp // self.resolution
        start = bin - self.size
        end = bin + self.size
        
        out_matrix = matrix[start:end+1, start:end+1]
        if dir == "left":
            out_matrix = out_matrix[::-1, ::-1]
        return out_matrix
    
    def load_breaks(self, breaks_path):
        breaks = []
        header = True
        with open(breaks_path) as breaks_file:
            for breakpoint in breaks_file:
                if header:
                    header = False
                    continue
                chr, bp, dir = breakpoint.split()[:3]
                #bin = int(bp)//self.resolution
                breaks.append((chr, int(bp), dir))
        
        return breaks
    
    def multiple_pos(self, positions, chr1_path):
        first = True
        
        for pos in positions:
            matrix = self.load_matrix(chr1_path.replace('chr1', pos[0]))
            if first:
                out_matrix = self.contacts_around_pos(pos, matrix)
                #print("First matrix...")
                first = False
            else:
                out_matrix = np.add(out_matrix, self.contacts_around_pos(pos, matrix))
        
        return out_matrix
    
    def write_out(self, matrix, out_path):
        with open(out_path, "w") as out_file:
            out_str = ''
            for x in matrix:
                next_line = []
                for y in x:
                    next_line.append(y)
                out_str += '\t'.join([str(y) for y in next_line]) + '\n'
            out_file.write(out_str)
    
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
    
    def main(self, in_matrix, breaks_path, out_path, size, resolution, bam_file):
        breaks = self.load_breaks(breaks_path)
        print(f"Initial breaks: {len(breaks)}")
        valid_breaks = []
        for b in breaks:
            if self.valid_locus(b, size, resolution, bam_file):
                valid_breaks.append(b)
        print(f"Valid breaks: {len(valid_breaks)}")
        breaks = self.multiple_pos(valid_breaks, in_matrix)
        breaks_out = out_path + '_breaks.txt'
        self.write_out(breaks, breaks_out)
        
        n = len(valid_breaks)
        control_breaks = self.n_random_loci(n, size, resolution, bam_file)
        print(f"Number of control breaks (should be the same as valid breaks): {len(control_breaks)}")
        control = self.multiple_pos(control_breaks, in_matrix)
        control_out = out_path + '_control.txt'
        self.write_out(control, control_out)

class Hg19(Assembly):
    
    assembly_name = 'hg19'

    def __init__(self, size, resolution, assembly, ref_path):
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
        super(Hg19, self).__init__(size, resolution, assembly, ref_path)

class Hg38(Assembly):
    
    assembly_name = 'hg38'

    def __init__(self, size, resolution, assembly, ref_path):
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
        super(Hg38, self).__init__(size, resolution, assembly, ref_path)

class T2T_CHM13v2_0(Assembly):
    
    assembly_name = 't2t_chm13'

    def __init__(self, size, resolution, assembly, ref_path):
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
        super(T2T_CHM13v2_0, self).__init__(size, resolution, assembly, ref_path)

def get_assembly(assembly_name):
    for cls in Assembly.__subclasses__():
        if cls.assembly_name == assembly_name:
            return cls()

hg38 = Hg38(args.size, args.resolution, args.assembly, args.ref_path)
hg38.main(args.in_matrix, args.breaks, args.output, args.size, args.resolution, args.bam_path)
