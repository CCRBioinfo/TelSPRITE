from collections import OrderedDict
import random
import sys
import pysam

__author__ = "David Wilson"

class Assembly:
    
    def __init__(self, exclude_Y=True, exclude_M=True):
        if exclude_Y:
            self.chr_sizes = OrderedDict([(chr, size) for (chr, size) in self.chr_sizes.items() if chr not in ['chrY', 'CP086569.2']])
        if exclude_M:
            self.chr_sizes = OrderedDict([(chr, size) for (chr, size) in self.chr_sizes.items() if chr not in ['chrM']])
        self.chr_names = [chr for (chr, size) in self.chr_sizes.items()]
        self.chr_weights = [size for (chr, size) in self.chr_sizes.items()]
    
    def random_locus(self):
        # Returns a randomly selected genomic locus
        random_chr = random.choices(self.chr_names, weights=self.chr_weights, k=1)[0]
        random_locus = random.choice(range(self.chr_sizes[random_chr]))
        
        return (random_chr, random_locus)
    
    def multiple_random_loci(self, n):
        # Returns a large number of random loci
        # Will run much faster than individually calling .random_locus() many times
        random_chrs = random.choices(self.chr_names, weights=self.chr_weights, k=n)
        
        random_loci = []
        
        for chr in self.chr_names:
            n_chr = random_chrs.count(chr)
            random_bps = random.choices(list(range(self.chr_sizes[chr])), k=n_chr)
            random_loci += [(chr, bp) for bp in random_bps]
        
        return random_loci
    
    def load_ref_seq(self, ref_seq_path):
        # Loads the entire sequence of the reference genome, provided path to .fa file
        NAMES_CONVERT = {
            'chr1': 'NC_000001.11',
            'chr2': 'NC_000002.12',
            'chr3': 'NC_000003.12',
            'chr4': 'NC_000004.12',
            'chr5': 'NC_000005.10',
            'chr6': 'NC_000006.12',
            'chr7': 'NC_000007.14',
            'chr8': 'NC_000008.11',
            'chr9': 'NC_000009.12',
            'chr10': 'NC_000010.11',
            'chr11': 'NC_000011.10',
            'chr12': 'NC_000012.12',
            'chr13': 'NC_000013.11',
            'chr14': 'NC_000014.9',
            'chr15': 'NC_000015.10',
            'chr16': 'NC_000016.10',
            'chr17': 'NC_000017.11',
            'chr18': 'NC_000018.10',
            'chr19': 'NC_000019.10',
            'chr20': 'NC_000020.11',
            'chr21': 'NC_000021.9',
            'chr22': 'NC_000022.11',
            'chrX': 'NC_000023.11'
        }
        NAMES_CONVERT = {value: key for key, value in NAMES_CONVERT.items()}
        
        NAMES_CONVERT_GB = {
            'chr1': 'CM000663.2',
            'chr2': 'CM000664.2',
            'chr3': 'CM000665.2',
            'chr4': 'CM000666.2',
            'chr5': 'CM000667.2',
            'chr6': 'CM000668.2',
            'chr7': 'CM000669.2',
            'chr8': 'CM000670.2',
            'chr9': 'CM000671.2',
            'chr10': 'CM000672.2',
            'chr11': 'CM000673.2',
            'chr12': 'CM000674.2',
            'chr13': 'CM000675.2',
            'chr14': 'CM000676.2',
            'chr15': 'CM000677.2',
            'chr16': 'CM000678.2',
            'chr17': 'CM000679.2',
            'chr18': 'CM000680.2',
            'chr19': 'CM000681.2',
            'chr20': 'CM000682.2',
            'chr21': 'CM000683.2',
            'chr22': 'CM000684.2',
            'chrX': 'CM000685.2'
        }
        NAMES_CONVERT_GB = {value: key for key, value in NAMES_CONVERT_GB.items()}
        
        self.ref_seq = OrderedDict([(chr, '') for chr in self.chr_names])
        
        ref_seq_process = open(ref_seq_path).read()
        ref_seq_process = [x for x in ref_seq_process.split('>') if len(x) > 0]
        for chr_seq in ref_seq_process:
            current_chr = False
            first_line = chr_seq[:chr_seq.index('\n')]
            i = 0
            while i < len(first_line):
                if first_line[:i+1] in self.chr_names:
                    current_chr = first_line[:i+1]
                elif first_line[:i+1] in NAMES_CONVERT:
                    current_chr = NAMES_CONVERT[first_line[:i+1]]
                elif first_line[:i+1] in NAMES_CONVERT_GB:
                    current_chr = NAMES_CONVERT_GB[first_line[:i+1]]
                i += 1
            
            if current_chr and len(self.ref_seq[current_chr]) == 0:
                #print("Loading "+current_chr+"...")
                self.ref_seq[current_chr] = chr_seq[chr_seq.index('\n'):] # Remove first line with chr name
                self.ref_seq[current_chr] = self.ref_seq[current_chr].upper() # Ensure all bases are same case
                #self.ref_seq[current_chr] = ''.join([x for x in self.ref_seq[current_chr] if x in ['A', 'T', 'C', 'G', 'N']]) # Remove all new line characters
                self.ref_seq[current_chr] = ''.join([x for x in self.ref_seq[current_chr] if x != '\n'])
        
        # Check that chromosome sizes match expectations
        if OrderedDict([(chr, len(seq)) for (chr, seq) in self.ref_seq.items()]) != self.chr_sizes:
            for chr, seq in self.ref_seq.items():
                print(chr, len(seq))
            sys.exit("Improper reference sequence loaded. Chromosome sizes are not correct for this assembly or headers are incorrect.")
        
        self.get_chr_ends()
    
    def get_chr_ends(self):
        # Determine the first and last non-"N" base pair of each chromosome
        # Intended to be called after loading reference sequence
        self.chr_ends = OrderedDict([(chr, '') for chr in self.chr_names])
        
        for chr in self.ref_seq:
            chr_str = self.ref_seq[chr]
            chr_len = len(chr_str)
            
            i_for = 0
            while chr_str[i_for] == 'N':
                i_for += 1
            
            i_rev = len(chr_str)-1
            while chr_str[i_rev] == 'N':
                i_rev -= 1
            
            self.chr_ends[chr] = (i_for, i_rev)
    
    def get_sequence_at_locus(self, locus, length, direction):
        # Return the DNA sequence of the reference just before or after a locus
        # Must load reference sequence first
        if not hasattr(self, 'ref_seq'):
            sys.exit('Error: cannot run .get_sequence_at_locus() without loading reference genome sequence first')
        
        chr, bp = locus
        
        if direction == "left":
            return self.ref_seq[chr][bp-length+1:bp+1]
        elif direction == "right":
            return self.ref_seq[chr][bp:bp+length]

    def valid_locus(self, locus, chr_end_prox_threshold, bam_path, locus_qual_threshold, locus_qual_frac):
        # Check if a genomic locus is "valid" for purposes of counting it as a telo insertion site
        # Loci not on a chromosome listed in the assembly are invalid
        # Loci very close to the end of the chromosome are invalid
        # Loci in regions that are very poorly mapped are invalid
        # Loci represented as "N" nucleotides are invalid
        # Must load reference sequence first
        if not hasattr(self, 'ref_seq'):
            sys.exit('Error: cannot run .valid_locus() without loading reference genome sequence first')

        # Test if chromosome is in assembly
        chr, bp = locus
        if chr not in self.chr_names:
            return False
        
        # Test if locus is too close to chromosome end
        chr_start, chr_end = self.chr_ends[chr]
        if (bp <= chr_start+chr_end_prox_threshold) or (bp >= chr_end-chr_end_prox_threshold):
            return False
        
        # Test if locus is poorly mapped
        overlapping_read_mapqs = []
        
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for read in bam_file.fetch(chr, bp-1, bp+1):
                overlapping_read_mapqs.append(read.mapping_quality)

        if len(overlapping_read_mapqs) == 0:
            return False
        
        below_threshold = [mapq for mapq in overlapping_read_mapqs if mapq <= locus_qual_threshold]
        if len(below_threshold)/len(overlapping_read_mapqs) >= locus_qual_frac:
            return False
        
        # Test if locus is an "N" nucleotide
        if self.ref_seq[chr][bp] == "N":
            return False
        
        # Get here if none of the checks failed
        return True

    def calc_offset(self, target_chr, resolution):
        offset = 0
        
        for chr in self.chr_names[:self.chr_names.index(target_chr)]:
            offset += -(-self.chr_sizes[chr]//resolution)
        
        return offset

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

