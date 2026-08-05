import argparse
import sys
from collections import OrderedDict
from statistics import mean
from scipy.stats import sem

parser = argparse.ArgumentParser()
parser.add_argument('input')
parser.add_argument('assembly', choices=['hg19', 'hg38', 'chm13', 'mm9', 'mm10'])
parser.add_argument('resolution')
parser.add_argument('output')
parser.add_argument('-e', '--exclude_chr', action='store_false')
parser.add_argument('-z', '--ignore_zeros', action='store_false')
parser.add_argument(
    '-t',
    '--output_type',
    default='all',
    choices=['all', 'telo_v_dist', 'telo_contacts', 'normalized_telo_contacts'])
parser.add_argument(
    '-s',
    '--output_style',
    default='all',
    choices=['all', 'table', 'bedgraph'])
parser.add_argument('-m', '--min_data_points', default=1)
args = parser.parse_args()

class Telo_bin:
    
    def __init__(self, chr, bin, contact_prob, assembly_chr_bins):
        self.chr = chr
        self.bin = bin
        self.contact_prob = contact_prob
        self.assembly_chr_bins = assembly_chr_bins
        self.distance = self.get_distance()
    
    def get_distance(self):
        # Returns distance to nearest telomere
        chr_bin = self.bin + 1
        chr_num_bins = self.assembly_chr_bins[self.chr]
        
        dist_telo1 = chr_bin # First bin has a distance of 1
        dist_telo2 = chr_num_bins - chr_bin + 1 # Last bin also has a distance of 1
        
        return min([dist_telo1, dist_telo2])
    
    def normalized_contact_prob(self, telo_vs_dist_dict):
        # Use external average telo contact probability vs distance dictionary to normalize contact value
        try:
            avg, std_err = telo_vs_dist_dict[self.distance]
            norm_contact = self.contact_prob / avg
        except:
            # Errors if the given distance is not in the dictionary
            norm_contact = 0
        
        return norm_contact

class Assembly:
    
    def __init__(self, resolution):
        self.res = resolution
        self.chr_names = list(self.chr_sizes.keys())
        self.chr_bins = OrderedDict([(chr, -(-size//self.res)) for (chr, size) in self.chr_sizes.items()])
        
        if self.chr_bins['chrT'] != 1:
            sys.exit("Make sure that the number of telomere bins equals 1")
        else:
            # Determine index of matrix position corresponding to chrT
            chrs = []
            for chr in self.chr_bins:
                chrs += [chr]*self.chr_bins[chr]
            self.telo_index = chrs.index('chrT')
    
    def main(self, input_path, n_min, output_path, output_type, output_style, exclude_chr=True, ignore_zeros=True):
        # Parent function that calls all functions in proper order
        self.load_matrix(input_path, exclude_chr)
        self.sort_by_distance()
        self.stats_for_all_telo(n_min, ignore_zeros)
        
        for type in output_type:
            if type == 'telo_v_dist':
                self.output_telo_by_distance(output_path)
            else:
                for style in output_style:
                    if type == 'telo_contacts':
                        self.output_contact_data(output_path, style, norm=False)
                    elif type == 'normalized_telo_contacts':
                        self.output_contact_data(output_path, style, norm=True)
    
    def load_matrix(self, input_path, exclude_chr):
        # Loads matrix file and generates a dictionary with telo contacts
        self.telo_contacts = OrderedDict([(chr, []) for chr in self.chr_names])
        
        chr = self.chr_names[0]
        bin = -1
        
        with open(input_path) as input_file:
            for row in input_file:
                bin += 1
                if bin == self.chr_bins[chr]:
                    chr = self.chr_names[self.chr_names.index(chr)+1]
                    bin = 0
                
                self.telo_contacts[chr].append(Telo_bin(chr, bin, float(row.split()[self.telo_index]), self.chr_bins))
        
        if exclude_chr:
            # Recommended to exclude chromosomes Y, M, and T
            self.telo_contacts = OrderedDict([(chr, telos) for (chr, telos) in self.telo_contacts.items() if chr not in ['chrY', 'chrM', 'chrT']])
    
    def sort_by_distance(self):
        # Sorts telo contacts by distance to nearest telomere
        # Must load a matrix file first
        self.telo_by_distance = {}
        
        for chr in self.telo_contacts:
            for telo_bin in self.telo_contacts[chr]:
                dist = telo_bin.distance
                if dist not in self.telo_by_distance:
                    self.telo_by_distance[dist] = []
                self.telo_by_distance[dist].append(telo_bin)
        
        sort_keys = list(self.telo_by_distance.keys())
        sort_keys.sort()
        self.telo_by_distance = OrderedDict([(dist, self.telo_by_distance[dist]) for dist in sort_keys])
    
    def stats_for_contacts_at_distance(self, dist, ignore_zeros):
        # Calculates average and standard error for telo contacts at a given distance
        # Must sort telo contacts by distance first
        contacts = [x.contact_prob for x in self.telo_by_distance[dist]]
        
        # Recommended to exclude 0's, which generally represent bins lacking mappable DNA
        if ignore_zeros:
            contacts = [x for x in contacts if x != 0]
        
        n = len(contacts)
        avg = mean(contacts)
        std_err = sem(contacts)
        
        return (n, avg, std_err)
    
    def stats_for_all_telo(self, n_min, ignore_zeros):
        # Calculates telo average and standard error at all distances
        # Must sort telo contacts by distance first
        self.stats_by_distance = OrderedDict([])
        
        for dist in self.telo_by_distance:
            try:
                n, avg, std_err = self.stats_for_contacts_at_distance(dist, ignore_zeros)
            except: # Errors if there are no nonzero data points at this distance
                continue
            if n >= n_min:
                self.stats_by_distance[dist] = (avg, std_err)
    
    def output_telo_by_distance(self, output_path):
        # Outputs telo statistics in multiple formats
        # Must calculate telo vs distance statistics first
        
        with open(output_path + '_telo_by_distance.txt', "w") as out_table:
            for dist in self.stats_by_distance:
                avg, std_err = self.stats_by_distance[dist]
                new_line = '\t'.join([str(x) for x in [dist, avg, std_err]])
                out_table.write(new_line + '\n')
    
    def output_contact_data(self, output_path, style, norm=False):
        # Outputs telo contacts (raw or normalized by distance) in each bin
        # Must calculate telo vs distance statistics first to output normalized contacts
        if norm:
            output_path += '_normalized_telo_contacts'
        else:
            output_path += '_telo_contacts'
        
        if style == 'table':
            output_path += '.txt'
        elif style == 'bedgraph':
            output_path += '.bedgraph'
        
        with open(output_path, "w") as out_telo:
            if style == 'bedgraph':
                out_telo.write('track type=bedGraph name="BedGraph Format" description="BedGraph format" visibility=full color=200,100,0 altColor=0,100,200 priority=20')
           
            for chr in self.telo_contacts:
                for telo_bin in self.telo_contacts[chr]:
                    bin = telo_bin.bin
                    
                    if norm:
                        contact_prob = telo_bin.normalized_contact_prob(self.stats_by_distance)
                    else:
                        contact_prob = telo_bin.contact_prob
                    
                    if style == 'table':
                        new_line = '\t'.join([str(x) for x in [chr, bin, contact_prob]])
                    elif style == 'bedgraph':
                        new_line = '\t'.join([str(x) for x in [chr, (bin*self.res), ((bin+1)*self.res-1), contact_prob]])
                    
                    out_telo.write(new_line + '\n')

class Hg19(Assembly):

    def __init__(self, resolution):
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
            ('chrM', 16571),
            ('chrT', 1000)])
        super(Hg19, self).__init__(resolution)

class Hg38(Assembly):

    def __init__(self, resolution):
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
            ('chrM', 16569),
            ('chrT', 1000)])
        super(Hg38, self).__init__(resolution)

class Chm13(Assembly):

    def __init__(self, resolution):
        self.chr_sizes = OrderedDict([
            ('chr1', 248387328),
            ('chr2', 242696752),
            ('chr3', 201105948),
            ('chr4', 193574945),
            ('chr5', 182045439),
            ('chr6', 172126628),
            ('chr7', 160567428),
            ('chr8', 146259331),
            ('chr9', 150617247),
            ('chr10', 134758134),
            ('chr11', 135127769),
            ('chr12', 133324548),
            ('chr13', 113566686),
            ('chr14', 101161492),
            ('chr15', 99753195),
            ('chr16', 96330374),
            ('chr17', 84276897),
            ('chr18', 80542538),
            ('chr19', 61707364),
            ('chr20', 66210255),
            ('chr21', 45090682),
            ('chr22', 51324926),
            ('chrX', 154259566),
            ('chrY', 62460029),
            ('chrM', 16569),
            ('chrT', 1000)])
        super(Chm13, self).__init__(resolution)

class Mm9(Assembly):

    def __init__(self, resolution):
        self.chr_sizes = OrderedDict([
            ('chr1', 197195432),
            ('chr2', 181748087),
            ('chr3', 159599783),
            ('chr4', 155630120),
            ('chr5', 152537259),
            ('chr6', 149517037),
            ('chr7', 152524553),
            ('chr8', 131738871),
            ('chr9', 124076172),
            ('chr10', 129993255),
            ('chr11', 121843856),
            ('chr12', 121257530),
            ('chr13', 120284312),
            ('chr14', 125194864),
            ('chr15', 103494974),
            ('chr16', 98319150),
            ('chr17', 95272651),
            ('chr18', 90772031),
            ('chr19', 61342430),
            ('chrX', 166650296),
            ('chrY', 15902555),
            ('chrM', 16299),
            ('chrT', 1000)])
        super(Mm9, self).__init__(resolution)

class Mm10(Assembly):

    def __init__(self, resolution):
        self.chr_sizes = OrderedDict([
            ('chr1', 195471971),
            ('chr2', 182113224),
            ('chr3', 160039680),
            ('chr4', 156508116),
            ('chr5', 151834684),
            ('chr6', 149736546),
            ('chr7', 145441459),
            ('chr8', 129401213),
            ('chr9', 124595110),
            ('chr10', 130694993),
            ('chr11', 122082543),
            ('chr12', 120129022),
            ('chr13', 120421639),
            ('chr14', 124902244),
            ('chr15', 104043685),
            ('chr16', 98207768),
            ('chr17', 94987271),
            ('chr18', 90702639),
            ('chr19', 61431566),
            ('chrX', 171031299),
            ('chrY', 91744698),
            ('chrM', 16299),
            ('chrT', 1000)])
        super(Mm10, self).__init__(resolution)

res = int(args.resolution)
min_data_points = int(args.min_data_points)

if args.output_type == 'all':
    output_type = ['telo_v_dist', 'telo_contacts', 'normalized_telo_contacts']
else:
    output_type = [args.output_type]

if args.output_style == 'all':
    output_style = ['table', 'bedgraph']
else:
    output_style = [args.output_style]

if args.assembly == 'hg19':
    assembly = Hg19(res)
elif args.assembly == 'hg38':
    assembly = Hg38(res)
elif args.assembly == 'chm13':
    assembly = Chm13(res)
elif args.assembly == 'mm9':
    assembly = Mm9(res)
elif args.assembly == 'mm10':
    assembly = Mm10(res)

assembly.main(
    args.input, 
    min_data_points, 
    args.output, 
    output_type,
    output_style,
    exclude_chr=args.exclude_chr,
    ignore_zeros=args.ignore_zeros)

