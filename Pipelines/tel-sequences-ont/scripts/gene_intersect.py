__author__ = "David Wilson"

class Gene_intersect:
    
    def __init__(self, gene_annotation_refseq_path):
        self.get_gene_locations(gene_annotation_refseq_path)
    
    def get_gene_locations(self, gene_annotation_refseq_path):
        self.genes = {}
        first_row = True
        
        with open(gene_annotation_refseq_path) as gene_annotation_file:
            for line in gene_annotation_file:
                if first_row:
                    # Ignore first row with header
                    first_row = False
                else:
                    bin, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, score, name2, cdsStartStat, cdsEndStat, exonFrames = line.split()
                    if chrom not in self.genes:
                        self.genes[chrom] = []
                    self.genes[chrom].append({'txStart': txStart, 'txEnd': txEnd, 'name2': name2})

    def gene_size(self, name):
        # Returns the size of the gene with a given name
        for gene in self.genes:
            if gene['name2'] == 'name':
                return abs(gene['txEnd'] - gene['txStart'])
    
    def locus_gene_intersect(self, locus):
        chr, start = locus
        gene_intersected = False
        
        for gene in self.genes[chr]:
            gene_start = int(gene['txStart'])
            gene_end = int(gene['txEnd'])
            if gene_start <= int(start) <= gene_end:
                gene_intersected = True
                name = gene['name2']
                size = gene_end-gene_start
        
        if gene_intersected:
            return {'name': name, 'size': size}
        else:
            return False

