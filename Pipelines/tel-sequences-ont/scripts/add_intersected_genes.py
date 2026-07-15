import argparse
from tel_table import Tel_table
from gene_intersect import Gene_intersect

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
    '-g', 
    '--gene_annotation',
    required=True
    )
args = parser.parse_args()

tel_table = Tel_table()
tel_table.load_tel_table(args.input)

gene_intersect = Gene_intersect(args.gene_annotation)

for tel in tel_table.tels:
    locus = (tel['Chromosome'], tel['Start position'])
    gene_intersected = gene_intersect.locus_gene_intersect(locus)
    tel['Genes intersected'] = gene_intersected['name'] if gene_intersected else ' '

tel_table.col_names.append('Genes intersected')
tel_table.write_tel_table(args.output)

