import argparse
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
    '-m', 
    '--min_read_count',
    required=True
    )
args = parser.parse_args()

tel_table = Tel_table()
tel_table.load_tel_table(args.input)

tel_table.tels = [tel for tel in tel_table.tels if tel['Read count'] >= args.min_read_count]

tel_table.write_tel_table(args.output)
