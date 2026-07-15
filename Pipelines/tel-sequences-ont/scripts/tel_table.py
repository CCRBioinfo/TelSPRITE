__author__ = "David Wilson"

class Tel_table:

    def __init__(self):
        self.tels = []

    def add_tel(self, tel):
        # Adds new telomere with list as input
        # Assumes that list elements are in the same order as column names from self.col_names
        self.tels.append({self.col_names[i]:tel[i] for i in range(len(tel))})
    
    def load_tel_table(self, in_table_path):
        # Load input telomere table and store telomere rows as well as column names
        first_line = True
        
        with open(in_table_path) as in_file:
            for line in in_file:
                if first_line:
                    first_line = False
                    self.col_names = [name.replace('\n', '') for name in line.split('\t')]
                else:
                    self.add_tel([name.replace('\n', '') for name in line.split('\t')])
    
    def write_tel_table(self, out_table_path):
        # Write telomere table as output
        with open(out_table_path, "w") as out_file:
            out_file.write('\t'.join(self.col_names))
            
            for tel in self.tels:
                new_line = [tel[col_name] for col_name in self.col_names]
                out_file.write('\n'+'\t'.join(new_line))
    
    def tel_pos(self, min_read_count=1):
        # Returns positions of telomere insertions from table for downstream analysis
        return [(tel['Chromosome'], int(tel['Start position'])) for tel in self.tels if int(tel['Read count']) >= min_read_count]

