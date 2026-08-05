import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument(
    '-i', 
    '--input',
    required=True
    )
parser.add_argument(
    '-c', 
    '--chr',
    required=True
    )
parser.add_argument(
    '-s', 
    '--start',
    required=True,
    type=int
    )
parser.add_argument(
    '-e', 
    '--end',
    required=True,
    type=int
    )
parser.add_argument(
    '-a', 
    '--hash_size',
    required=True,
    type=int
    )
parser.add_argument(
    '-t', 
    '--size_threshold',
    required=True,
    type=int
    )
parser.add_argument(
    '-o', 
    '--output',
    required=True
    )
args = parser.parse_args()

chr = args.chr
start = args.start
end = args.end
hash_size = args.hash_size
size_threshold = args.size_threshold

length = end-start
print('Plot length from start to end point: '+str(length))
if length%hash_size != 0:
    sys.exit('Exiting because plot length is not a multiple of hash size')
n_hashes = int(length/hash_size)
print('Number of hashes to cover plot length: '+str(n_hashes))

print('Maximum cluster size: '+str(size_threshold))

total = 0
used_clusters = 0

with open(args.input) as f:
    with open(args.output, "w") as out_file:
        cluster_count = 0
        for cluster in f:
            if 'chrT' not in cluster:
                continue
            
            reads = cluster.split()[1:]
            size = len(reads)
            total += 1/size
            
            #try:
            if size > size_threshold:
                continue
            #except:
            #    assert 1 == 1
            
            used_clusters += 1/size
            row_list = [0 for x in range(n_hashes)]
            
            for read in reads:
                _, coord = read.split('_',1)
                r_chr, r_start, r_end = coord.replace('-', ':').split(':')
                if r_chr == 'chrT':
                    continue
                r_start = int(r_start)
                r_end = int(r_end)
                if r_chr==chr and r_start>=start and r_end<=end:
                    row_list[(r_start-start)//hash_size] = 1
            
            if 1 in row_list:
                row_str = '\t'.join([str(x) for x in row_list])
                out_file.write(row_str + '\n')
            
                cluster_count += 1
        print("Clusters written: "+str(cluster_count))
        print("Clusters account for "+str(used_clusters/total*100)+"% of the contact probability")
