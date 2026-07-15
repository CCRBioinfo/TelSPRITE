import pysam
import argparse
import sys
from collections import OrderedDict
import count_tel_repeats as ct

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
    '-w', 
    '--window_size',
    required=True,
    type=int,
    default='300'
    )
parser.add_argument(
    '-t',
    '--threshold',
    required=True,
    type=float,
    default='0.9'
    )
args = parser.parse_args()

repeat_threshold = int(args.window_size*args.threshold/6)

def process_alignment(read_seq, cigar_tuples):
    # Takes read sequence and cigar tuples as input
    # Returns tuples with operations and corresponding DNA sequence
    processed_tuples = []
    i = 0
    
    for cigar_tuple in cigar_tuples:
        operation, length = cigar_tuple
        if operation in [0,1,4,7,8]:
            processed_tuples.append((operation, read_seq[i:i+length]))
            i += length
    
    read_len = sum([len(x[1]) for x in processed_tuples])

    if read_len != len(read_seq):
        # This error should never trigger if this function is written correctly
        sys.exit("Output tuples are not the same length as input read")
    
    return processed_tuples

def get_tel_windows(seq, g_repeats, c_repeats, g_shared, c_shared, window_size, repeat_threshold):
    seq_len = len(seq)
    
    # Immediately disqualify sequences that do not have enough repeats
    if not ct.seq_is_tel(seq, ct.g_repeats, ct.c_repeats, ct.g_shared, ct.c_shared, repeat_threshold):
        return False
    
    i = 0
    tel_windows = []
    
    while seq_len-i >= window_size:
        window = seq[i:i+window_size]
        is_tel = ct.seq_is_tel(window, ct.g_repeats, ct.c_repeats, ct.g_shared, ct.c_shared, repeat_threshold)
        if is_tel == 'TTAGGG':
            i_start = i
            while ct.enough_repeats(window, ct.g_repeats, repeat_threshold, shared=ct.g_shared) and seq_len-i >= window_size:
                # Keep stepping window until it is no longer in telomere sequence
                i += 1
                window = seq[i:i+window_size]
            i += window_size-1
            i_end = i-1
            tel_windows.append([i_start, i_end, "TTAGGG"])
        elif is_tel == 'CCCTAA':
            i_start = i
            while ct.enough_repeats(window, ct.c_repeats, repeat_threshold, shared=ct.c_shared) and seq_len-i >= window_size:
                # Keep stepping window until it is no longer in telomere sequence
                i += 1
                window = seq[i:i+window_size]
            i += window_size-1
            i_end = i-1
            tel_windows.append([i_start, i_end, "CCCTAA"])
        else:
            i += 1
    
    if len(tel_windows) == 0:
        return False
    
    all_windows = []
    
    i_max = len(tel_windows)-1
    i = 0
    
    while i<i_max:
        prev_end = tel_windows[i][1]
        next_start = tel_windows[i+1][0]
        all_windows.append(tel_windows[i])
        if next_start-prev_end > 1:
            all_windows.append([prev_end+1, next_start-1, 'non_tel'])
        i += 1
    all_windows.append(tel_windows[-1])
    
    if all_windows[0][0] > 0:
        all_windows = [[0, tel_windows[0][0]-1, 'non_tel']] + all_windows
    
    if all_windows[-1][1] < seq_len-1:
        all_windows.append([all_windows[-1][1]+1, seq_len-1, 'non_tel'])
    
    return all_windows

def get_tel_info(read, windows, dir):
    windows = [window[:3] for window in windows]
    
    G_total = 0
    C_total = 0
    
    for window in windows:
        start, end, type = window
        length = end-start+1
        if type == 'TTAGGG':
            G_total += length
        elif type == 'CCCTAA':
            C_total += length
    
    if G_total >= C_total:
        strand = 'TTAGGG'
    else:
        strand = 'CCCTAA'
    
    if dir == 'right':
        start = read.reference_end
        if strand == 'TTAGGG':
            orientation = 'standard'
        elif strand == 'CCCTAA':
            orientation = 'inverted'
    elif dir == 'left':
        start = read.reference_start
        if strand == 'TTAGGG':
            orientation = 'inverted'
        elif strand == 'CCCTAA':
            orientation = 'standard'
    
    chr = read.reference_name
    
    return [chr, start, dir, orientation]

def add_seq_to_windows(seq, windows):
    for i, window in enumerate(windows):
        start = window[0]
        end = window[1]
        window_seq = seq[start:end+1]
        windows[i].append(window_seq)
    
    return windows

tel_list = []

with pysam.AlignmentFile(args.input, "rb") as in_file:
    for read in in_file:
        processed_tuples = process_alignment(read.query_sequence, read.cigartuples)
        
        if len(processed_tuples) == 1:
            # Must have both soft-clipped and aligned sequence
            continue
        first = processed_tuples[0]
        last = processed_tuples[-1]
        
        for i, segment in enumerate([first, last]):
            operation, seq = segment
            if operation == 4:
                # Only test soft-clipped segments at ends for telomere sequence
                tel_windows = get_tel_windows(seq, ct.g_repeats, ct.c_repeats, ct.g_shared, ct.c_shared, args.window_size, repeat_threshold)
                if tel_windows:
                    if i == 0:
                        dir = 'left'
                    elif i == 1:
                        dir = 'right'
                    with_seq = add_seq_to_windows(seq, tel_windows)
                    tel_list.append(get_tel_info(read, tel_windows, dir) + with_seq)

with open(args.output, "w") as out_tel_table:
    for tel in tel_list:
        tel_info = tel[:4]
        window_coords = []
        window_seqs = []
        for window in tel[4:]:
            window_coords.append(str(window[0])+'-'+str(window[1])+' '+window[2])
            window_seqs.append(window[3])
        out_tel_table.write('\t'.join([str(x) for x in tel_info]) + '\t')
        out_tel_table.write('; '.join(window_coords)+'\t')
        out_tel_table.write('\t'.join(window_seqs)+'\t \n')

