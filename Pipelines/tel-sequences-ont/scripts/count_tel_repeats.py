__author__ = "David Wilson"

g_repeats = [
    "TTAGGG",
    "TCAGGG",
    "TGAGGG",
    "TTGGGG",
    "TTCGGG",
    "TTTGGG",
    "ATAGGG",
    "CATGGG",
    "CTAGGG",
    "GTAGGG",
    "TAAGGG",
]

c_repeats = [
    "CCCTAA", 
    "CCCTGA", 
    "CCCTCA", 
    "CCCCAA", 
    "CCCGAA", 
    "CCCAAA", 
    "CCCTAT", 
    "CCCATG", 
    "CCCTAG", 
    "CCCTAC", 
    "CCCTTA"
]

def shared_pattern(repeats):
    # Returns maximum shared pattern between variant tel repeats
    # Searching for this first can be used to increase the speed of the tel filter
    longest = ''
    r_len = len(repeats[0])
    repeats_ext = [repeat*2 for repeat in repeats]
    
    start = 0
    while start < r_len:
        for x in range(r_len):
            end = start+x+1
            patterns = set()
            for repeat in repeats_ext:
                patterns.add(repeat[start:end])
            if len(patterns) == 1:
                if len(list(patterns)[0]) > len(longest):
                    longest = list(patterns)[0]
        start += 1
    
    return longest

g_shared = shared_pattern(g_repeats)
c_shared = shared_pattern(c_repeats)

def enough_repeats(seq, repeats, repeat_threshold, shared=False):
    # Determines if sequence meets telomere repeat threshold
    # Replacing counted repeats with a dummy character ('-') prevents overlap counting
    # Ex. 'CCCCAAA' is counted only once
    # Checking for shared pattern first can increase speed in situations where many sequences are expected to fail
    if shared:
        if not seq.count(shared) >= repeat_threshold:
            return False
    
    n_repeats = 0
    for repeat in repeats:
        n_repeats += seq.count(repeat)
        if n_repeats >= repeat_threshold:
            return True
        seq = seq.replace(repeat, '-')
    
    return False

def seq_is_tel(seq, g_repeats, c_repeats, g_shared, c_shared, repeat_threshold):
    # Determines if the sequence has enough G or C repeats
    if enough_repeats(seq, g_repeats, repeat_threshold, shared=g_shared):
        return 'TTAGGG'
    elif enough_repeats(seq, c_repeats, repeat_threshold, shared=c_shared):
        return 'CCCTAA'
    else:
        return False

