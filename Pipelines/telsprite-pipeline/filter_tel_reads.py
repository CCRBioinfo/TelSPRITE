#!/usr/bin/env python

import sys
import os
import gzip

__author__ = "Sarah Clatterbuck Soper, David Wilson"

threshold_total = 4
threshold_canonical = 2
n = 4


def process(lines=None):
    ks = ["name", "sequence", "optional", "quality"]
    return {k: v for k, v in zip(ks, lines)}


def formatRecord(myRecord):
    nr = (
        myRecord["name"]
        + "\n"
        + myRecord["sequence"]
        + "\n"
        + myRecord["optional"]
        + "\n"
        + myRecord["quality"]
        + "\n"
    )
    return nr

def reverse_complement(sequence):
    base_pairing = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    complement = [base_pairing[base] for base in sequence]
    complement.reverse()
    return ''.join(complement)
    
def non_overlapping_count(read, CTR, VTR_list):
    n_CTR = read.count(CTR)
    read = read.replace(CTR, '-')
    
    n_VTR = 0
    for VTR in VTR_list:
        n_VTR += read.count(VTR)
        read = read.replace(VTR, '-')
    
    return (n_CTR, n_VTR)

try:
    fn = sys.argv[1]
except IndexError as ie:
    raise SystemError("Error: Specify file name\n")

if not os.path.exists(fn):
    raise SystemError("Error: File does not exist\n")

try:
    outdir = sys.argv[2]
except IndexError as ie:
    raise SystemError("Error: Specify output directory\n")

if not os.path.exists(outdir):
    raise SystemError("Error: Output directory does not exist\n")

try:
    threshold_canonical = int(sys.argv[3])
    threshold_variant = int(sys.argv[4])
    threshold_total = int(sys.argv[5])
    threshold_frac_canonical = int(sys.argv[6])
except IndexError as ie:
    raise SystemError("Error: Did not include all necessary telomere filtering parameters\n")

if os.path.basename(fn).split(".")[-1] == "gz":
    print("Opening gzipped file ", fn)
    bfn = ".".join(os.path.basename(fn).split(".")[:-2])
    fh = gzip.open(fn, "rt")
else:
    print("Opening plain text file", fn)
    bfn = ".".join(os.path.basename(fn).split(".")[:-1])
    fh = open(fn, "r")

noTelo = outdir + "/" + bfn + ".noTelo.fastq.gz"
Telo = outdir + "/" + bfn + ".Telo.fastq.gz"

try:
    noTeloFH = gzip.open(noTelo, "wt")
except OSError:
    print("Could not open/read file: ", noTelo)
    sys.exit()

try:
    TeloFH = gzip.open(Telo, "wt")
except OSError:
    print("Could not open/read file: ", Telo)
    sys.exit()


repeats = [
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
telSeqs = []
telSeqs_rc = []

for repeat in repeats:
    telSeqs.append(repeat)
    telSeqs_rc.append(reverse_complement(repeat))

lines = []

for line in fh:
    line = str(line)
    lines.append(line.rstrip())
    if len(lines) == n:
        record = process(lines)
        mySeq = record["sequence"]
        canonical, variant = non_overlapping_count(mySeq, "TTAGGG", telSeqs)
        canonical_rc, variant_rc = non_overlapping_count(mySeq, "CCCTAA", telSeqs_rc)
        if (canonical + variant) < (canonical_rc + variant_rc):
            canonical = canonical_rc
            variant = variant_rc
        total = canonical + variant
        record["name"] += " CTR="
        record["name"] += str(canonical)
        record["name"] += " VTR="
        record["name"] += str(variant)
        newRecord = formatRecord(record)
        if (
            (canonical >= threshold_canonical) and
            (variant >= threshold_variant) and
            (total >= threshold_total) and
            (canonical / total >= threshold_frac_canonical)
            ):
            TeloFH.write(newRecord)
        else:
            noTeloFH.write(newRecord)
        lines = []

