#!/usr/bin/env python

import sys
import os

# from Bio import AlignIO
# from Bio.Align import AlignInfo
from Bio import SeqIO

import csv
from tqdm import tqdm

# Usage: python silva_filter_fasta.py <fasta alignments file> <csv SILVA output>

# Settings
out_filename = os.path.splitext(sys.argv[1])[0] + "_id_filt.fasta"
filter_threshold = sys.argv[3]  # 100 = 100%, not 1 = 100%

print("FASTA file read")

reads_wanted = []

with open(sys.argv[2], newline='') as csvfile:
    silva_csv = csv.DictReader(csvfile)
    print("CSV file read")
    for read in tqdm(silva_csv):
        if float(read['identity']) > filter_threshold:
            reads_wanted.append(read['sequence_identifier'])

print(str(len(reads_wanted)) + " reads above threshold")
output_reads = []

for record in tqdm(SeqIO.parse(sys.argv[1], 'fasta')):
    if record.id in reads_wanted:
        output_reads.append(record)

SeqIO.write(output_reads, out_filename, "fasta")
