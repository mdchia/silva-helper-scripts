#!/usr/bin/env python

import sys
import os

from Bio import AlignIO
from Bio.Align import AlignInfo

# Usage: python silva_consensus.py <fasta alignments file>

# Settings
ambiguous_char = "X"
out_filename = os.path.splitext(sys.argv[1])[0] + "_consensus.fasta"
consensus_threshold = 0.6  # 0.5 = 50% of reads needed to call consensus
min_reads = 3  # minimum number of reads required to assign the presence of a base

# read in alignments
alignment = AlignIO.read(sys.argv[1], 'fasta')
print("File read")
summary_align = AlignInfo.SummaryInfo(alignment)

print("Generating consensus")

# Generate consensus
consensus = summary_align.gap_consensus(consensus_threshold, ambiguous=ambiguous_char, require_multiple=min_reads)
consensus_seq = str(consensus)

# get rid of bases that can't be called at confidence threshold
consensus_seq = consensus_seq.replace(ambiguous_char, "")

# get rid of gaps
consensus_seq = consensus_seq.replace("-", "")

# for SILVA - get back to DNA
consensus_seq = consensus_seq.replace("U", "T")

print("Consensus for " + sys.argv[1] + " is " + str(len(consensus_seq)) + " bases long")

print("Writing file")

# write out consensus
with open(out_filename, 'w') as fasta:
    # create a fasta contig entry
    fasta.write(">" + sys.argv[1] + "_consensus\n")
    fasta.write(consensus_seq)
