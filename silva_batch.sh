#!/bin/bash

for f in *.csv; do
  filename=$(basename $f .csv)
  silva_filter_fasta.py ${filename}.fasta $f 70
  silva_consensus.py ${filename}_id_filt.fasta
  fasta2png.py ${filename}_id_filt.fasta ${filename}_id_filt.fasta.png
done
