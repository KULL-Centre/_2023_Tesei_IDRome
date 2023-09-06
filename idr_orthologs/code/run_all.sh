#!/bin/bash

wget https://omabrowser.org/All/oma-uniprot.txt.gz

wget https://omabrowser.org/All/oma-protein-annotations.txt.gz

wget https://omabrowser.org/All/oma-seqs.fa.gz
gunzip oma-seqs.fa.gz

wget https://omabrowser.org/All/oma-pairs.txt.gz

mkdir human_homologs
python3 oma_homologs.py > oma_homologs.out

mkdir human_alignments
ls human_homologs/*.fasta | parallel -P 48 ./run_mafft.sh {} human_alignments

python3 idr_homologs.py > idr_homologs.out
