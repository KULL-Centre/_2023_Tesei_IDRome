#!/usr/bin/bash

usage="usage: ./run_mafft <input.fasta> [out_dir]"

echo "cmd:$0 $@"

[[ $1 ]] || { echo "ERROR: No input"; echo $usage; exit 2; }
[[ -f $1 ]] || { echo "ERROR: input '$1' does not excist"; exit 2; }
[[ $1 == *.fasta ]] || { echo "ERROR: input '$1' does not have .fasta extension"; exit 2; }

if [ $2 ]; then
    [[ -d $2 ]] || { echo "ERROR: Output dir '$2' does not excist (second argument)"; exit 2; }
    outdir=$2
else
    outdir="."
fi

infile=$(basename $1)
fileid=${infile%.fasta}
outfile=${outdir}/${fileid}_aln.fasta

mafft --maxiterate 1000 --quiet --thread 0 $1 1> $outfile 2> ${outdir}/${fileid}.out

[[ $? == 0 ]] || { echo "WARNING: Bad exit statur for '$1'"; exit 1; }
exit 0
