#!/bin/sh
#SBATCH --job-name=PED_Sic1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=sbinlab_ib2
#SBATCH --mem=30GB
#SBATCH -t 200:00:00
#SBATCH -o PED_Sic1/out
#SBATCH -e PED_Sic1/err

source /groups/sbinlab/giulio/.bashrc

conda activate hoomd

echo $SLURM_CPUS_PER_TASK

echo $SLURM_CPUS_ON_NODE

python ./simulate.py --seq_name PED_Sic1 --path PED_Sic1