#!/bin/bash

#SBATCH -q tempboost
#SBATCH -p fn1 
#SBATCH --job-name=MDS_plots
#SBATCH -o slurm.%j.out                
#SBATCH -e slurm.%j.err                
#SBATCH --mail-type=END,FAIL           
#SBATCH --mail-user=mjpete11@asu.edu 
#SBATCH -t 10:0:00

# Export functions so they will be available in subshells
source ~/miniconda3/etc/profile.d/conda.sh

conda activate Differential_Expression

Rscript multidimensional_scaling.R
