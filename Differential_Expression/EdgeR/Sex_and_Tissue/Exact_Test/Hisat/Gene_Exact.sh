#!/bin/bash 

#SBATCH -p fn1                         # fat node; more memmory                                        
#SBATCH -n 28                          # number of cores                
#SBATCH -t 2:00:00                     # wall time                                       
#SBATCH --job-name=Hisat_Gene_Exact    # job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=mjpete11@asu.edu   # send-to address

module purge

module load r/3.5.2

Rscript Gene_Matched_Exact.R
Rscript Gene_AgeMatched_Exact.R
