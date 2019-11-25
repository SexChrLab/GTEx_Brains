#!/bin/bash
#SBATCH --job-name=Salmon_Amygdala # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=mjpete11@asu.edu # send-to address
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 5:00:00
#SBATCH --qos=normal

source activate salmon_environment

snakemake --snakefile Quantification.snakefile -j 20 --keep-target-files --rerun-incomplete --cluster "sbatch -n 8 -c 1 -t 5:00:00"

