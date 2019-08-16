#!/bin/bash

#SBATCH -t 48:00:00                                                           
#SBATCH --job-name=GLM_Ratio_Test_Sex_Tissue         
#SBATCH -o slurm.%j.out               
#SBATCH -e slurm.%j.err               
#SBATCH --mail-type=END,FAIL           
#SBATCH --mail-user=mjpete11@asu.edu  

module purge

module load r/3.5.2

Rscript Age_Matched.R 
