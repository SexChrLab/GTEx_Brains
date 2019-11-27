#!/usr/bin/Rscript  

#SBATCH -p phi                                                                 
#SBATCH -n 28                          # Max is 256 on phi node                
#SBATCH -t 48:00:00                                                           
#SBATCH --job-name=Exact_DGX        # Job name
#SBATCH -o slurm.%j.out                # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                # STDERR (%j = JobId)
#SBATCH --mail-type=END,FAIL           # notifications for job done & fail
#SBATCH --mail-user=mjpete11@asu.edu   # send-to address

module purge

module load r/3.5.2

Rscript Gene_Matched_Exact.R
Rscript Gene_AgeMatched_Exact.R
Rscript Trans_Matched_Exact.R
Rscript Trans_AgeMatched_Exact.R
