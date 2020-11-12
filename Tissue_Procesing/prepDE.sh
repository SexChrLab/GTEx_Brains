#!/bin/bash

 #SBATCH -N 1
 #SBATCH -n 28
 #SBATCH -t 00:02:00
 #SBATCH --job-name=prepDE
 #SBATCH -o slurm.%j.out
 #SBATCH -e slurm.%j.err
 #SBATCH --mail-type=END,FAIL
 #SBATCH --mail-user=mjpete11@asu.edu

# iterate through dirs and apply prepDE.py to file in dir begining with dir name
for i in $(ls /scratch/mjpete11/GTEx/Tissue_Procesing/) # this will try to iteerate through all objects in dir
do
    dir=$(readlink -f $i) # iterate through symbolic dir links
    file=${i%/*} # match until backslash
    python2.7 prepDE.py -i $dir/"$file"_Files.txt -g $dir/"$file"_gene_matrix.csv -t $dir/"$file"_transcript_matrix.csv # gene and transcript output matrix paths
done
