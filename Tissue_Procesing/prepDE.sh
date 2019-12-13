#!/bin/bash

#SBATCH -N 1
#SBATCH -n 28
#SBATCH -t 2:00:00 
#SBATCH --job-name=prepDE
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mjpete11@asu.edu

python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Anterior/Anterior_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Caudate/Caudate_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Cerebellar/Cerebellar_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Cerebellum/Cerebellum_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Cortex/Cortex_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Frontal_Cortex/Frontal_Cortex_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Hippocampus/Hippocampus_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Hypothalamus/Hypothalamus_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Nucleus_Accumbens/Nucleus_Accumbens_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Putamen/Putamen_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Spinal_Cord/Spinal_Cord_Files.txt
python2.7 prepDE.py -i /scratch/mjpete11/GTEx/Tissue_Procesing/Substantia_Nigra/Substantia_Nigra_Files.txt
