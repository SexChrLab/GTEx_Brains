# Use tximport to convert stringtie transcript counts to gene counts 
# and to store counts as single matrix
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Hisat/")

# Constants
PATHS <- c('/scratch/mjpete11/GTEx/Amygdala/Hisat_Stringtie/gene_count_matrix.csv', 
           '/scratch/mjpete11/GTEx/Anterior/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Caudate/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Cerebellar/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Cerebellum/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Cortex/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Frontal_Cortex/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Hippocampus/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Hypothalamus/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Nucleus_Accumbens/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Putamen/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Spinal_Cord/Hisat_Stringtie/gene_count_matrix.csv',
           '/scratch/mjpete11/GTEx/Substantia_Nigra/Hisat_Stringtie/gene_count_matrix.csv')

FILES <- c("Amygdala_Gene_Hisat_CountMatrix.tsv", "Anterior_Gene_Hisat_CountMatrix.tsv", 
                "Caudate_Gene_Hisat_CountMatrix.tsv", "Cerebellar_Gene_Hisat_CountMatrix.tsv",
                "Cerebellum_Gene_Hisat_CountMatrix.tsv", "Cortex_Gene_Hisat_CountMatrix.tsv",
                "FrontalCortex_Gene_Hisat_CountMatrix.tsv", "Hippocampus_Gene_Hisat_CountMatrix.tsv",
                "Hypothalamus_Gene_Hisat_CountMatrix.tsv", "NucleusAccumbens_Gene_Hisat_CountMatrix.tsv",
                "Putamen_Gene_Hisat_CountMatrix.tsv", "SpinalCord_Gene_Hisat_CountMatrix.tsv",
                "SubstantiaNigra_Gene_Hisat_CountMatrix.tsv")

FILES.1 <- c("Amygdala_Transcript_Hisat_CountMatrix.tsv", "Anterior_Transcript_Hisat_CountMatrix.tsv", 
             "Caudate_Transcript_Hisat_CountMatrix.tsv", "Cerebellar_Transcript_Hisat_CountMatrix.tsv",
             "Cerebellum_Transcript_Hisat_CountMatrix.tsv", "Cortex_Transcript_Hisat_CountMatrix.tsv",
             "FrontalCortex_Transcript_Hisat_CountMatrix.tsv", "Hippocampus_Transcript_Hisat_CountMatrix.tsv",
             "Hypothalamus_Transcript_Hisat_CountMatrix.tsv", "NucleusAccumbens_Transcript_Hisat_CountMatrix.tsv",
             "Putamen_Transcript_Hisat_CountMatrix.tsv", "SpinalCord_Transcript_Hisat_CountMatrix.tsv",
             "SubstantiaNigra_Transcript_Hisat_CountMatrix.tsv")

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr)

# Read Metadata CSV.                                                            
samples = read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Metadata.csv"), header = TRUE)

# Samples missing t_data.ctab file
# Drop GTEX-ZAB4-0011-R4a-SM-4SOKB; Amydgala male (53) 
# Drop GTEX-13N2G-0011-R2a-SM-5MR4Q Substantia_Nigra Male (1288)
samples <- samples[-c(53, 1288),]

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample      

# List of samples by tissue type
Tissues <- as.character(unique(samples$Tissue))

Tissue_Lst <- list()
for (i in Tissues){
  Tissue_Lst[[i]] <- as.character(samples$Sample[samples$Tissue==i])
}

# Set path to dirs with t_data.ctab files and check that all files are there.   
Paths <- lapply(Tissue_Lst, function(x) file.path('/scratch/mjpete11/GTEx/All_Hisat_Quants', x, "t_data.ctab") )
all(lapply(Paths, function(x) all(file.exists(x)))==TRUE) # TRUE

# Set the names of the file paths object equal to the sample names.   
# Note: This renames obj in global env; Return var must have same name as original var (e.g. obj <- setNames(obj, x))
Name_Vec <- function(w, x){
  names(w) <- samples$Sample[samples$Tissue==x]
  return(w)
}
Paths <- Map(Name_Vec, w=Paths, x=Tissues)

# Read transcript ID and gene ID cols from one tsv from one sample from each tissue type.
# Map of transcript to gene IDs will be diff for each tissue since they were run seperately.
# e.g. different novel transcripts/diff number of transcripts will be detected in diff tissues.
tmp <- lapply(Paths, function(x) read_tsv(x[[1]]))
tx2gene <- lapply(tmp, function(x) x[, c("t_name", "gene_name")])
head(tx2gene)

# Import Hisat output files.
# If you want the list of genes, add: txOut = FALSE. To get list of isoforms, add: txOut = TRUE
Tximport_Func <- function(w, x){
  res <- tximport(w, type = "stringtie", tx2gene = x, txOut = TRUE)
  return(res)
}
Txi <- Map(w=Paths, x=tx2gene, Tximport_Func)

cts <- lapply(Txi, function(x) x$counts)
df_lst <- lapply(cts, function(x) data.frame(x, check.names=FALSE))

# Write to TSV
Write_Func <- function(w, x){
  write.table(w, file=x, sep = "\t", row.names = TRUE, quote = FALSE)
}
Map(Write_Func, w=df_lst, x=FILES.1)
