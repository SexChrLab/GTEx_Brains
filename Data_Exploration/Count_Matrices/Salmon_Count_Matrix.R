# Create count matrix from salmon results. 
# Note: if complaint about select(), unload dplyr
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices")

FILE = "/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/Salmon/Transcript_Salmon_CountMatrix.tsv"

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures) 
library(AnnotationDbi)
library(edgeR)   
library(readr)

# Read Metadata CSV.                                                            
samples = read.csv(file.path("/scratch/mjpete11/GTEx/Metadata/", "Metadata.csv"), header = TRUE)

# GTEX-13N2G-0011-R2a-SM-5MR4Q (male: substantia nigra) is missing the quant file. 
# Remove from list of samples
samples <- samples[-c(1288),]

# Set rownames of metadata object equal to sample names.                        
rownames(samples) <- samples$Sample   

# Set path to quant.sf files and check that all files are there.                
files = file.path("/scratch/mjpete11/GTEx/All_Salmon_Quants", samples$Sample, "quant.sf")
head(files)                                                                     
all(file.exists(files)) # TRUE

# Set the names of the file paths object equal to the sample names.             
names(files) <- samples$Sample                                                  

# Make tx2gene table to map gene IDs to transcript IDs.                         
TxDb <- makeTxDbFromGFF(file = "/scratch/mjpete11/GTEx/gencode.v29.annotation.gtf")
k <- keys(TxDb, keytype = "TXNAME")                                             
tx2gene <- select(TxDb, k, "GENEID", "TXNAME")                                  
head(tx2gene)                                                                   

# Import salmon output files.
# If you want the list of genes, add: txOut = FALSE. To get list of isoforms, add: txOut = TRUE
Txi <- tximport(files, type = "salmon", tx2gene = tx2gene, txOut = TRUE)      
cts <- Txi$counts
df_cts <- data.frame(cts, check.names = FALSE)

# Write to TSV
write.table(df_cts, file = FILE, sep = "\t", row.names = TRUE, quote = FALSE)

