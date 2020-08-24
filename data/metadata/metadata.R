#!/usr/bin/env Rscript

# Combine GTEx v8 metadata into single metadata file with information relevant to DE w/ limma
# Sample ID, tissue type, sex, age, RIN (RNA integrity number), post-mortem interval

# Future questions:
# How many samples are replicates? (same tissue samples multiple times from one individual)
# Do we want to drop reps?

# Snakemake constants
# Input
ATTRIBUTES <- snakemake@input[[1]]
PHENOTYPES <- snakemake@input[[2]]
COUNTS <- snakemake@input[[3]] 

# Output
TABLE1 <- snakemake@output[[1]]
TABLE2 <- snakemake@output[[2]]
META <- snakemake@output[[3]]
OLD_META <- snakemake@output[[4]]

# Libraries
library(tidyverse)
library(data.table)
library(stringr)

# Read in files
# Contains: sample ID, tissue type, RIN, post-mortem info
attribs <- read_tsv(ATTRIBUTES, col_names=TRUE) 
# Contains sex info
phenotypes <- read_tsv(PHENOTYPES, col_names=TRUE) 
gene_counts <- data.frame(fread(COUNTS))

#______________________________________________________________________________________________________
# metadata preprocessing 
#______________________________________________________________________________________________________
# Which columns contain relevant info in attribs df?
attribs[1,] # SAMPID = 1, SMRIN = 5, SMTSD (detailed tissue type) = 7, SMSISCH (post-mort) = 9

# Subset sample ID, tissue type, RIN, and ischemic time
meta <- attribs[,c(1,5,7,9)]

# Rename columns
colnames(meta) <- c("Sample_ID", "RIN", "Tissue", "Ischemic_Time")

# Add individual ID col;
# grep pattern won't work for the leukemia cell line samples, but I will drop all the cell lines
meta[["Individual_ID"]] <- str_extract(meta$Sample_ID, "GTEX-[0-9A-Z]+")

# Reformat tissue names to be contain only '_' between words 
meta$Tissue <- str_replace_all(meta$Tissue, c("-" = "", " " = "_", "__" = "_"))

# Replace . to - in colnames
colnames(gene_counts) <- str_replace_all(colnames(gene_counts), pattern = "\\.","-")

# Get list of female IDs
fems <- phenotypes$SUBJID[which(phenotypes$SEX==2)]

# Make list of sex of each individual
sex <- with(meta['Individual_ID'], ifelse(Individual_ID %in% fems, "Female", "Male"))

# Add column containing sex
meta <- cbind(meta, "Sex"=sex)

# Add column containing age (only decade intervals are publically accessible)
meta$Age <- phenotypes$AGE[match(meta$Individual_ID, phenotypes$SUBJID)]

# Rearrange column order
meta <- meta %>% select(Individual_ID, Sex, Age, Tissue, Sample_ID, Ischemic_Time, RIN)

# Drop samples in metadata that do not have count data
select_samples <- colnames(gene_counts)[colnames(gene_counts) %in% meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples, ]

# Subset sample replicates
Reps <- meta %>%
		    group_by(Individual_ID, Tissue) %>%
		    filter(n() > 1) %>%
			ungroup()

# Samples minus replicates and non-brain_meta tissue 
brain_meta <- meta[grepl("Brain", meta$Tissue),]

# Remove the tissue replciates (Cerebellum/Cortex are replicates of Cerebellar/Frontal_Cortex)
brain_meta <- brain_meta[!grepl("Cerebellum|Brain_Cortex", brain_meta$Tissue),]

# Subset samples >= 50; re-do once I get the metadata with the exact ages and not just decade intervals
old_brain_meta <- brain_meta %>% filter(Age >= 50)

# Number of samples in meta
nrow(brain_meta) # 2,146
nrow(old_brain_meta) # 1,831

# Function to summarise percent missing values
Missing <- function(x){
		   x %>%
		 	 select(everything()) %>%
		  	 summarise_all(list(~sum(is.na(.))))/nrow(x) * 100
}

# How many samples have missing values?
Missing(meta)

# Drop rows missing values
meta <- meta %>% drop_na()

# Check that all samples missing values were removed
Missing(meta)

# Quick view
head(brain_meta)
tail(brain_meta)

# Summary stats: how many tissues per sex
Summary_Stat <- function(x){
	x %>% group_by(Tissue, Sex) %>% tally()
}

# In the full sample set
print(Summary_Stat(brain_meta), n=22)
print(Summary_Stat(old_brain_meta), n=22)

#______________________________________________________________________________________________________
# Write to file 
#______________________________________________________________________________________________________
# summary stats
write.csv((as.data.frame(Summary_Stat(brain_meta))), TABLE1)
write.csv((as.data.frame(Summary_Stat(old_brain_meta))), TABLE2)

# metadata files
write.csv(brain_meta, META, row.names=FALSE)
write.csv(old_brain_meta, OLD_META, row.names=FALSE)
