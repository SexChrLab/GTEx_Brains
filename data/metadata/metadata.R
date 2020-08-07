#!/usr/bin/env Rscript

# Combine GTEx v8 metadata into single metadata file with information relevant to DE w/ limma
# Sample ID, tissue type, sex, age, RIN (RNA integrity number), post-mortem interval

# Future questions:
# How many samples are replicates? (same tissue samples multiple times from one individual)
# Do we want to drop reps?

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain/data"
COUNTS <- file.path(BASE, "counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
ATTRIBUTES <- file.path(BASE, "metadata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
PHENOTYPES <- file.path(BASE, "metadata/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# Libraries
library(tidyverse)

# Read in files
# Contains: sample ID, tissue type, RIN, post-mortem info
Attributes <- read_tsv(ATTRIBUTES, col_names=TRUE) 
# Contains sex info
Phenotypes <- read_tsv(PHENOTYPES, col_names=TRUE) 

#______________________________________________________________________________________________________
# Metadata preprocessing 
#______________________________________________________________________________________________________
# Which columns contain relevant info in attributes df?
Attributes[1,] # SAMPID = 1, SMRIN = 5, SMTSD (detailed tissue type) = 7, SMSISCH (post-mort) = 9

# Subset sample ID, tissue type, RIN, and ischemic time
Meta <- Attributes[,c(1,5,7,9)]

# Rename columns
colnames(Meta) <- c("Sample_ID", "RIN", "Tissue", "Ischemic_Time")

# Add individual ID col;
# grep pattern won't work for the leukemia cell line samples, but I will drop all the cell lines
Meta[["Individual_ID"]] <- str_extract(Meta$Sample_ID, "GTEX-[0-9A-Z]+")

# Reformat tissue names to be contain only '_' between words 
Meta$Tissue <- str_replace_all(Meta$Tissue, c("-" = "", " " = "_", "__" = "_"))

# Get list of female IDs
Fems <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Make list of sex of each individual
Sex <- with(Meta['Individual_ID'], ifelse(Individual_ID %in% Fems, "Female", "Male"))

# Add column containing sex
Meta <- cbind(Meta, "Sex"=Sex)

# Add column containing age (only decade intervals are publically accessible)
Meta$Age <- Phenotypes$AGE[match(Meta$Individual_ID, Phenotypes$SUBJID)]

# Rearrange column order
Meta <- Meta %>% select(Individual_ID, Sex, Age, Tissue, Sample_ID, Ischemic_Time, RIN)

# Drop rows missing values
Meta <- Meta %>% drop_na()

# Function to summarise percent missing values
Missing <- function(x){
		   x %>%
		 	 select(everything()) %>%
		  	 summarise_all(list(~sum(is.na(.))))/nrow(x) * 100
}

# Check that all samples missing values were removed
Missing(Meta)

# Quick view
head(Meta)
tail(Meta)

# Subset sample replicates
Reps <- Meta %>%
		    group_by(Individual_ID, Tissue) %>%
		    filter(n() > 1) %>%
			ungroup()

# Samples minus replicates and non-brain tissue 
Brain <- Meta[grepl("Brain", Meta$Tissue),]

# Remove the tissue replciates (Cerebellum/Cortex are replicates of Cerebellar/Frontal_Cortex)
Brain <- Brain[!grepl("Cerebellum|Brain_Cortex", Brain$Tissue),]

# Subset samples >= 50; re-do once I get the metadata with the exact ages and not just decade intervals
Old_Brain <- Brain %>% filter(Age >= 50)

# 3,192 samples
nrow(Brain)

# 2,701
nrow(Old_Brain)

# Summary stats: how many tissues per sex
Summary_Stat <- function(x){
	x %>% group_by(Tissue, Sex) %>% tally()
}

# In the full sample set
print(Summary_Stat(Brain), n=22)
print(Summary_Stat(Old_Brain), n=22)

#______________________________________________________________________________________________________
# Write to file 
#______________________________________________________________________________________________________
write.csv(Brain, file.path(BASE, "metadata/metadata.csv"), row.names=FALSE)
write.csv(Old_Brain, file.path(BASE, "metadata/older_metadata.csv"), row.names=FALSE)
