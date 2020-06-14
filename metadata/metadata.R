#!/usr/bin/env Rscript

# Combine GTEx v8 metadata into single metadata file with information relevant to DE w/ limma
# Sample ID, tissue type, sex, age, RIN (RNA integrity number), post-mortem interval

# Future questions:
# How many samples are replicates? (same tissue samples multiple times from one individual)
# Do we want to drop reps?

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain/metadata"
COUNTS <- file.path(BASE, "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
ATTRIBUTES <- file.path(BASE, "GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
PHENOTYPES <- file.path(BASE, "GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

# Libraries
library(tidyverse)

# Read in files
# Contains: sample ID, tissue type, RIN, post-mortem info
Attributes <- read_tsv(ATTRIBUTES, col_names=TRUE) 
# Contains sex info
Phenotypes <- read_tsv(PHENOTYPES, col_names=TRUE) 

#______________________________________________________________________________________________________
# Combine relevant info into a csv
#______________________________________________________________________________________________________
# Which columns contain relevant info in attributes df?
Attributes[1,] # SAMPID = 1, SMRIN = 5, SMTSD (detailed tissue type) = 7, SMSISCH (post-mort) = 9

# Subset sample ID, tissue type, RIN, and ischemic time
Meta <- Attributes[,c(1,5,7,9)]

# Rename columns
colnames(Meta) <- c("Sample_ID", "RIN", "Tissue", "Ischemic_Time")

# How many samples are in the unfiltered metadata?
nrow(Meta) # 22,951

# Quick view
head(Meta)
tail(Meta)

# How many items are missing in each column from the unfiltered metadata?
# 12.9% samples missing RIN, 1.3% missing ischemic time data
Meta %>%
    select(everything()) %>%
    summarise_all(list(~sum(is.na(.))))/nrow(Meta) * 100

# Add individual ID col;
# grep pattern won't work for the leukemia cell line samples, but I will drop all the cell lines
Meta[["Individual_ID"]] <- str_extract(Meta$Sample_ID, "GTEX-[0-9A-Z]+")

# Drop cell line samples; expression is unstable
Meta <- Meta[!grepl("Cells", Meta$Tissue),]

# How many samples are left?
nrow(Meta) # 22,015

# Reformat tissue names to be contain only '_' between words 
Meta$Tissue <- str_replace_all(Meta$Tissue, c("-" = "", " " = "_", "__" = "_"))
head(Meta)
tail(Meta)

# What are the tissue type categories remaining?
levels(as_factor(Meta$Tissue))

# Get list of female IDs
Fems <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Make list of sex of each individual
Sex <- with(Meta['Individual_ID'], ifelse(Individual_ID %in% Fems, "Female", "Male"))

# Add column containing sex
Meta <- cbind(Meta, "Sex"=Sex)

# Write to file
write.csv(Meta, file.path(BASE, "metadata.csv"))
