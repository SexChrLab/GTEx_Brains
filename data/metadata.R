#!/usr/bin/env Rscript

# Combine GTEx v8 metadata into single metadata file with information relevant to DE w/ limma
# Sample ID, tissue type, sex, age, RIN (RNA integrity number), post-mortem interval

# Future questions:
# How many samples are replicates? (same tissue samples multiple times from one individual)
# Do we want to drop reps?

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain/data"
COUNTS <- file.path(BASE, "input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
ATTRIBUTES <- file.path(BASE, "input/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
PHENOTYPES <- file.path(BASE, "input/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

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
head(Meta)
tail(Meta)

# Get list of female IDs
Fems <- Phenotypes$SUBJID[which(Phenotypes$SEX==2)]

# Make list of sex of each individual
Sex <- with(Meta['Individual_ID'], ifelse(Individual_ID %in% Fems, "Female", "Male"))

# Add column containing sex
Meta <- cbind(Meta, "Sex"=Sex)

# Add column containing age (only decade intervals are publically accessible)
Meta$Age <- Phenotypes$AGE[match(Meta$Individual_ID, Phenotypes$SUBJID)]

# Rearrange column order
Meta <- Meta %>% select(Individual_ID, Sex, Age, Tissue, Ischemic_Time, RIN)

# Quick view
head(Meta)
tail(Meta)

#______________________________________________________________________________________________________
# Full samples set (no filtering)
#______________________________________________________________________________________________________
# How many samples are in the unfiltered metadata?
nrow(Meta) # 22,951

# Function to summarise percent missing values
Missing <- function(x){
		   x %>%
		 	 select(everything()) %>%
		  	 summarise_all(list(~sum(is.na(.))))/nrow(x) * 100
}

#______________________________________________________________________________________________________
# Samples minus cell lines
#______________________________________________________________________________________________________
# Drop cell line samples; expression is unstable
Meta <- Meta[!grepl("Cells", Meta$Tissue),]

# 22,015 non-unique samples
nrow(Meta) 

# What are the tissue type categories remaining?
levels(as_factor(Meta$Tissue))

# ISC: 0.28%, RIN: 13.33%
Missing(Meta)

#______________________________________________________________________________________________________
# Samples replicates minus cell lines
#______________________________________________________________________________________________________
# Subset sample replicates
Reps <- Meta %>%
		    group_by(Individual_ID, Tissue) %>%
		    filter(n() > 1) %>%
			ungroup()

# ISC: 0.90%, RIN: 43.48% 
Missing(Reps)

# Summarize how many f/m samples have replicates per tissue
Summary_Reps <- Reps %>% 
					group_by(Individual_ID, Tissue, Sex) %>%
					tally()

#______________________________________________________________________________________________________
# Samples minus cell lines and sample replicates
#______________________________________________________________________________________________________
Meta <- Meta %>% distinct(Individual_ID, Tissue, .keep_all=TRUE)

# 17,827 unique samples
nrow(Meta)

# ISC: 0.08%, RIN: 5.85%
Missing(Meta)

#______________________________________________________________________________________________________
# Samples minus cell lines, sample replicates, and individuals < 55 years old
#______________________________________________________________________________________________________
Old_Meta <- Meta %>% filter(Age >= 50)

# 12,413 unique samples
nrow(Old_Meta)

# ISC: 0.04%, RIN: 5.76%
Missing(Old_Meta)

#______________________________________________________________________________________________________
# Write to file the last two sample sets
#______________________________________________________________________________________________________
write.csv(Meta, file.path(BASE, "output/metadata.csv"), row.names=FALSE)
write.csv(Old_Meta, file.path(BASE, "output/older_metadata.csv"), row.names=FALSE)
