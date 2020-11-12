# Write table to see what types of tissue samples are available per individuals, including whole blood samples.

# Paths
METADATA <- file.path("/scratch/mjpete11/GTEx/", "Metadata.csv")
GTEX_SRA_CONFIG <- "/scratch/mjpete11/GTEx/Spinal_Cord/Male_Spinal_Cord.config.json"
WHOLE_BLOOD <- "/data/storage/public/dbgap-8834/whole_blood_rna/whole_blood"
FILE_NAME <- "Blood_Brain_Table.csv"

library(jsonlite)
library(stringr)

# Read in Metadata.csv
Meta <- read.csv(METADATA, header = TRUE)

# Read in config containing dictionary mapping GTEx IDs to Illumina SRA IDs 
Map_Lst <- fromJSON(GTEX_SRA_CONFIG)

# Read in list of whole blood SRA IDs
SRA_Whole_Blood <- list.files(path = WHOLE_BLOOD, pattern = "*_1.fastq", all.files = TRUE)

# Keep only run ID
SRA_Whole_Blood <- str_remove_all(SRA_Whole_Blood, "_1.fastq")

# Make named character vector: Names are GTEx IDs, values are SRA IDs
All_GTEx_IDs <- sapply(Map_Lst, '[[', 9)

# Check that the whole blood samples on agave are in the GTEx sample list
SRA_Whole_Blood %in% All_GTEx_IDs # TRUE

# Get whole blood GTEx IDs
# Order matters for %in% command. The returned boolean vector will be relative to the first in order 
GTEx_Whole_Blood_Wrong_Index <- names(All_GTEx_IDs[(which(SRA_Whole_Blood %in% All_GTEx_IDs))]) # Returns SRA_Whole_Blood index

GTEx_Whole_Blood <- names(All_GTEx_IDs[(which(All_GTEx_IDs %in% SRA_Whole_Blood))]) # Returns All_GTEx_IDs index
 
# Make table
Tmp_Mat <- rbind(Meta[,1:2], cbind(Sample=GTEx_Whole_Blood,Tissue="Blood")) # Drop 'Age' and 'Sex' columns

Tmp_Mat$Sample <- str_extract(Tmp_Mat$Sample,"GTEX-[0-9A-Z]+") # Filter for only individual ID

Sample_Tble <- as.data.frame.matrix(table(Tmp_Mat))

# Write to file
write.csv(Sample_Tble, file = FILE_NAME)


