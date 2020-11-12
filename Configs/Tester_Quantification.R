# This script is to write a table sumamrizing how many tissue samples (of different regions) come from the same individuals.  

CONFIG = "/scratch/mjpete11/GTEx/Configs/Quantification.config.json"
METADATA = "/scratch/mjpete11/GTEx/Metadata.csv"
FILE_NAME = "/scratch/mjpete11/GTEx/Configs/Sample_Table.csv"

library(plyr)

# Keep only tissue, and drop Cerebellar and Frontal Cortex.
Quant <- fromJSON(CONFIG)

Samples_Quant <- Quant[-c(1:3)]; # Keep only lists of tissue samples

Samples_Unique <- Samples_Quant[-c(4,7)] # drop Cerebellar and Frontal Cortex

names(Samples_Unique)

# Make vector of unique individual IDs by sex.
Fem_Lst <- unlist(unique(str_extract(Quant$Females,"GTEX-[0-9A-Z]+"))) # 59 unique individuals

Male_Lst <- unlist(unique(str_extract(Quant$Males,"GTEX-[0-9A-Z]+")))  # 143 unique individuals

# Function that returns how many individuals are present in x tissues
# Only 4 females have samples for all 11 tissues
lapply(Fem_Lst, function(x){
  tmp <- lapply(Samples_Unique, function(y){
    unique(str_extract(y, "GTEX-[0-9A-Z]+"))==x
  })
  if(sum(sapply(tmp, sum))==length(Samples_Unique)){return(x)} 
})

# No male individuals have samples for all 11 tissues
lapply(Male_Lst, function(x){
  tmp <- lapply(Samples_Unique, function(y){
    unique(str_extract(y, "GTEX-[0-9A-Z]+"))==x
  })
  if(sum(sapply(tmp, sum))==length(Samples_Unique)){return(x)} 
})

# Make a table of how frequent each female individual is across tissues
Res_Fem <- sapply(Fem_Lst, function(x){
           sum(str_extract(unlist(Samples_Unique),"GTEX-[0-9A-Z]+")==x)
         })

# Same for males 
Res_Male <- sapply(Male_Lst, function(x){
  sum(str_extract(unlist(Samples_Unique),"GTEX-[0-9A-Z]+")==x)
})

# Combine results into single dataframe
Female <- t(data.frame(lapply(Res_Fem, type.convert), stringsAsFactors=FALSE))

Male <- t(data.frame(lapply(Res_Male, type.convert), stringsAsFactors=FALSE))

# Add column indicating sex
Append_df <- function(df_Names) {
  do.call(rbind, lapply(df_Names, function(x) {
    cbind(get(x), source = x)
  }))
}

Res_df <- as.data.frame(Append_df(c("Female", "Male")))

rownames(Res_df) <- str_replace_all(rownames(Res_df),pattern = "\\.","-") # Replace . to - in sample names

# Rename columns
colnames(Res_df) <- c("Count", "Sex")

# Write to file
write.csv(Res_df, FILE_NAME)
