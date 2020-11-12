# Create config of samples that passed QC without sample replicates.
# Results sorted by tissue only.

library(jsonlite)
library(stringr)

Quant_Files <- read_json("/scratch/mjpete11/GTEx/Configs/Dup_Quantification.config.json")

# Check if samples have same individual ID
# If same ID, keep only one.
No_Reps <- lapply(Quant_Files, function(y){
    c <- y[-which(duplicated(str_extract(y,"GTEX-[0-9A-Z]+")))]
    if (length(c)==0) {y} else {c}
})

# Write config
write_json(No_Reps,"/scratch/mjpete11/GTEx/Configs/Quantification.config.json")
