# Write configuration file 

library(jsonlite)
library(stringr)

# Read in orginal config 
Quant_Files <- read_json("/scratch/mjpete11/GTEx/Configs/Dup_Quantification.config.json")

# Identify sample replicates
Filt_Quant <- Quant_Files[-c(1:3)];

Reps <- lapply(Filt_Quant, function(y){
  c <- y[which(duplicated(str_extract(y,"GTEX-[0-9A-Z]+")))]
  if (length(c)==0) {NA} else {c}
})

# Write config
write_json(Reps,"/scratch/mjpete11/GTEx/Configs/Rep_Only_Quantification.config.json")
