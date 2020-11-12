# Create configuration file where post-QC samples are sorted by sex and tissue.

library(jsonlite)
library(stringr)

# Load JSON file
Post_QC_Quants <- read_json("/scratch/mjpete11/GTEx/Configs/Quantification.config.json")

# Sort post-QC samples by tissue to both by tissue and by sex
# e.g. have list of Amygdala samples, want lists of female/male amygdala samples.
Out <- lapply(Post_QC_Quants[-c(1:3)], function(x){
  a <- unlist(Post_QC_Quants$Females)
  b <- unlist(Post_QC_Quants$Males)
  females <- x[which(x %in% a)]
  males <- x[which(x %in% b)]
  list(Female=females, Male=males)
})

# Convert list-of-lists to json object
Sex_Sorted <- jsonlite::toJSON(Out, pretty = TRUE, auto_unbox = TRUE)

# Write to JSON file 
write(Sex_Sorted, "/scratch/mjpete11/GTEx/Configs/Sex_Sorted_Quantification.config.json")