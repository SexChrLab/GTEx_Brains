#-----------------------------------------------------------------------------------------------------
# Create config with age-matched, equal sample sizes across tissues
#-----------------------------------------------------------------------------------------------------
setwd("/scratch/mjpete11/GTEx/Metadata")

# Subset only older individuals
df.Samples <- subset(Samples, Age >= 55)

# Find which tissue has smallest sample size before matching
table(df.Samples$Tissue[which(df.Samples$Sex=='Female')]) # 18 females above 55 in Spinal_Cord/Substantia_Nigra

# Add Group column with sex as logical
df.Samples$Group <- as.logical(df.Samples$Sex == 'Female')

# Match male and females exactly (i.e equal ratios of male and females) exactly for 
# tissue and age; include race and ethincity for propensity score
match.it <- matchit(Group ~ Tissue + Age + Race + Ethnicity, data = df.Samples, method = "nearest", exact = c("Age","Tissue"))
df.Samples <- match.data(match.it)[1:ncol(df.Samples)]

# Sort by sex and age
df.Samples <- df.Samples[with(df.Samples, order(Tissue, Sex, Age)),]

# Find which tissue has smallest sample size 
table(df.Samples$Tissue[which(df.Samples$Sex=='Female')]) # 11 females spinal cord samples

# For each df in list, subset 11 female/male samples; descending by age
df.Samples <- do.call(rbind, lapply(unique(df.Samples$Tissue), function(w){
  tmp <- df.Samples[which(df.Samples$Tissue==w), ]
  tmp <- tmp[order(tmp$Age,decreasing = T), ]
  rbind(tmp[which(tmp$Sex=="Male"),][1:11,], tmp[which(tmp$Sex=="Female"), ][1:11, ])
}))

# Check sample ages
min(df.Samples$Age) # 55
max(df.Samples$Age) # 69

# Write to file
write.csv(df.Samples, file='/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes/SexAgeRaceEth_Matched/Matched_Metadata.csv', row.names=FALSE)
