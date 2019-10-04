# Histograms of sample phenotypes labelled by sex
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes")

library(dplyr)
library(ggplot2)
library(MatchIt)

# Constants
PLOT_DIR <- "/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes"
AGE_HIST <- "/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes/SexAgeRaceEth_Matched/Age_Histograms.pdf"
RACE_HIST <- "/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes/SexAgeRaceEth_Matched/Race_Histograms.pdf"
ETH_HIST <- "/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes/SexAgeRaceEth_Matched/Ethnicity_Histograms.pdf"

# Read sample metadata.                                                            
Samples <- read.csv(file.path("/scratch/mjpete11/GTEx/Metadata", "Metadata.csv"), header = TRUE)

#-----------------------------------------------------------------------------------------------------
# Create config with age-matched, equal sample sizes across tissues
#-----------------------------------------------------------------------------------------------------
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

# Organize samples by tissue type into list of dfs for plotting
Meta_Samples <- list()
for(i in 1:length(levels(df.Samples$Tissue))){
  Meta_Samples[[i]] <- df.Samples[df.Samples$Tissue == levels(df.Samples$Tissue)[i],]
}
names(Meta_Samples) <- levels(df.Samples$Tissue)

#-----------------------------------------------------------------------------------------------------
# Density plots
#-----------------------------------------------------------------------------------------------------
# Overlaid density plots for all samples
ggplot(df.Samples, aes(x=df.Samples$Age, fill=df.Samples$Sex)) +
  geom_density(alpha=.3) + ggtitle("Density Plot of Sample Age Distribution") +
  xlab("Age") + ylab("Density") + scale_fill_manual(name = "Sex", values=c("blue", "green"))

#-----------------------------------------------------------------------------------------------------
# Histograms
#-----------------------------------------------------------------------------------------------------
# Function to plot histograms on top of each other
x_axis_labels <- seq(min(df.Samples[,'Age']), max(df.Samples[,'Age']), 5)
y_axis_labels <- seq(0, 30, 1)

Hist_Func <- function(META, NAMES, TITLE){
  Hist_Plots <- ggplot(META, aes(x=Age, fill=Sex)) + 
    geom_histogram(position = position_dodge(), alpha=0.4, binwidth=1, colour='gray50') +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Age") + ylab("Count") 
}

# Plot sample histograms
pdf(AGE_HIST)
Map(Hist_Func, META = Meta_Samples, NAMES = names(Meta_Samples), TITLE = "Sample Age Histogram")
dev.off()

# Density plots by tissue
Density_Func <- function(META, NAMES, TITLE){
  Density_Plots <- ggplot(META, aes(x=Age, fill=Sex)) +
    geom_density(data=subset(META), alpha=.3) + 
    scale_fill_manual(name = "Sex", values=c("blue", "green")) +  
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Age") + ylab("Density") 
}

# Plot samples
Map(Density_Func, META = Meta_Samples, NAMES = names(Meta_Samples), TITLE = "Sample Age Density Plot")

# Race histogram plots
y_axis_labels <- seq(0, 100, 2)

Hist_Func <- function(a, b){
  a$Race<-as.character(a$Race)
  Hist_Plots <-ggplot(a, aes(x=Race, fill=Sex)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    ggtitle(paste(b, "Sample Race Histogram")) +
    xlab("Race") + ylab("Count")
}

pdf(RACE_HIST)
Map(Hist_Func, a = Meta_Samples, b = names(Meta_Samples)) 
dev.off()

# Ethnicity histogram plots
y_axis_labels <- seq(0, 100, 2) 

Hist_Func <- function(a, b){
  a$Ethnicity<-as.character(a$Ethnicity)
  Hist_Plots <-ggplot(a, aes(x=Ethnicity, fill=Sex)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    scale_y_continuous(labels=y_axis_labels, breaks=y_axis_labels) +
    ggtitle(paste(b, "Sample Ethnicity Histogram")) +
    xlab("Ethnicity") + ylab("Count")
}

pdf(ETH_HIST)
Map(Hist_Func, a = Meta_Samples, b = names(Meta_Samples))
dev.off()

