# Histograms for data used in mashr analysis
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Phenotypes")

library(dplyr)
library(ggplot2)
library(MatchIt)

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/Metadata.csv"
PLOT_ALL <- "All_Samples/"
PLOT_55 <- "Above_55/"
AGE_HIST <- "Age_Histograms.pdf"
RACE_HIST <- "Race_Histograms.pdf"
ETH_HIST <- "Ethnicity_Histograms.pdf"

# Organize samples by tissue type into list of dfs for plotting
# Read Metadata CSV.                                                            
Samples <- read.csv(METADATA, header = TRUE, stringsAsFactors=FALSE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample   

# Set tissues column as factor to add levels attribute
Samples$Tissue <- as.factor(Samples$Tissue)

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
      Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Make a subset samples with individuals  >= 55 years old
Old_Meta <- lapply(Meta, function(x) x[x$Age >= 55,])

# List of tissue types
Tissues <- names(Meta)

#-----------------------------------------------------------------------------------------------------
# Density plots
#-----------------------------------------------------------------------------------------------------
# Overlaid density plots for all samples
Overlaid_Density <- function(meta, lab){
    ggplot(meta, aes(x=meta$Age, fill=meta$Sex)) +
    geom_density(alpha=.3) + ggtitle(y) +
    xlab("Age") + ylab("Density") + scale_fill_manual(name = "Sex", values=c("blue", "green"))
}
pdf("test.pdf")
Map(Overlaid_Density, x=Meta, y="Age Distribution of All Samples")
dev.off()

#-----------------------------------------------------------------------------------------------------
# Histograms
#-----------------------------------------------------------------------------------------------------
# Age range across samples: 20-70
min(sapply(Meta, function(x) min(x[['Age']]))) 
max(sapply(Meta, function(x) max(x[['Age']])))

# x- and y-axis labels
x_labels <- lapply(Meta, function(x) seq(min(x[,'Age']), max(x[,'Age']), 5))
y_labels <- seq(0, 30, 1)
# Race and ethinicity y-axis labels
eth_ylab <- seq(0, 100, 2)

# Age histograms
Hist_Func <- function(META, NAMES, TITLE){
  Hist_Plots <- ggplot(META, aes(x=Age, fill=Sex)) + 
    geom_histogram(position = position_dodge(), alpha=0.4, binwidth=1, colour='gray50') +
    scale_y_continuous(labels=y_labels, breaks=y_labels) +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Age") + ylab("Count") 
}

# Plot sample histograms
pdf(paste0(PLOT_ALL,AGE_HIST))
Map(Hist_Func, META = Meta, NAMES = names(Meta), TITLE = "Age Histogram of All Samples")
dev.off()

pdf(paste0(PLOT_55,AGE_HIST))
Map(Hist_Func, META = Old_Meta, NAMES = names(Old_Meta), TITLE = "Age Histogram of Samples >= 55")
dev.off()

# Race histogram plots
Hist_Func <- function(META, NAMES, TITLE){
  META$Race<-as.character(META$Race)
  Hist_Plots <-ggplot(META, aes(x=Race, fill=Sex)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    scale_y_continuous(labels=eth_ylab, breaks=eth_ylab) +
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Race") + ylab("Count")
}

pdf(paste0(PLOT_ALL,RACE_HIST))
Map(Hist_Func, META = Meta, NAMES = names(Meta), TITLE="Race Histogram of All Samples") 
dev.off()

pdf(paste0(PLOT_55,RACE_HIST))
Map(Hist_Func, META = Old_Meta, NAMES = names(Old_Meta), TITLE="Race Histograms of Sample >=55") 
dev.off()

# Ethnicity histogram plots
Hist_Func <- function(META, NAMES, TITLE){
  META$Ethnicity<-as.character(META$Ethnicity)
  Hist_Plots <-ggplot(META, aes(x=Ethnicity, fill=Sex)) +
    geom_bar(stat = "count", position = "dodge") +
    scale_fill_manual(name="Sex", values=c("blue", "green"), labels=c("Female", "Male")) +
    scale_y_continuous(labels=eth_ylab, breaks=eth_ylab) +
    ggtitle(paste(NAMES, TITLE)) +
    xlab("Ethnicity") + ylab("Count")
}

pdf(paste0(PLOT_ALL,ETH_HIST))
Map(Hist_Func, META = Meta, NAMES = names(Meta), TITLE="Ethnicity Histograms of All Samples")
dev.off()

pdf(paste0(PLOT_55,ETH_HIST))
Map(Hist_Func, META = Old_Meta, NAMES = names(Old_Meta), TITLE="Ethnicity Histograms of >=55")
dev.off()
