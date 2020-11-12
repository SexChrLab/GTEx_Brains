# This script is for plotting PCA dimensions 1:4 for the top 1000 most variable genes in a tissue vs tissue approach.

METADATA = file.path("/scratch/mjpete11/GTEx/", "Toy_Metadata.csv")
COUNTS = file.path("/scratch/mjpete11/GTEx/Data_Exploration/", "Toy_Count_Matrix.tsv")
PLOT_DIR = "PCA_Tissue_Vs_Tissue_Plots/"
FILE_NAME = "PCA_Tissue_Vs_Tissue.pdf"

# Load packages                                                                 
library(readr)
library(stringr)
library(ggplot2)
library(dplyr)

# Read Metadata CSV.                                                            
Samples = read.csv(METADATA, header = TRUE)

# Create every combination of tissue matrices
# Make list of lists of samples for each tissue
Tissue_Lst <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Tissue_Lst[[i]] <- as.vector(Samples[Samples$Tissue == levels(Samples$Tissue)[i], "Sample"])
}
# Rename lists in list as tissue names
names(Tissue_Lst) <- levels(Samples$Tissue)

# Read in counts 
cts <- read.csv(COUNTS, sep = "\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Filter for most variable genes
Var_Genes <- apply(cts, 1, var) # 1 is to  set margin to rows

# Select 1000 most variable genes
Select_Var <- names(sort(Var_Genes, decreasing=TRUE)[1:1000]) 

# Subset the most variable genes
Highly_Variable_CPM <- cts[Select_Var,]
dim(Highly_Variable_CPM)

# Create list of data frames containing counts for each tissue
Var_Count_Lst <- list()
for(i in 1:length(Tissue_Lst)){
  Var_Count_Lst[[i]] <- Highly_Variable_CPM[,which(colnames(Highly_Variable_CPM) %in% Tissue_Lst[[i]])]
}

# Rename lists in list
names(Var_Count_Lst) <- paste("count_",names(Tissue_Lst),sep="")

# Create matrix of all combos of matrices
Combos <- combn(1:length(Var_Count_Lst), m = 2, simplify = T) # m = number of elements to choose

# Create list of matrices of combinations of all tissue counts 
Combo_Count <- list()
for(i in 1:ncol(Combos)){
  Combo_Count[[i]] <- t(cbind(Var_Count_Lst[[Combos[,i][1]]], Var_Count_Lst[[Combos[,i][2]]])) # Create combo of matrices by iterating through all possible combinations
  names(Combo_Count)[i] <- paste(names(Tissue_Lst)[c(Combos[,i])], collapse = "-") # Rename by combining all tissue names
}

# Euclidean distances between the rows
Combo_Dist <- lapply(Combo_Count, dist)

# Calculate PCA
Fit <- lapply(Combo_Dist, function(y){
  cmdscale(y, eig=T, k=4) # Set number of dimensions
})

# Add metadata to PCA data
# Get metadata for each tissue and store dataframes in list
Meta <- lapply(levels(Samples$Tissue), function(x){
  data.frame(Samples[which(Samples$Tissue==x),])
})
names(Meta) <- levels(Samples$Tissue)

# Create list of metadata of each tissue combination
Combos <- combn(1:length(levels(Samples$Tissue)),m = 2,simplify = T)
Meta_lst <- list()
for(i in 1:ncol(Combos)){
  Meta_Lst[[i]] <- rbind(Meta[[Combos[,i][1]]], Meta[[Combos[,i][2]]])
  names(Meta_Lst)[i] <- paste(names(Meta)[c(Combos[,i])], collapse = "-")
}

# Merge metadata and PCA data into single dataframe
# e.g. combine two lists of dataframes into a list of merged pairs of dfs
PCA <- list()
for(i in 1:length(Fit)){
  # make data frame with samples as a column
  tmpA <- data.frame(Sample=rownames(Fit[[i]][[1]]), Fit[[i]][[1]])
  # save meta as a tmp dataframe
  tmpB <- Meta[[i]]
  # refactor tmpB data frame to avoid errors
  tmpB$Sample <- factor(tmpB$Sample)
  # join tables together by samples name
  PCA[[i]] <- left_join(tmpA,tmpB)
  # match unnamed columns
  tmpC <- str_which(colnames(PCA[[i]]) ,pattern = "^[A-Z][0-9]+")
  # rename columns
  colnames(PCA[[i]])[tmpC] <- paste("Dim_", 1:length(tmpC), sep = "")
}
names(PCA) <- names(Meta)

# PCA plots
setwd(PLOT_DIR)
pdf(FILE_NAME)

# Function to plot PCA for dimensions 1 and 2, color by tissue and shape by sex
colors <- c("blue", "darkgreen")
shapes <-  c(21,22)

PCA_k2_Sex <- function(a, b){
  Plot_k2_Sex <- ggplot(a, aes(x=Dim_1, y=Dim_2, shape=Tissue, color=Sex, label=Tissue)) +
                 geom_point(size=2, aes(fill = Sex)) + ggtitle(paste("PCA plot", b, sep = " : " )) +
                 scale_shape_manual(values=shapes) + scale_color_manual(values=colors) +
                 scale_fill_manual(values=colors)  
}

Map(PCA_k2_Sex, a = PCA, b = names(PCA))

# Function to plot PCA for dimensions 1 and 2, color by age and shape by tissue
PCA_k2_Age <- function(a, b){
    Plot_k2_Age <- ggplot(a, aes(x=Dim_1, y=Dim_2, color=Age, shape=Tissue, label=Tissue)) +
                   geom_point(size=2) + ggtitle(paste("PCA plot", b, sep = " : " )) + 
                   scale_color_gradient(low="black", high="lightblue")
}

Map(PCA_k2_Age, a = PCA, b = names(PCA))

# Function to plot PCA for dimensions 3 and4, color by tissue and shape by sex
PCA_k4_Sex <- function(a, b){
  Plot_k4_Sex <- ggplot(a, aes(x=Dim_3, y=Dim_4, shape=Tissue, color=Sex, label=Tissue)) +
                 geom_point(size=2, aes(fill = Sex)) + ggtitle(paste("PCA plot", b, sep = " : " )) +
                 scale_shape_manual(values=shapes) + scale_color_manual(values=colors) +
                 scale_fill_manual(values=colors)  
}

Map(PCA_k4_Sex, a = PCA, b = names(PCA))

# Function to plot PCA for dimensions 1 and 2, color by age and shape by tissue
PCA_k4_Age <- function(a, b){
  Plot_k4_Age <- ggplot(a, aes(x=Dim_3, y=Dim_4, color=Age, shape=Tissue, label=Tissue)) +
                 geom_point(size=2) + ggtitle(paste("PCA plot", b, sep = " : " )) + 
                 scale_color_gradient(low="black", high="lightblue")
}

Map(PCA_k4_Age, a = PCA, b = names(PCA))

dev.off()


