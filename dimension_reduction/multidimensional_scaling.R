#!/usr/bin/env Rscript

# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold changes between each pair of samples.

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"
METADATA <- file.path(BASE, "data/output/metadata.csv")
COUNTS <- file.path(BASE, "data/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
SEX_DIM12 <- file.path(BASE, "dimension_reduction/MDS_plots/Sex_Dim12.pdf")
RIN_DIM12 <- file.path(BASE, "dimension_reduction/MDS_plots/RIN_Dim12.pdf")
RIN_DIM34 <- file.path(BASE, "dimension_reduction/MDS_plots/RIN_Dim34.pdf")
ISC_DIM12 <- file.path(BASE, "dimension_reduction/MDS_plots/Isc_Dim12.pdf")
ISC_DIM34 <- file.path(BASE, "dimension_reduction/MDS_plots/Isc_Dim34.pdf")

# Load packages
library(data.table)
library(stringr)
library(edgeR)

# Read in files                                                            
meta <- read.table(METADATA, header = TRUE, sep=",", stringsAsFactors=FALSE)
counts <- data.frame(fread(COUNTS))

# Drop gene name and ID from counts df
counts <- counts[,3:ncol(counts)]

# Replace . to - in colnames
colnames(counts) <- str_replace_all(colnames(counts),pattern = "\\.","-")

# Drop samples in metadata that do not have count data
select_samples <- colnames(counts)[colnames(counts)%in%meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples,]

# Set rownames of metadata object equal to sample names
rownames(meta) <- meta$Sample_ID

# Subset counts to only samples present in metadata
counts <- counts[,select_samples]

# Check that the count and meta data have the same samples in the same order
identical(colnames(counts), rownames(meta)) # TRUE

# Make factor indicating sex of samples for filtering
sex <- factor(meta$Sex)

# Make design matrix to use with voom
design <- model.matrix(~meta$Sex)
rownames(design) <- colnames(counts)

# How many genes are there before filtering?
nrow(counts)

# Remove genes with cpm < 1 in each sex
keep <- filterByExpr(counts, design=design, min.count=1, min.prop=0.5)
counts <- counts[keep, ]

# How many genes are left after filtering?
nrow(counts) # 35,610

# limma-voom normalization
counts <- voom(counts, design=design)

# Make list of lists of samples for each tissue
meta$Tissue <- factor(meta$Tissue)

tissue_lst <- list()
for(i in 1:length(levels(meta$Tissue))){
		  tissue_lst[[i]] <- as.vector(meta[meta$Tissue == levels(meta$Tissue)[i], "Sample_ID"])
}

# Rename lists in list as tissue names
names(tissue_lst) <- levels(meta$Tissue)

# Split counts into list of dfs by tissue
tissue_count <- list()
for(i in 1:length(tissue_lst)){
		  tissue_count[[i]] <- counts[,which(colnames(counts) %in% tissue_lst[[i]])]
}
names(tissue_count) <- levels(meta$Tissue)

# metadata split into list of dfs by tissue
meta_lst <- list()
for(i in 1:length(levels(meta$Tissue))){
		  meta_lst[[i]] <- meta[meta$Tissue == levels(meta$Tissue)[i],]
}
names(meta_lst) <- levels(meta$Tissue)

# Check that rownames equals colnames 
Check <- function(a, b){
	all(rownames(a) %in% colnames(b))
}
Res_1 <- Map(Check, a=meta, b=tissue_count)
all(Res_1==TRUE)

# Check that order matches
Match_Check <- function(a, b){
	 Match <- all(rownames(a) == colnames(b))
}
Res_2 <- Map(Match_Check, a=meta, b=tissue_count)
all(Res_2==TRUE)

# MDS plots on one page
# Color samples by sex
sex_colors <- c("blue", "darkgreen")

Sex_Dim12 <- function(DGE, NAME, META, TITLE) {
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(1,2), 
                 col = sex_colors[as.factor(META[['Sex']])],
                 main = NAME)
 		 mtext(TITLE, side=3, outer=TRUE, line=3)
 		 mtext('Dimension 1', side = 1, outer = TRUE, line=1)
  	 	 mtext('Dimension 2', side = 2, outer = TRUE, line=2)
  		 return(plt)
}

pdf(SEX_DIM12)
# margins: c(bottom, left, top, right)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(Sex_Dim12, 
	DGE = tissue_count, 
	NAME = names(tissue_count), 
	META = meta_lst, 
	TITLE='GTEx v8 MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes')
legend(6.0, 2.5, inset=0, legend=c("female", "male"), pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

# Color plots by RIN; dim 1 and 2
RIN_Dim12 <- function(DGE, NAME, META){
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(1,2), 
                 col = 1:length(as.factor(META[['RIN']])),
                 main = NAME)
 		 mtext('GTEx v8 MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes', side=3, outer=TRUE, line=3)
 		 mtext('Dimension 1', side = 1, outer = TRUE, line=1)
  	 	 mtext('Dimension 2', side = 2, outer = TRUE, line=2)
  		 return(plt)
}

# Plot
pdf(RIN_DIM12)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(RIN_Dim12, 
	DGE = tissue_count, 
	NAME = names(tissue_count), 
	META = meta_lst) 
dev.off()

# Color plots by RIN; dim 3 and 4
RIN_Dim34 <- function(DGE, NAME, META){
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(3,4), 
                 col = 1:length(as.factor(META[['RIN']])),
                 main = NAME)
 		 mtext('GTEx v8 MDS Plots: Dimensions 3 and 4; Top 100 Most Variable Genes', side=3, outer=TRUE, line=3)
 		 mtext('Dimension 3', side = 1, outer = TRUE, line=1)
  	 	 mtext('Dimension 4', side = 2, outer = TRUE, line=2)
  		 return(plt)
}

# Plot
pdf(RIN_DIM34)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(RIN_Dim34, 
	DGE = tissue_count, 
	NAME = names(tissue_count), 
	META = meta_lst) 
dev.off()

# Color plots by ischemic time; dim 1 and 2
Isc_Dim12 <- function(DGE, NAME, META){
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(1,2), 
                 col = 1:length(as.factor(META[['Ischemic_Time']])),
                 main = NAME)
 		 mtext('GTEx v8 MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes',
				side=3, outer=TRUE, line=3)
 		 mtext('Dimension 1', side = 1, outer = TRUE, line=1)
  	 	 mtext('Dimension 2', side = 2, outer = TRUE, line=2)
  		 return(plt)
}

# Plot
pdf(ISC_DIM12)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(Isc_Dim12, 
	DGE = tissue_count, 
	NAME = names(tissue_count), 
	META = meta_lst)
dev.off()

# Color plots by ischemic time; dim 3 and 4
Isc_Dim34 <- function(DGE, NAME, META){
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(3,4), 
                 col = 1:length(as.factor(META[['Ischemic_Time']])),
                 main = NAME)
 		 mtext('GTEx v8 MDS Plots: Dimensions 3 and 4; Top 100 Most Variable Genes',
				side=3, outer=TRUE, line=3)
 		 mtext('Dimension 3', side = 1, outer = TRUE, line=1)
  	 	 mtext('Dimension 4', side = 2, outer = TRUE, line=2)
  		 return(plt)
}

# Plot
pdf(ISC_DIM34)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(Isc_Dim34, 
	DGE = tissue_count, 
	NAME = names(tissue_count), 
	META = meta_lst)
dev.off()
