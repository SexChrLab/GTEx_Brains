#!/usr/bin/env Rscript

# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold changes between each pair of samples.

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"
METADATA <- file.path(BASE, "data/output/metadata.csv")
COUNTS <- file.path(BASE, "data/input/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
MDS <- file.path(BASE, "dimension_reduction/MDS_plots/MDS.pdf")

# Load packages
library(data.table)
#library(readr)
library(stringr)
#library(dplyr)
#library(limma)
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

# check that the count and meta data have the same samples in the same order
identical(colnames(counts), rownames(meta)) # TRUE

# make factor indicating sex of samples for filtering
sex <- factor(meta$Sex)

# make design matrix to use with voom
design <- model.matrix(~meta$Sex)
rownames(design) <- colnames(counts)

# Remove genes with cpm < 1 in each sex
counts <- DGEList(counts, group=sex)
keep <- filterByExpr(counts, design=design, min.count=1, min.prop=0.5)
counts <- counts[keep, ,keep.lib.sizes=FALSE]

# limma-voom normalization
counts <- voom(counts, design=design)

# Make list of lists of samples for each tissue
meta$Tissue <- factor(meta$Tissue)

Tissue_Lst <- list()
for(i in 1:length(levels(meta$Tissue))){
		  Tissue_Lst[[i]] <- as.vector(meta[meta$Tissue == levels(meta$Tissue)[i], "Sample_ID"])
}

# Rename lists in list as tissue names
names(Tissue_Lst) <- levels(meta$Tissue)

# Split counts into list of dfs by tissue
Tissue_Count <- list()
for(i in 1:length(Tissue_Lst)){
		  Tissue_Count[[i]] <- counts[,which(colnames(counts) %in% Tissue_Lst[[i]])]
}
names(Tissue_Count) <- levels(meta$Tissue)

# metadata split into list of dfs by tissue
meta_Lst <- list()
for(i in 1:length(levels(meta$Tissue))){
		  meta_Lst[[i]] <- meta[meta$Tissue == levels(meta$Tissue)[i],]
}
names(meta_Lst) <- levels(meta$Tissue)

# Check that rownames equals colnames 
Check <- function(a, b){
	all(rownames(a) %in% colnames(b))
}
Res_1 <- Map(Check, a=meta, b=Tissue_Count)
all(Res_1==TRUE)

# Check that order matches
Match_Check <- function(a, b){
	 Match <- all(rownames(a) == colnames(b))
}
Res_2 <- Map(Match_Check, a=meta, b=Tissue_Count)
all(Res_2==TRUE)

# MDS plots on one page
colors <- c("blue", "darkgreen")

MDS_Plot <- function(DGE, NAME, META, TITLE) {
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(1,2), 
                 col = colors[META[['Sex']]],
                 main = NAME)
 		 mtext(TITLE, side=3, outer=TRUE, line=3)
 		 mtext('Dimension 1', side = 1, outer = TRUE, line=1)
  	 	 mtext('Dimension 2', side = 2, outer = TRUE, line=2)
  		 return(plt)
}
pdf(MDS)
# margins: c(bottom, left, top, right)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(MDS_Plot, 
	DGE = Tissue_Count, 
	NAME = names(Tissue_Count), 
	META = meta_Lst, 
	TITLE='GTEx v8 MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes')
# female will be blue because color in plot/legend assigned alphabetically
legend(10.0, 2.5, inset=0, legend=levels(meta$Sex), pch=16, cex=2.0, col=colors, xpd=NA)
dev.off()

