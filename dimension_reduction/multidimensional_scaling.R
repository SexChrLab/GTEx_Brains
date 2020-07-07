#!/usr/bin/env Rscript

# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-fold changes between each pair of samples.

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
CHRX <- file.path(BASE, "data/gene_annotation/gencodeGenes_Xchr.txt")
CHRY <- file.path(BASE, "data/gene_annotation/gencodeGenes_Ychr.txt")

# Output
SEX_DIM12 <- file.path(BASE, "dimension_reduction/MDS_plots/Sex_Dim12.pdf")
RIN_DIM12 <- file.path(BASE, "dimension_reduction/MDS_plots/RIN_Dim12.pdf")
RIN_DIM34 <- file.path(BASE, "dimension_reduction/MDS_plots/RIN_Dim34.pdf")
ISC_DIM12 <- file.path(BASE, "dimension_reduction/MDS_plots/Isc_Dim12.pdf")
ISC_DIM34 <- file.path(BASE, "dimension_reduction/MDS_plots/Isc_Dim34.pdf")

# Packages
library(data.table)
library(stringr)
library(edgeR)
library(colorRamps)

# Read in files                                                            
xchr <- read.table(CHRX, sep = "")
ychr <- read.table(CHRY, sep = "")
meta <- read.table(METADATA, header = TRUE, sep=",", stringsAsFactors=FALSE)
gene_counts <- data.frame(fread(COUNTS))

#_______________________________________________________________________________
# Drop X and Y-linked genes to check if clustering is driven by sex-linked loci
#_______________________________________________________________________________
# Original number of genes
nrow(gene_counts) # 56,200

# Remove X and Y-linked genes
rows_to_drop <- intersect(xchr$V6, gene_counts$Name)
rows_to_drop <- c(rows_to_drop, intersect(ychr$V6, gene_counts$Name))

gene_counts <- gene_counts[-which(gene_counts[,1] %in% rows_to_drop),]

# Are the expected number of genes remaining?; Yes
nrow(gene_counts) # 53,325 
nrow(xchr) + nrow(ychr) # 2,875

#_______________________________________________________________________________
# Sample pre-processing
#_______________________________________________________________________________
# Drop gene name and ID from gene_counts df
gene_counts <- gene_counts[,3:ncol(gene_counts)]

# Replace . to - in colnames
colnames(gene_counts) <- str_replace_all(colnames(gene_counts),pattern = "\\.","-")

# Drop samples in metadata that do not have count data
select_samples <- colnames(gene_counts)[colnames(gene_counts)%in%meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples,]

# Set rownames of metadata object equal to sample names
rownames(meta) <- meta$Sample_ID

# Subset gene_counts to only samples present in metadata
gene_counts <- gene_counts[,select_samples]

# Check that the count and meta data have the same samples in the same order
identical(colnames(gene_counts), rownames(meta)) # TRUE

#_______________________________________________________________________________
# Filter genes by expression in each sex and voom normalize 
#_______________________________________________________________________________
# Make factor indicating sex of samples for filtering
sex <- factor(meta$Sex)

# Make design matrix to use with voom
design <- model.matrix(~meta$Sex)
rownames(design) <- colnames(gene_counts)

# How many genes are there before filtering?
nrow(gene_counts) # 53,325

# Remove genes with cpm < 1 in each sex
keep <- filterByExpr(gene_counts, design=design, min.count=1, min.prop=0.5)
gene_counts <- gene_counts[keep, ]

# How many genes are left after filtering?
nrow(gene_counts) # 34,341 

# limma-voom normalization
gene_counts <- voom(gene_counts, design=design)

#_______________________________________________________________________________
# Organize gene_counts and metadata by tissue type 
#_______________________________________________________________________________
# Make list of lists of samples for each tissue
meta$Tissue <- factor(meta$Tissue)

tissue_lst <- list()
for(i in 1:length(levels(meta$Tissue))){
	tissue_lst[[i]] <- as.vector(meta[meta$Tissue == 
								 levels(meta$Tissue)[i], "Sample_ID"])
}

# Rename lists in list as tissue names
names(tissue_lst) <- levels(meta$Tissue)

# Split gene_counts into list of dfs by tissue
tissue_count <- list()
for(i in 1:length(tissue_lst)){
	tissue_count[[i]] <- gene_counts[,which(colnames(gene_counts) %in% tissue_lst[[i]])]
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

#_______________________________________________________________________________
# Plot the PCs 1 through 4 vs RIN and ischemic time  
#_______________________________________________________________________________
# Get the PCs 1 through 4 by generating an MDS object with plotMDS

# I don't think I need any of the plot related parameters...
MDS_obj <- plotMDS(tissue_count[[1]],
               gene.selection = "common",
               top = 100, 
               dim.plot = c(1,2), 
			   plot = FALSE) # Do not write to graphics device

class(MDS_obj)

# Extract PC 1
head(MDS_obj$x)

# Make matrix: rownames are the samples and cols are their corresponding RIN and PC1 
# Get the RIN values and sample IDs for only the Amygdala samples
mat <- cbind(MDS_obj$x, meta_lst[[1]]$RIN)
rownames(mat) <- meta_lst[[1]]$Sample_ID
colnames(mat) <- c("PC1", "RIN")
head(mat)

sex_colors <- c("blue", "darkgreen")
pdf("test.pdf")
plot(mat,
   	pch = 16, 
	cex = 1, 
    col = sex_colors[as.factor(meta[['Sex']])],
    main = "Do I need a title?")
dev.off()

# Repeat, but scale up 
# Function to generate MDS objects containing PCs
do_MDS <- function(DGE){
	obj <- plotMDS(DGE,
				   gene.selection = "common",
				   top = 100,
				   dim.plot = c(1,2),
				   plot = FALSE) # Do not write to graphics device
	return(obj)
}
MDS_lst <- lapply(tissue_count, do_MDS)

length(MDS_lst)
head(MDS_lst[[1]])

# Function to make matrices for plotting
# make list of matrices; rows are samples from same tissue type;
# cols are their correspinding RIN and PC1 values
make_matrix <- function(MDS_OBJ, META){
	mat <- cbind(MDS_OBJ$x, META$RIN)
	rownames(mat) <- META$Sample_ID
	colnames(mat) <- c("PC1", "RIN")
	return(mat)
}
mat_lst <- Map(make_matrix, MDS_OBJ=MDS_lst, META=meta_lst)

length(mat_lst)
head(mat_lst[[11]])

# Function to plot matrices
sex_colors <- c("blue", "darkgreen")
pdf("test.pdf")
plot_func <- function(){
	plot(mat,
	   	pch = 16, 
		cex = 1, 
   	    col = sex_colors[as.factor(meta[['Sex']])],
    main = "Do I need a title?")
dev.off()



#_______________________________________________________________________________
# MDS plots on one page
#_______________________________________________________________________________
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
	TITLE='Sex MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes')
legend(6.0, 2.5, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

# Color plots by RIN; dim 1 and 2
plot_colors <- colorRampPalette(c('green','blue'))

RIN_Dim12 <- function(DGE, NAME, META){
  plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(1,2), 
				 col = plot_colors(length(unique((META[['RIN']])/max(META[['RIN']])))),
                 main = NAME)
 		 mtext('RIN Value MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes', 
				side=3, outer=TRUE, line=3)
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
				 col = plot_colors(length(unique((META[['RIN']])/max(META[['RIN']])))),
                 main = NAME)
 		 mtext('RIN value MDS Plots: Dimensions 3 and 4; Top 100 Most Variable Genes', 
				side=3, outer=TRUE, line=3)
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
				 col = plot_colors(length(unique((META[['Ischemic_Time']])/max(META[['Ischemic_Time']])))),
                 main = NAME)
 		 mtext('Ischemic Time MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes',
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
				 col = plot_colors(length(unique((META[['Ischemic_Time']])/max(META[['Ischemic_Time']])))),
                 main = NAME)
 		 mtext('Ischemic Time MDS Plots: Dimensions 3 and 4; Top 100 Most Variable Genes',
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
