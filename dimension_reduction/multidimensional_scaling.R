#!/usr/bin/env Rscript

# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the 
# largest absolute log-fold changes between each pair of samples.

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"
PLOT <- "dimension_reduction/New_MDS_Plots"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
CHRX <- file.path(BASE, "data/gene_annotation/gencodeGenes_Xchr.txt")
CHRY <- file.path(BASE, "data/gene_annotation/gencodeGenes_Ychr.txt")

# Output
SEX_DIM12 <- file.path(BASE, PLOT, "Sex_Dim12.pdf")
RIN_DIM1 <- file.path(BASE, PLOT, "RIN_Dim1.pdf")
RIN_DIM2 <- file.path(BASE, PLOT, "RIN_Dim2.pdf")
RIN_DIM3 <- file.path(BASE, PLOT, "RIN_Dim3.pdf")
RIN_DIM4 <- file.path(BASE, PLOT, "RIN_Dim4.pdf")
ISC_DIM1 <- file.path(BASE, PLOT, "Isc_Dim1.pdf")
ISC_DIM2 <- file.path(BASE, PLOT, "Isc_Dim2.pdf")
ISC_DIM3 <- file.path(BASE, PLOT, "Isc_Dim3.pdf")
ISC_DIM4 <- file.path(BASE, PLOT, "Isc_Dim4.pdf")

# Packages
library(data.table)
library(stringr)
library(edgeR)
library(colorRamps)
library(purrr)

# Read in files                                                           
xchr <- read.table(CHRX, sep = "")
ychr <- read.table(CHRY, sep = "")
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
gene_counts <- data.frame(fread(COUNTS))

#_______________________________________________________________________________
# Drop X and Y-linked genes to check if clustering is driven by sex-linked loci
#_______________________________________________________________________________
# Original number of genes
nrow(gene_counts) # 56,200

# Remove X and Y-linked genes
rows_to_drop <- intersect(xchr$V6, gene_counts$Name)
rows_to_drop <- c(rows_to_drop, intersect(ychr$V6, gene_counts$Name))

gene_counts <- gene_counts[-which(gene_counts[, 1] %in% rows_to_drop), ]

# Are the expected number of genes remaining?; Yes
nrow(gene_counts) # 53,325
nrow(xchr) + nrow(ychr) # 2,875

#_______________________________________________________________________________
# Sample pre-processing
#_______________________________________________________________________________
# Drop gene name and ID from gene_counts df
gene_counts <- gene_counts[, 3:ncol(gene_counts)]

# Replace . to - in colnames
colnames(gene_counts) <- str_replace_all(colnames(gene_counts), pattern = "\\.","-")

# Drop samples in metadata that do not have count data
select_samples <- colnames(gene_counts)[colnames(gene_counts) %in% meta$Sample_ID]
meta <- meta[meta$Sample_ID %in% select_samples, ]

# Set rownames of metadata object equal to sample names
rownames(meta) <- meta$Sample_ID

# Subset gene_counts to only samples present in metadata
gene_counts <- gene_counts[, select_samples]

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
keep <- filterByExpr(gene_counts, design = design, min.count = 1, min.prop = 0.5)
gene_counts <- gene_counts[keep, ]

# How many genes are left after filtering?
nrow(gene_counts) # 34,341

# limma-voom normalization
gene_counts <- voom(gene_counts, design = design)

#_______________________________________________________________________________
# Organize gene_counts and metadata by tissue type 
#_______________________________________________________________________________
# Make list of lists of samples for each tissue
meta$Tissue <- factor(meta$Tissue)

tissue_lst <- list()
for (i in 1:seq_len(levels(meta$Tissue))) {
    tissue_lst[[i]] <- as.vector(meta[meta$Tissue ==
                                 levels(meta$Tissue)[i], "Sample_ID"])
}

# Rename lists in list as tissue names
names(tissue_lst) <- levels(meta$Tissue)

# Split gene_counts into list of dfs by tissue
tissue_count <- list()
for (i in 1:seq_len(tissue_lst)) {
    tissue_count[[i]] <- gene_counts[, which(colnames(gene_counts)
                                             %in% tissue_lst[[i]])]
}
names(tissue_count) <- levels(meta$Tissue)

# metadata split into list of dfs by tissue
meta_lst <- list()
for (i in 1:seq_len(levels(meta$Tissue))) {
    meta_lst[[i]] <- meta[meta$Tissue == levels(meta$Tissue)[i], ]
}
names(meta_lst) <- levels(meta$Tissue)

# Check that rownames equals colnames
check <- function(a, b) {
    all(rownames(a) %in% colnames(b))
}
res_1 <- Map(check, a=meta, b=tissue_count)
all(res_1= = TRUE)

# Check that order matches
match_check <- function(a, b){
	match <- all(rownames(a) == colnames(b))
}
res_2 <- Map(match_check, a = meta, b = tissue_count)
all(res_2 = TRUE)

#_______________________________________________________________________________
# MDS by tissue; samples labeled by sex 
#_______________________________________________________________________________
sex_colors <- c("blue", "darkgreen")

Sex_Dim12 <- function(DGE, NAME, META){
	plt <- plotMDS(DGE,
                 gene.selection = "common",
                 top = 100, 
                 pch = 16, 
                 cex = 1, 
                 dim.plot = c(1,2), 
                 col = sex_colors[as.factor(META[['Sex']])],
                 main = NAME)
 	       mtext('Sex MDS Plots: Dimensions 1 and 2; Top 100 Most Variable Genes', 
				  side=3, outer=TRUE, line=3)
  	   	   mtext('Dimension 1', side = 1, outer = TRUE, line=1)
   	 	   mtext('Dimension 2', side = 2, outer = TRUE, line=2)
	return(plt)
}

pdf(SEX_DIM12)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
Map(Sex_Dim12, 
	DGE = tissue_count, 
	NAME = names(tissue_count), 
	META = meta_lst) 
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
   	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

#_______________________________________________________________________________
# RIN and Ischemic time vs PC1:4
#_______________________________________________________________________________
# Function to generate MDS objects containing PCs
make_MDS <- function(DGE, DIM){
	obj <- plotMDS(DGE,
       			   gene.selection = "common",
				   top = 100,
				   dim.plot = DIM,
				   plot = FALSE) # Do not write to graphics device
	return(obj)
}
# MDS object with PCs 1 and 2
PC_12<- pmap(list(DGE=tissue_count), make_MDS, DIM=c(1,2))
# MDS object with PCs 3 and 4
PC_34<- pmap(list(DGE=tissue_count), make_MDS, DIM=c(3,4))

# Function to make matrices for plotting
# make list of matrices; rows are samples from same tissue type;
# cols are their correspinding RIN and PC1 values
make_matrix <- function(MDS_OBJ, META, COLS, PC, VAL){
	mat <- cbind(MDS_OBJ[[PC]], META[[VAL]])
	rownames(mat) <- META$Sample_ID
	colnames(mat) <- COLS 
	return(mat)
}
# List of matrices of PC1 and RIN value for each sample of same tissue type
mat_lst1 <- pmap(list(MDS_OBJ=PC_12, META=meta_lst), make_matrix, 
				 COLS=c("PC1", "RIN"), PC="x", VAL="RIN")
# List of matrices of PC2 and RIN value for each sample of same tissue type
mat_lst2 <- pmap(list(MDS_OBJ=PC_12, META=meta_lst), make_matrix, 
				 COLS=c("PC2", "RIN"), PC="y", VAL="RIN")
# List of matrices of PC3 and RIN value for each sample of same tissue type
mat_lst3 <- pmap(list(MDS_OBJ=PC_34, META=meta_lst), make_matrix, 
				 COLS=c("PC3", "RIN"), PC="x", VAL="RIN")
# List of matrices of PC4 and RIN value for each sample of same tissue type
mat_lst4 <- pmap(list(MDS_OBJ=PC_34, META=meta_lst), make_matrix, 
				 COLS=c("PC4", "RIN"), PC="y", VAL="RIN")

# List of matrices of PC1 and ischemic time for each sample of same tissue type
mat_lst5 <- pmap(list(MDS_OBJ=PC_12, META=meta_lst), make_matrix, 
				 COLS=c("PC1", "RIN"), PC="x", VAL="Ischemic_Time")
# List of matrices of PC2 and ischemic time for each sample of same tissue type
mat_lst6 <- pmap(list(MDS_OBJ=PC_12, META=meta_lst), make_matrix, 
				 COLS=c("PC2", "RIN"), PC="y", VAL="Ischemic_Time")
# List of matrices of PC3 and ischemic time for each sample of same tissue type
mat_lst7 <- pmap(list(MDS_OBJ=PC_34, META=meta_lst), make_matrix, 
				 COLS=c("PC3", "RIN"), PC="x", VAL="Ischemic_Time")
# List of matrices of PC4 and ischemic time for each sample of same tissue type
mat_lst8 <- pmap(list(MDS_OBJ=PC_34, META=meta_lst), make_matrix, 
				 COLS=c("PC4", "RIN"), PC="y", VAL="Ischemic_Time")

# Function to plot matrices
plot_func <- function(MAT, META, NAME, TITLE, XAXIS, YAXIS){
	plt <- plot(MAT,
			   	pch = 16, 
				cex = 1, 
   	   	   	    col = sex_colors[as.factor(META[['Sex']])],
	       	    main = NAME)
			mtext(TITLE, side=3, outer=TRUE, line=3)
		    mtext(XAXIS, side = 1, outer = TRUE, line=1)
			mtext(YAXIS, side = 2, outer = TRUE, line=2)
	return(plt)
}

# Plot
pdf(RIN_DIM1)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst1,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 1 vs RIN by tissue type",
    XAXIS = "Dimension 1",
	YAXIS = "RIN")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(RIN_DIM2)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst2,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 2 vs RIN by tissue type",
    XAXIS = "Dimension 2",
	YAXIS = "RIN")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(RIN_DIM3)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst3,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 3 vs RIN by tissue type",
    XAXIS = "Dimension 3",
	YAXIS = "RIN")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(RIN_DIM4)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst4,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 4 vs RIN by tissue type",
    XAXIS = "Dimension 4",
	YAXIS = "RIN")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(ISC_DIM1)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst5,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 1 vs ischemic time by tissue type",
    XAXIS = "Dimension 1",
	YAXIS = "Ischemic time")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(ISC_DIM2)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst6,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 2 vs ischemic time by tissue type",
    XAXIS = "Dimension 2",
	YAXIS = "Ischemic time")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(ISC_DIM3)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst7,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 3 vs ischemic time by tissue type",
    XAXIS = "Dimension 3",
	YAXIS = "Ischemic time")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()

pdf(ISC_DIM4)
par(mfrow = c(4, 3), cex=0.4, mar = c(3, 2, 2, 2), oma =c(5, 5, 6, 2), xpd=TRUE) 
pmap(list(MAT = mat_lst8,
	META = meta_lst,
	NAME = names(tissue_count)),
	plot_func, 
	TITLE = "Dimension 4 vs ischemic time by tissue type",
    XAXIS = "Dimension 4",
	YAXIS = "Ischemic time")
legend(4.0, 4.0, inset=0, legend=c("female", "male"), 
	   pch=16, cex=2.0, col=sex_colors, xpd=NA)
dev.off()
