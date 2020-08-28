#!/usr/bin/env Rscript

# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the 
# largest absolute log-fold changes between each pair of samples.

# Snakemake constants
# Input
METADATA <- snakemake@input[[1]]
COUNTS <- snakemake@input[[2]]

# Output
SEX_DIM12 <- snakemake@output[[1]]
RIN_DIM1 <- snakemake@output[[2]]
RIN_DIM2 <- snakemake@output[[3]]
RIN_DIM3 <- snakemake@output[[4]]
RIN_DIM4 <- snakemake@output[[5]]
ISC_DIM1 <- snakemake@output[[6]]
ISC_DIM2 <- snakemake@output[[7]]
ISC_DIM3 <- snakemake@output[[8]]
ISC_DIM4 <- snakemake@output[[9]]

# Packages
library(data.table)
library(stringr)
library(edgeR)
library(purrr)

# Read in files                                                           
meta <- read.csv(METADATA, header = TRUE, stringsAsFactors = FALSE)
counts <- read.csv(COUNTS, header = TRUE, stringsAsFactors = FALSE)

#_______________________________________________________________________________
# Organize counts and metadata by tissue type 
#_______________________________________________________________________________
# Replaces '.' to '-' in the sample IDs for the projections
colnames(counts) <- str_replace_all(colnames(counts), pattern = "\\.", "-")

# Make list of lists of samples for each tissue
meta$Tissue <- factor(meta$Tissue)

tissue_lst <- list()
for (i in seq_len(length(levels(meta$Tissue)))) {
    tissue_lst[[i]] <- as.vector(meta[meta$Tissue == levels(meta$Tissue)[i], "Sample_ID"])
}

# Rename lists in list as tissue names
names(tissue_lst) <- levels(meta$Tissue)

# Split counts into list of dfs by tissue
tissue_count <- list()
for (i in seq_len(length(tissue_lst))) {
    tissue_count[[i]] <- counts[, which(colnames(counts)
                                             %in% tissue_lst[[i]])]
}
names(tissue_count) <- levels(meta$Tissue)

# metadata split into list of dfs by tissue
meta_lst <- list()
for (i in seq_len(length(levels(meta$Tissue)))) {
    meta_lst[[i]] <- meta[meta$Tissue == levels(meta$Tissue)[i], ]
}
names(meta_lst) <- levels(meta$Tissue)

# Check that rownames equals colnames
check <- function(a, b) {
    all(rownames(a) %in% colnames(b))
}
res_1 <- Map(check, a=meta, b=tissue_count)
all(res_1 = TRUE)

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
