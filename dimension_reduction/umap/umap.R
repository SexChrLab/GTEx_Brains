#!/usr/bin/env Rscript

# Script to plot the umap projections

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/expression_matrices/output/processed_counts_with_sex_chr.csv")
PROJECTION <- file.path(BASE, "dimension_reduction/umap/projections/umap_with_sex_chr.txt")

# Output
PLOT1 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_tissue.pdf")
PLOT2 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_sex.pdf")
PLOT3 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_age.pdf")
PLOT4 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_isc.pdf")
PLOT5 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_rin.pdf")

# Read in files
meta <- read.table(METADATA, header = TRUE, sep = ",", stringsAsFactors = FALSE)
#counts <- read.csv(COUNTS)
a <- read.table(PROJECTION, sep = " ")

#_______________________________________________________________________________
# Sanity check 
#_______________________________________________________________________________
# Number of samples in meta
nrow(meta) # 2,146

# Number of samples in projection
nrow(counts) # 13,955

# Convert gene_counts (EList object) to matrix for use w/ umap
counts <- data.matrix(counts)

#_______________________________________________________________________________
# umap
#_______________________________________________________________________________
# gene_counts is an expression matrix: genes as rows and samples as columns.
start_t <- Sys.time()
a <- umap(t(counts), n_neighbors = 50, min_dist = 0.5)
end_t <- Sys.time()
print(paste("Elapsed time:", end_t - start_t))

# Write umap object to file
write.table(a$layout, file = PROJECTION, sep = " ")

# Read in projection
#a <- read.table(PROJECTION, sep = " ")
head(a)
nrow(a)
colnames(a) <- c('X1', 'X2')

# extract the projections (a$layout) and join it with metadata (particularly the 'meta.Tissue' column)
#full.covariates is a vector of model covariate names (in case you want to plot by those as well)
a.umap = data.frame(a, meta$Tissue, meta$Sex, meta$Age, meta$RIN, meta$Ischemic_Time)

#_______________________________________________________________________________
# Plots 
#_______________________________________________________________________________
# Function to plot color points by one feature
umap_plot <- function(x, y){
	p = ggplot(a.umap, aes_string('X1','X2', color=x)) +
		geom_point(size=0.5) +
	#	scale_color_manual(values=region.colors) +
		theme_classic(base_size=18) +
		xlab('UMAP 1') +
		ylab('UMAP 2') +
		coord_fixed() +
		theme(axis.ticks=element_blank(), axis.text=element_blank())
	p <- p + labs(fill = y) 
	p <- p + guides(shape = guide_legend(override.aes = list(size = 3))) 
	p <- p + guides(color = guide_legend(override.aes = list(size = 3))) 
	p <- p + theme(legend.title = element_text(size = 8), 
				   legend.text = element_text(size = 6))
	return(p)
}
# Color by tissue
ggsave(umap_plot(x = 'meta.Tissue', y = 'Tissue'), file = PLOT1, useDingbats = FALSE)

# Color by sex
ggsave(umap_plot(x = 'meta.Sex', y= 'Sex'), file = PLOT2, useDingbats = FALSE)

# Color by age
ggsave(umap_plot(x = 'meta.Age', y = 'Age'), file = PLOT3, useDingbats = FALSE)

# Color by ischemic time
ggsave(umap_plot(x = 'meta.Ischemic_Time', y = 'Ischemic time'), file = PLOT4, useDingbats = FALSE)

# Color by RIN
ggsave(umap_plot(x = 'meta.RIN', y = 'RIN'), file = PLOT5, useDingbats = FALSE)

##################### Kenny's example code ##################################
## Redo the plot, this time without highlighting batch
#p = ggplot(a.umap,aes_string('X1','X2',color='meta.Tissue')) +
#	geom_point(size=1.25) +
##	scale_color_manual(values=region.colors) +
#	theme_classic(base_size=18) +
#	xlab('UMAP 1') +
#	ylab('UMAP 2') +
#	coord_fixed() +
#	theme(axis.ticks=element_blank(),axis.text=element_blank())
#ggsave(p, file = PLOT2, useDingbats = FALSE)
#
## Redo the plot, this time showing color AND shape to better differentiate points since many colors look similar
#p = ggplot(a.umap,aes_string('X1','X2',color='meta.Tissue',shape='meta.Tissue')) +
#	geom_point(size=1.25) +
##	scale_color_manual(values=region.colors) +
##	scale_shape_manual(values=region.shapes) +
#	theme_classic(base_size=18) +
#	xlab('UMAP 1') +
#	ylab('UMAP 2') +
#	coord_fixed() +
#	guides(shape = guide_legend(override.aes = list(size = 3))) +
#	theme(axis.ticks=element_blank(),axis.text=element_blank())
#ggsave(p, file = PLOT3, useDingbats = FALSE)
