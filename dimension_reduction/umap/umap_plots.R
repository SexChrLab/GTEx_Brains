#!/usr/bin/env Rscript

# Plot umap projection with and without sex chr

# Load libraries
library(umap)
library(ggplot2)
library(stringr)
library(edgeR)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"

# Input
METADATA <- snakemake@input[[1]] 
SEX_PROJ <- snakemake@input[[2]]
NO_SEX_PROJ <- snakemake@input[[3]]
 
# Output
PLOT1 <- snakemake@output[[1]] 
PLOT2 <- snakemake@output[[2]]
PLOT3 <- snakemake@output[[3]]
PLOT4 <- snakemake@output[[4]]
PLOT5 <- snakemake@output[[5]]
PLOT6 <- snakemake@output[[6]]
PLOT7 <- snakemake@output[[7]]
PLOT8 <- snakemake@output[[8]]
PLOT9 <- snakemake@output[[9]]
PLOT10 <- snakemake@output[[10]]
PLOT11 <- snakemake@output[[11]]
PLOT12 <- snakemake@output[[12]]

# Read in files
meta <- read.csv(METADATA, header = TRUE, stringsAsFactors = FALSE)
sex_proj <- read.csv(SEX_PROJ, header = TRUE, stringsAsFactors = FALSE)
no_sex_proj <- read.csv(NO_SEX_PROJ, header = TRUE, stringsAsFactors = FALSE)

#______________________________________________________________________________
# Plots 
#______________________________________________________________________________
# Rename cols in projection dfs
colnames(sex_proj) <- c("Sample_ID", "UMAP_1", "UMAP_2")
colnames(no_sex_proj) <- c("Sample_ID", "UMAP_1", "UMAP_2")

# Add metadata to projections for plotting; bind using sample ID
sex_proj <- merge(sex_proj, meta, by = "Sample_ID")
no_sex_proj <- merge(no_sex_proj, meta, by = "Sample_ID")

# Function to plot color points by one feature
umap_plot <- function(PROJ, FEATURE){
	p = ggplot(PROJ, aes_string("UMAP_1", "UMAP_2", color = FEATURE)) +
		geom_point(size = 0.5) +
		theme_classic(base_size = 18) +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		coord_fixed() +
		theme(axis.ticks=element_blank(), axis.text=element_blank())
	p <- p + labs(fill =  FEATURE) 
	p <- p + guides(shape = guide_legend(override.aes = list(size = 3))) 
	p <- p + guides(color = guide_legend(override.aes = list(size = 3))) 
	p <- p + theme(legend.title = element_text(size = 8), 
				   legend.text = element_text(size = 6))
	return(p)
}

# Print plots
# +/- sex chr, color by tissue
ggsave(umap_plot(PROJ = sex_proj, FEATURE = "Tissue"), file = PLOT1) 
ggsave(umap_plot(PROJ = no_sex_proj, FEATURE = "Tissue"), file = PLOT2) 

# +/- sex chr, color by sex
ggsave(umap_plot(PROJ = sex_proj, FEATURE = "Sex"), file = PLOT3) 
ggsave(umap_plot(PROJ = no_sex_proj, FEATURE = "Sex"), file = PLOT4) 

# +/- sex chr, color by age
ggsave(umap_plot(PROJ = sex_proj, FEATURE = "Age"), file = PLOT5) 
ggsave(umap_plot(PROJ = no_sex_proj, FEATURE = "Age"), file = PLOT6) 

# +/- sex chr, color by ischemic time
ggsave(umap_plot(PROJ = sex_proj, FEATURE = "Ischemic_Time"), file = PLOT7) 
ggsave(umap_plot(PROJ = no_sex_proj, FEATURE = "Ischemic_Time"), file = PLOT8) 

# +/- sex chr, color by RIN
ggsave(umap_plot(PROJ = sex_proj, FEATURE = "RIN"), file = PLOT9) 
ggsave(umap_plot(PROJ = no_sex_proj, FEATURE = "RIN"), file = PLOT10) 

# Get rid of legengs when coloring by individual id
umap_plot <- function(PROJ, FEATURE){
	p = ggplot(PROJ, aes_string("UMAP_1", "UMAP_2", color = FEATURE)) +
		geom_point(size = 0.5) +
		theme_classic(base_size = 18) +
		xlab("UMAP 1") +
		ylab("UMAP 2") +
		coord_fixed() +
		theme(legend.title = element_blank(), legend.position = "none")
	return(p)
}

# +/- sex chr, color by individual
ggsave(umap_plot(PROJ = sex_proj, FEATURE = "Individual_ID"), file = PLOT11) 
ggsave(umap_plot(PROJ = no_sex_proj, FEATURE = "Individual_ID"), file = PLOT12) 
