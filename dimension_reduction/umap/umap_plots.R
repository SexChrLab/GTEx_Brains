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
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
SEX_PROJ <- file.path(BASE, "dimension_reduction/umap/sex_chr/sex_projection.csv")
NO_SEX_PROJ <- file.path(BASE, "dimension_reduction/umap/no_sex/no_sex_projection.csv")
 
# Output
PLOT1 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_tissue.pdf")
PLOT2 <- file.path(BASE, "dimension_reduction/umap/no_sex/umap_tissue.pdf")
PLOT3 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_sex.pdf")
PLOT4 <- file.path(BASE, "dimension_reduction/umap/no_sex/umap_sex.pdf")
PLOT5 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_age.pdf")
PLOT6 <- file.path(BASE, "dimension_reduction/umap/no_sex/umap_age.pdf")
PLOT7 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_isc.pdf")
PLOT8 <- file.path(BASE, "dimension_reduction/umap/no_sex/umap_isc.pdf")
PLOT9 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_rin.pdf")
PLOT10<- file.path(BASE, "dimension_reduction/umap/no_sex/umap_rin.pdf")
PLOT11 <- file.path(BASE, "dimension_reduction/umap/sex_chr/umap_indiv.pdf")
PLOT12 <- file.path(BASE, "dimension_reduction/umap/no_sex/umap_indiv.pdf")

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
