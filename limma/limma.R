#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Differential gene expression with limma-voom 
#------------------------------------------------------------------------------

#Load packages                                                                 
library(limma)
library(edgeR) 
library(readr)
library(stringr)

# Constants
BASE <- "/scratch/mjpete11/human_monkey_brain"
METADATA <- file.path(BASE, "data/metadata/metadata.csv")
COUNTS <- file.path(BASE, "data/expression_matrices/output/filtered_counts_with_sex_chr.csv")

# Read in files
meta <- read.csv(METADATA, header = TRUE, stringsAsFactors=FALSE)
counts <- read.csv(COUNTS, header = TRUE, stringsAsFactors = FALSE)

#------------------------------------------------------------------------------
# Split counts by region  
#------------------------------------------------------------------------------
# Set rownames of metadata object equal to sample names.                        
rownames(meta) <- meta$Sample   

# Set tissues column as factor to add levels attribute
meta$Tissue <- as.factor(meta$Tissue)

# Replace . to - in colnames
colnames(counts) <- str_replace_all(colnames(counts),pattern = "\\.","-")

# Metadata split into list of dfs by tissue
meta_lst <- list()
for(i in 1:length(levels(meta$Tissue))){
      meta_lst[[i]] <- meta[meta$Tissue == levels(meta$Tissue)[i],]
}
names(meta_lst) <- levels(meta$Tissue)

# Function to split count matrix into list of dfs
tissues <- names(meta_lst)

Split_Cols <- function(w, z){
  res <- which(colnames(w) %in% subset(meta[['Sample_ID']], meta[['Tissue']] == z))
  counts_lst <- counts[, res]
  return(counts_lst)
}

# Split count matrix into list of dfs
counts_lst <- list()
for (i in tissues){
      counts_lst[[i]] <- Split_Cols(w=counts, z=i)
}

#------------------------------------------------------------------------------
# Sanity check; make sure samples are in same order in count dfs and metadata dfs
#------------------------------------------------------------------------------
Check_Order <- function(x, z){
    counts_lst <- identical(colnames(x), rownames(z))
    return(counts_lst)
}
all(as.logical(Map(Check_Order, x=counts_lst, z=meta_lst))) # TRUE

#------------------------------------------------------------------------------
# Prep list of vectors for design matrices 
#------------------------------------------------------------------------------
# distribution of RIN values
summary(meta$RIN)
quantile(meta$RIN, seq(0,1,.1))

# distribution of ischemic time values
summary(meta$Ischemic_Time)
quantile(meta$Ischemic_Time, seq(0,1,.1))

# Encode covariates as binary values; first factor categoricical vars
sex_groups <- lapply(meta_lst, function(x) factor(x[['Sex']]))

# Function to replace age decade intervals with median of interval
median_age <- function(x){
   x[['Age']][x[['Age']] >= 20 & x[['Age']] < 29] <- 25
   x[['Age']][x[['Age']] >= 30 & x[['Age']] < 39] <- 35
   x[['Age']][x[['Age']] >= 40 & x[['Age']] < 49] <- 45
   x[['Age']][x[['Age']] >= 50 & x[['Age']] < 59] <- 55
   x[['Age']][x[['Age']] >= 60 & x[['Age']] < 69] <- 65
   x[['Age']][x[['Age']] >= 70 & x[['Age']] < 79] <- 75
   x[['Age']][x[['Age']] >= 80 & x[['Age']] < 89] <- 85
   x[['Age']][x[['Age']] >= 90 & x[['Age']] < 99] <- 95
   return(x)
}
meta_lst <- Map(median_age, meta_lst)
age_groups <- lapply(meta_lst, function(x) as.numeric(x[['Age']]))

# Store RIN/ischemic time values as numeric vectors per region
rin_groups <- lapply(meta_lst, function(x) as.numeric(x[['RIN']]))
isc_groups <- lapply(meta_lst, function(x) as.numeric(x[['Ischemic_Time']]))

#------------------------------------------------------------------------------
# TEST: Try limma voom on one region
#------------------------------------------------------------------------------
#tmp.dge <- DGEList(counts_lst[[1]], group = sex_groups[[1]])
#head(tmp.dge)
#tmp.mat <- model.matrix(~group + age_groups[[1]] + sex_groups[[1]]:age_groups[[1]] + isc_groups[[1]] + rin_groups[[1]], data = tmp.dge[['samples']])
#head(tmp.mat)
#write.csv(tmp.mat, "/scratch/mjpete11/human_monkey_brain/limma/design_matrix_amygdala.csv")
#tmp.v <- voom(tmp.dge, tmp.mat)
#head(tmp.v)
#tmp.fit <- lmFit(tmp.v, tmp.mat)
#head(tmp.fit)
#tmp.fit <- eBayes(tmp.fit)
#head(tmp.fit)

#------------------------------------------------------------------------------
# Function to 1) makde DGEList objects, deign matric, then limma-voom per region
#------------------------------------------------------------------------------
dge_by_region <- function(dat, sex, age, isc, rin) {
	dge <- DGEList(dat, group = sex)
	design <- model.matrix(~group + age + sex:age + isc + rin, data = dge[['samples']])
	v <- voom(dge, design)
	fit <- lmFit(v, design)
	fit <- eBayes(fit) # compute t-stat, f-stat, and log-odds of DE
}

# Apply to each region
dge_res <- Map(dge_by_region, dat = counts_lst, sex = sex_groups, 
               age = age_groups, isc = isc_groups, rin = rin_groups)
exists('dge_res')

#------------------------------------------------------------------------------
# Write csv of genes with sig p-vals in at least one region
#------------------------------------------------------------------------------
# Join the resulting p-vals for the male covariate from each region into one df
p_vals <- do.call(cbind, lapply(dge_res, function(x) x[['p.value']][, 2]))
colnames(p_vals) <- levels(factor(meta$Tissue))

# Apply DFR correction column-wise, then drop p-vals with > 0.05
p_vals <- apply(p_vals, 1, p.adjust, n = ncol(p_vals))
dim(p_vals)
summary(t(p_vals))
p_vals[, 1:10]

# Append gene names
p_vals <- cbind(counts[, 1:2], t(p_vals))
head(p_vals)
summary(p_vals)

# Keep only genes with at least one sig value in one region 
rows_to_subset <- which(apply(p_vals[, -c(1:2)], 1, function(x) sum(x <= 0.05) >= 1))
p_vals <- p_vals[rows_to_subset, ]
nrow(p_vals)
rownames(p_vals) <- NULL 

# Write genes with sig p-vals to file
write.csv(p_vals, "/scratch/mjpete11/human_monkey_brain/limma/tables/male_beta_pvals.csv")

#------------------------------------------------------------------------------
# Write the male coefficient and the square root of the standard error
#------------------------------------------------------------------------------
# Combine the male betas (second col in coeff matrix) with the standed dev 
make_tables <- function(x){
    counts_lst <- data.frame(beta=x[['coefficients']][,2], st_dev=x[['stdev.unscaled']][, 2])
    return(counts_lst)
}
table_lst <- lapply(dge_res, make_tables)

# Append the gene names to the output
table_lst <- lapply(table_lst, function(x) cbind(gene_ID = counts[, 1], x))

# Write tables 
sapply(names(table_lst), function(x) write.csv(table_lst[[x]],
       file = paste0("/scratch/mjpete11/human_monkey_brain/limma/tables/", x, ".csv"), 
       row.names = FALSE))

