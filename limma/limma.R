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
# Make list of DGEList objecounts with corresponding design matrix 
#------------------------------------------------------------------------------
# distribution of RIN values
summary(meta$RIN)
quantile(meta$RIN, seq(0,1,.1))

# distribution of ischemic time values
summary(meta$Ischemic_Time)
quantile(meta$Ischemic_Time, seq(0,1,.1))

# Encode covariates as binary values; first factor categoricical vars
sex_groups <- lapply(meta_lst, function(x) factor(x[['Sex']]))
age_groups <- lapply(meta_lst, function(x) factor(x[['Age']]))

# Bin the RIN and ischemic time values so 20% of samples are in 5 bins.
bin_data <- function(x, bins) {
	thres <- quantile(x, probs = seq(0, 1, 1 / bins))[1:bins]
    sapply(x, function(i) {
				           max(which(i >= thres))
					      })
}
rin_groups <- lapply(meta_lst, function(x) factor(bin_data(x$RIN, bin = 5)))
isc_groups <- lapply(meta_lst, function(x) factor(bin_data(x$Ischemic_Time, bin = 5)))

### TEST
### Try limma voom on one region
#tmp.dge <- DGEList(counts_lst[[1]], group = sex_groups[[1]])
#head(tmp.dge)
#tmp.mat <- model.matrix(~0 + sex_groups[[1]] + age_groups[[1]] + sex_groups[[1]]:age_groups[[1]] + isc_groups[[1]] + rin_groups[[1]], data = tmp.dge[['samples']])
## Skipt RIN and ischemic time for now --> leads to 'coefficients not estimable' errors
##tmp.mat <- model.matrix(~0 + sex_groups[[1]] + age_groups[[1]] + sex_groups[[1]]:age_groups[[1]], data = tmp.dge[['samples']])
##colnames(mat) <- c("female", "male")
#head(tmp.mat)
#tmp.v <- voom(tmp.dge, tmp.mat)
#head(tmp.v)
#tmp.fit <- lmFit(tmp.v, tmp.mat)
#head(tmp.fit)
#tmp.fit <- eBayes(tmp.fit)
#head(tmp.fit)

# Function to 1) makde DGEList objects, deign matric, then limma-voom per region
dge_by_region <- function(dat, sex, age, isc, rin) {
	dge <- DGEList(dat, group = sex)
	design <- model.matrix(~0 + sex + age + sex:age + isc + rin, data = dge[['samples']])
#	design <- model.matrix(~0 + sex + age + sex:age, data = dge[['samples']])
	v <- voom(dge, design)
	fit <- lmFit(v, design)
	fit <- eBayes(fit) # compute t-stat, f-stat, and log-odds of DE
}

# Apply to each region
dge_res <- Map(dge_by_region, dat = counts_lst, sex = sex_groups, 
               age = age_groups, isc = isc_groups, rin = rin_groups)
exists('dge_res')

# Combine the betas (second col in coeff matrix) with the sqrt standard error (var^2)
# Write the female coefficient and the square root of the standard error
#test <- data.frame(beta=dge_res[[1]][['coefficients']][,1], sqrt_SE=sqrt(dge_res[[1]][['s2.post']]))
make_tables <- function(x){
    counts_lst <- data.frame(beta=x[['coefficients']][,2], sqrt_SE=sqrt(x[['s2.post']]))
    return(counts_lst)
}
table_lst <- lapply(dge_res, make_tables)

# Write tables 
sapply(names(table_lst), function(x) write.csv(table_lst[[x]],
       file = paste0("/scratch/mjpete11/human_monkey_brain/limma/tables/", x, ".csv"), 
       row.names = FALSE))

