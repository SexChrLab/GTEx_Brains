#!/usr/bin/Rscript

# Purpose: Summarize mashr results
library(mashr)
library(dplyr)

# Input
source('/scratch/mjpete11/human_monkey_brain/mashr/kennys_example/_include_options.R')
mash_results = readRDS('/scratch/mjpete11/human_monkey_brain/mashr/output/mashr_results.rds')
mash_beta = get_pm(mash_results)
mash_sbet = get_pm(mash_results) / get_psd(mash_results)
mash_lfsr = get_lfsr(mash_results)

# Output
TABLE1 <- '/scratch/mjpete11/human_monkey_brain/mashr/output/total_sig.csv'
KEEP_GENES <- '/scratch/mjpete11/human_monkey_brain/mashr/output/keep_genes.csv'
NREGION <- '/scratch/mjpete11/human_monkey_brain/mashr/output/sig_n_regions.csv'
LONELY_SIG <- '/scratch/mjpete11/human_monkey_brain/mashr/output/sig_1_region.csv'
BIAS <- '/scratch/mjpete11/human_monkey_brain/mashr/output/bias_per_region.csv'
TOTAL_BIAS <- '/scratch/mjpete11/human_monkey_brain/mashr/output/total_bias.csv'
REGION_PROP <- '/scratch/mjpete11/human_monkey_brain/mashr/output/region_prop.csv'

#_____________________________________________________________________________ 
# Summary of absolute number of sDEGs
#_____________________________________________________________________________ 
# Significance threshold 
fsr_cutoff <- 0.05

# Table: total number of significant genes in each region at fsr 0.05
tab1 <- apply(mash_lfsr, 2, function(x) sum(x < fsr_cutoff))
tab1 <- stack(tab1)[,2:1]
colnames(tab1) <- c('region', 'sig_genes_0.05')

# Does a higher thresold result in more genes called sig across regions?
tab2 <- apply(mash_lfsr, 2, function(x) sum(x < 0.1))
tab2 <- stack(tab2)[,2:1]
colnames(tab2) <- c('region', 'sig_genes_0.1')

tab3 <- apply(mash_lfsr, 2, function(x) sum(x < 0.01))
tab3 <- stack(tab3)[,2:1]
colnames(tab3) <- c('region', 'sig_genes_0.01')

# combine into one table
tab0 <- merge(tab2, tab1, by = 'region')
tab0 <- merge(tab0, tab3, by = 'region')
write.csv(tab0, TABLE1)

# Are the total number of sDEGs at a couple thresholds subsets of each other?
tmp1 <- which(apply(mash_lfsr, 2, function(x) x < 0.1))
tmp2 <- which(apply(mash_lfsr, 2, function(x) x < 0.05))
tmp3 <- which(apply(mash_lfsr, 2, function(x) x < 0.01))
length(tmp1); length(tmp2); length(tmp3) # 8544, 6480, 3628
length(intersect(tmp1, tmp2)) # 6,480
length(intersect(tmp2, tmp3)) # 3,628

# Table: genes that passed filtering in each region
res <- apply(mash_lfsr, 2, function(x) x < fsr_cutoff)
str(res)
tmp <- split(res, rep(1:ncol(res), each = nrow(res)))
# Get row index of sig genes per region
tmp <- lapply(tmp, function(x) which(x))
# Get list of genes that pass filtering in each region
keep_genes <- lapply(tmp,function(i){row.names(mash_lfsr[i,])})
names(keep_genes) <- colnames(mash_lfsr)
str(keep_genes)

# Write to file list of genes that past filtering in each region
# First, reshape object (list of vectors of different lengths) into df
n_obs <- sapply(keep_genes, length)
seq_max <- seq_len(max(n_obs))
mat <- sapply(keep_genes, "[", i=seq_max)
dim(mat)
mat[1:5,1:5]
write.csv(mat, KEEP_GENES)

# How many genes are significant in n regions
nregion <- table(apply(mash_lfsr, 1, function(x) sum(x < fsr_cutoff)))
nregion <- as.data.frame(nregion)
colnames(nregion) <- c('regions', 'sig_genes')
write.csv(nregion, NREGION, row.names=FALSE)

# Genes that are significant in just one region
lonely_sig <- apply(mash_lfsr[apply(mash_lfsr, 1, function(x) sum(x < fsr_cutoff)) == 1,], 2, function(x) sum(x < fsr_cutoff))
lonely_sig <- setNames(stack(lonely_sig)[2:1], c('region', 'sig_genes'))
write.csv(lonely_sig, LONELY_SIG)

# How many genes are in the union?
# Combine list of lists  into one list and count the number of unique genes
lst <- unlist(keep_genes, recursive = FALSE)
length(unique(lst)) # 1,224

#_____________________________________________________________________________ 
# Summary of direction of sDEGs 
#_____________________________________________________________________________ 
# Code upregulated genes as 1 and downregulated genes as -1 
tissues <- colnames(mash_lfsr)
genes <- rownames(mash_lfsr)
out <- do.call(cbind, lapply(tissues, function(i) {
	tmp_tissue <- do.call("c", lapply(genes, function(j) {
		if (mash_lfsr[j, i] < fsr_cutoff & mash_beta[j, i] > 0) {
			1
		} else if (mash_lfsr[j, i] < fsr_cutoff & mash_beta[j, i] < 0) {
			-1
		} else {
		0
}
}))
data.frame(i = tmp_tissue)
}))
rownames(out) <- genes
colnames(out) <- tissues

# Table summarizing sDEGs upset plot
# Count the number of 1s and -1s
res <- as.data.frame(t(apply(out, 2, table)))
res <- cbind(rownames(res), res)
rownames(res) <- NULL
colnames(res) <- c('region', 'female_upreg', 'no_diff', 'male_upreg')
write.csv(res, BIAS)

# Is the proportion of sex biased genes the same across regions or is it
# region specific?
prop <- res %>%
		rowwise() %>%
		mutate(fem_prop = female_upreg/(female_upreg + no_diff + male_upreg)) %>%
		mutate(male_prop = male_upreg/(female_upreg + no_diff + male_upreg)) %>%
		mutate(fem_ratio = female_upreg/male_upreg) %>%
		mutate(male_ratio = male_upreg/female_upreg) %>%
		mutate_if(is.numeric, ~round(., 3))
write.csv(prop, REGION_PROP)

# Is there more bias in one sex compared to the other?
# Which genes are female biased across all genes?
which(apply(out, 1, function(x) sum(x < 0) == 11))

# Which genes are male biased across all genes?
which(apply(out, 1, function(x) sum(x > 0) == 11))

# Number of genes upregulated in females across n regions
fem <- as.numeric(table(apply(out, 1, function(x) sum(x < 0))))

# Number of genes upregulated in males across n regions
male <- as.numeric(table(apply(out, 1, function(x) sum(x > 0))))

# Combine into one df
total_bias <- data.frame(n_tissues = 0:11, female_upreg = fem, male_upreg = male)
total_bias <- total_bias %>%
		rowwise() %>%
		mutate(fem_prop = female_upreg/nrow(mash_lfsr)) %>% # total genes tested
		mutate(male_prop = male_upreg/nrow(mash_lfsr)) %>% # 13,468
		mutate(fem_ratio = female_upreg/male_upreg) %>%
		mutate(male_ratio = male_upreg/female_upreg) %>%
		mutate_if(is.numeric, ~round(., 3))
write.csv(total_bias, TOTAL_BIAS, row.names = FALSE)
