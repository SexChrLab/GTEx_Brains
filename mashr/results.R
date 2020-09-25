#!/usr/bin/Rscript

# Purpose: Summarize mashr results

library(mashr)

# Input
#_____________________________________________________________________________
# Read in files/ generate summary stats
#_____________________________________________________________________________
source('/scratch/mjpete11/human_monkey_brain/mashr/kennys_example/_include_options.R')

mash_results = readRDS('/scratch/mjpete11/human_monkey_brain/mashr/output/mashr_results.rds')
mash_beta = get_pm(mash_results)
mash_sbet = get_pm(mash_results) / get_psd(mash_results)
mash_lfsr = get_lfsr(mash_results)

# Significant genes per region
fsr_cutoff <- 0.2

# Table: total number of significant genes in each region
apply(mash_lfsr, 2, function(x) sum(x < fsr_cutoff))

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
write.csv(mat, '/scratch/mjpete11/human_monkey_brain/mashr/output/keep_genes.csv')

# How many genes are significant in n regions
nregion <- table(apply(mash_lfsr, 1, function(x) sum(x < fsr_cutoff)))
nregion <- as.data.frame(nregion)
colnames(nregion) <- c('regions', 'sig_genes')
write.csv(nregion, '/scratch/mjpete11/human_monkey_brain/mashr/output/sig_n_regions.csv', row.names=FALSE)

# Genes that are significant in just one region
lonely_sig <- apply(mash_lfsr[apply(mash_lfsr, 1, function(x) sum(x < fsr_cutoff)) == 1,], 2, function(x) sum(x < fsr_cutoff))
lonely_sig <- setNames(stack(lonely_sig)[2:1], c('region', 'sig_genes'))
write.csv(lonely_sig, '/scratch/mjpete11/human_monkey_brain/mashr/output/sig_1_region.csv')

# How many genes are in the union?
sum(as.numeric(lapply(keep_genes, function(x) length(x)))) # 11,598

# Is there more bias in one sex compared to the other?
# Is the proportion of sex biased genes the same across regions or is it
# region specific?
# Code upregulated genes as 1 and downregulated genes as -1 
tissues <- colnames(mash_lfsr)
genes <- rownames(mash_lfsr)
out <- do.call(cbind, lapply(tissues, function(i) {
	tmp_tissue <- do.call("c", lapply(genes, function(j) {
		if (mash_lfsr[j, i] < 0.02 & mash_beta[j, i] > 0) {
			1
		} else if (mash_lfsr[j, i] < 0.02 & mash_beta[j, i] < 0) {
			-1
		} else {
		0
}
}))
data.frame(genes = genes, i = tmp_tissue)
}))


