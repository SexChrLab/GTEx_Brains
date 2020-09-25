#!/usr/bin/env Rscript

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(mashr)

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

apply(mash_lfsr, 2, function(x) sum(x < fsr_cutoff))
 
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
table(apply(mash_lfsr, 1, function(x) sum(x < fsr_cutoff)))

# Genes that are significant in just one region
apply(mash_lfsr[apply(mash_lfsr, 1, function(x) sum(x < fsr_cutoff)) == 1,], 2, function(x) sum(x < fsr_cutoff))

#_____________________________________________________________________________
# Plots 
#_____________________________________________________________________________
# Put together longform data frame of mashr results
b_mash = data.frame(
	expand.grid(rownames(mash_beta),colnames(mash_beta)),
	qval=as.numeric(mash_lfsr),
	beta=as.numeric(mash_beta),
	stringsAsFactors=FALSE
)

p = ggplot(subset(b_mash,qval < fsr_cutoff),aes(beta,color=Var2)) +
	scale_color_manual(name='Region',values=region.colors) +
	geom_density() +
	theme_classic(base_size=5) +
	guides(color = guide_legend(ncol = 2)) +
	xlab(expression(italic(beta))) +
	ylab('Density') 
ggsave(p,file=paste0('/scratch/mjpete11/human_monkey_brain/mashr/output/plot1.pdf'),width=7,height=3,useDingbats=FALSE)

p = ggplot(melt(tapply(b_mash$qval,b_mash$Var2,function(x) sum(x < fsr.cutoff,na.rm=TRUE))),aes(Var1,value,fill=Var1)) +
	geom_bar(stat='identity') +
	scale_fill_manual(name='Region',values=region.colors) +
	theme_classic(base_size=5) +
	theme(legend.position='none',axis.text.x=element_text(angle=-45,hjust=0,vjust=1)) +
	xlab('Regions') +
	ylab('Number of genes')
ggsave(p,file=paste0('/scratch/mjpete11/human_monkey_brain/mashr/output/plot2.pdf'),width=7,height=3,useDingbats=FALSE)


#_____________________________________________________________________________
# Plot for poster 
#_____________________________________________________________________________
# Metadata counting of significant effects
b_mash_split_genes = split(b_mash,b_mash$Var1)

# Tally up sig region-combos, as well as up-regulated and down-regulated
region_combinations = table(unlist(lapply(b_mash_split_genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta != 0)$Var2,collapse='-')
})))

# Get rid of blanks
region_combinations = region_combinations[as.logical(nchar(names(region_combinations)))]

# Sort in same order (by total)
region_combinations = region_combinations[order(region_combinations,decreasing=TRUE)]

region_combinations_inc = region_combinations_dec = integer(length(region_combinations))
names(region_combinations_inc) = names(region_combinations_dec) = names(region_combinations)

# Ensure all have the same order
region_combinations_inc[names(region_combinations)] = table(unlist(lapply(b_mash_split_genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta > 0)$Var2,collapse='-')
})))[names(region_combinations)]

region_combinations_dec[names(region_combinations)] = table(unlist(lapply(b_mash_split_genes,function(x) {
	paste(subset(x,qval < fsr.cutoff & beta < 0)$Var2,collapse='-')
})))[names(region_combinations)]

# Set NAs to 0
region_combinations_inc[is.na(region_combinations_inc)] = 0
region_combinations_dec[is.na(region_combinations_dec)] = 0

# The sum of up-regulated genes and down-regulated genes do not add up to the sum of significant genes
# This is because genes can have different directions, causing their regions to be tallied differently.
# Thus, report based on the sum of up- and down-regulated region combos
region_combination_sum = region_combinations_inc + region_combinations_dec

region_combination_sum = sort(region_combination_sum,decreasing=TRUE)
region_combinations = region_combinations[names(region_combination_sum)]
region_combinations_inc = region_combinations_inc[names(region_combination_sum)]
region_combinations_dec = region_combinations_dec[names(region_combination_sum)]

# Set NAs to 0
region_combinations_inc[is.na(region_combinations_inc)] = 0
region_combinations_dec[is.na(region_combinations_dec)] = 0

# plot parameters
fraction_shared_cutoff <- 0.85
region_levels <- colnames(mash_lfsr)

width_of_bars = 0.8
region_combinations_results = do.call(rbind,lapply(1:length(region_combinations),function(i) {
	x = names(region_combinations)[i]
	n = unlist(lapply(strsplit(names(region_combinations),'-'),length))[i]
	count_all = as.integer(region_combinations[x])
	count_inc = as.integer(region_combinations_inc[x])
	count_dec = as.integer(region_combinations_dec[x])
	out = integer(length(keep_genes))
	names(out) = c(names(keep_genes))
	out[unlist(strsplit(x,split='-'))] = 1
	out = rbind(
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction_shared_cutoff * length(keep_genes),
			region=factor(region_levels,levels=region_levels),
			region_sig = factor(ifelse(as.logical(out),names(out),NA),levels=region_levels),
			value = 1,
			xmin = seq(1,length(keep_genes)) - (width_of_bars)/2,
			xmax = seq(1,length(keep_genes)) + (width_of_bars)/2,
			ymin = i - (width_of_bars)/2,
			ymax = i + (width_of_bars)/2,
			chart = 'meta'
		),
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction_shared_cutoff * length(keep_genes),
			region=NA,
			region_sig=NA,
			value=count_all,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_all'
		),
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction_shared_cutoff * length(keep_genes),
			region=NA,
			region_sig=NA,
			value=count_inc,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_increase'
		),
		data.frame(
			combination=i,
			n_regions = n,
			share_region = n > fraction_shared_cutoff * length(keep_genes),
			region=NA,
			region_sig=NA,
			value=count_dec,
			xmin = NA,
			xmax = NA,
			ymin = NA,
			ymax = NA,
			chart='count_decrease'
		)
	)
	out
}))

max.plot = length(keep_genes) * 2
base.color = '#000000'
fade.color = '#000000'
highlight.color = '#ff0000'
base.size = 14

meta.combinations.results.plot = subset(region_combinations_results,combination < max.plot)

p1 = ggplot(subset(meta.combinations.results.plot,chart=='meta')) +
	geom_blank(aes(x = region,y=combination)) +
	geom_rect(
		aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=region,alpha=region_sig,color=region_sig),
		linetype=1
	) +
	scale_fill_manual(values=region.colors) +
	scale_alpha_manual(values=rep(1,length(region.colors)),na.value=0) +
	scale_color_manual(values=rep('black',length(region.colors)),na.value=0) +
	scale_y_continuous(trans='reverse',expand=c(0,0)) +
	scale_x_discrete(position='top',expand=c(0,0)) +
	coord_equal() +
	theme_classic(base_size=base.size) + 
	theme(
		legend.position='none',
		axis.line=element_blank(),
		axis.title=element_blank(),
		axis.text.x=element_text(angle=45,vjust=0,hjust=0,margin=margin(t=-1)),
		axis.text.y=element_blank(),
		axis.ticks=element_blank()
	)
p2 = ggplot(subset(meta.combinations.results.plot,chart=='count_increase')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width_of_bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(expand=c(0,0),limits=c(0,ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100)) +
	coord_flip() +
	ylab('Upregulated') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
p3 = ggplot(subset(meta.combinations.results.plot,chart=='count_decrease')) +
	geom_bar(aes(combination,value),stat='identity',fill=base.color,width=width_of_bars) +
	scale_x_continuous(trans='reverse',expand=c(0,0)) +
	scale_y_continuous(trans='reverse',expand=c(0,0),limits=c(ceiling(max(subset(meta.combinations.results.plot,chart %in% c('count_increase','count_decrease'))$value,na.rm=TRUE)/100)*100,0)) +
	coord_flip() +
	ylab('Downregulated') +
	theme_classic(base_size=base.size) +
	theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())

library(egg)

pdf(file='/scratch/mjpete11/human_monkey_brain/mashr/output/plot3.pdf',useDingbats=FALSE,height=7,width=11)
	ggarrange(p3,p1,p2,ncol=3,nrow=1,widths=c(1,1,1),newpage=FALSE)
dev.off()
