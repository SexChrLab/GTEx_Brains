# Inputs
# mash.lfsr: matrix with regions as columns, genes as rows, and mashr LFSRs as values
# mash.beta: matrix with regions as columns, genes as rows, and mashr betas as values
# region.levels: vector with regions in preferred order
# region.colors: vector with color choices in hex RGB.
library(mashr)
library(ggplot2)
DIR <- '/scratch/mjpete11/human_monkey_brain/mashr/threshold_plots/'

source('/scratch/mjpete11/human_monkey_brain/mashr/kennys_example/_include_options.R')
mash_results = readRDS('/scratch/mjpete11/human_monkey_brain/mashr/output/mashr_results.rds')
mash.beta = get_pm(mash_results)
mash.lfsr = get_lfsr(mash_results)

threshold.range = seq(0.01,0.2,0.01)

# Loop through a series of cutoffs
dataset1 = do.call(rbind,lapply(threshold.range,function(threshold) {
	# raw counts
	# unique counts
	# Create matrices on only those genes that are significant in that direction in only one region
	this.inc.lfsr = mash.lfsr[rowSums(mash.lfsr < threshold & mash.beta > 0) == 1,]
	this.inc.beta = mash.beta[rowSums(mash.lfsr < threshold & mash.beta > 0) == 1,]
	this.dec.lfsr = mash.lfsr[rowSums(mash.lfsr < threshold & mash.beta < 0) == 1,]
	this.dec.beta = mash.beta[rowSums(mash.lfsr < threshold & mash.beta < 0) == 1,]
	out = rbind(
		# Counts/proportions of total significant genes with positive betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold & mash.beta > 0)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold & mash.beta > 0) /
					sum(mash.lfsr < threshold & mash.beta > 0)
			),
			direction = 'increase',
			type = 'total'
		),
		# Counts/proportions of total significant genes with negative betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(mash.lfsr < threshold & mash.beta < 0)),
			proportion = as.numeric(
				colSums(mash.lfsr < threshold & mash.beta < 0) /
					sum(mash.lfsr < threshold & mash.beta < 0)
			),
			direction = 'decrease',
			type = 'total'
		),
		# Counts/proportions of unique significant genes with positive betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold & this.inc.beta > 0)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold & this.inc.beta > 0) /
					sum(this.inc.lfsr < threshold & this.inc.beta > 0)
			),
			direction = 'increase',
			type = 'unique'
		),
		# Counts/proportions of unique significant genes with negative betas in each region
		data.frame(
			threshold,
			region = colnames(mash.lfsr),
			count = as.integer(colSums(this.inc.lfsr < threshold & this.inc.beta < 0)),
			proportion = as.numeric(
				colSums(this.inc.lfsr < threshold & this.inc.beta < 0) /
					sum(this.inc.lfsr < threshold & this.inc.beta < 0)
			),
			direction = 'decrease',
			type = 'unique'
		)
	)
	out$region = factor(out$region,levels=region.levels)
	out
}))
# The y-axis below can be switched between proportion and count below. The rank order will stay the same
# Plot for total number of significant genes per region
p = ggplot(subset(dataset1,type=='total'),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic() +
	theme(legend.title=element_blank()) +
	xlab('LFSR threshold') +
	ylab('Count')
pdf(paste0(DIR, "plot1.pdf"))
p
dev.off()

# Plot of unique significant genes per region
p = ggplot(subset(dataset1,type=='unique'),aes(threshold,count,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic() +
	theme(legend.title=element_blank()) +
	xlab('LFSR threshold') +
	ylab('Count')
pdf(paste0(DIR, "plot2.pdf"))
p
dev.off()

# Plot for total number of significant genes per region
p = ggplot(subset(dataset1,type=='total'),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic() +
	theme(legend.title=element_blank()) +
	xlab('LFSR threshold') +
	ylab('Proportion')
pdf(paste0(DIR, "plot3.pdf"))
p
dev.off()

# Plot of unique significant genes per region
p = ggplot(subset(dataset1,type=='unique'),aes(threshold,proportion,color=region)) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
	scale_color_manual(values=region.colors) +
	theme_classic() +
	theme(legend.title=element_blank()) +
	xlab('LFSR threshold') +
	ylab('Proportion')
pdf(paste0(DIR, "plot4.pdf"))
p
dev.off()

# Number of genes in n regions
dataset2 = do.call(rbind,lapply(threshold.range,function(threshold) {
	# raw counts
	table(rowSums(mash.lfsr < threshold & mash.beta > 0))
	table(rowSums(mash.lfsr < threshold & mash.beta < 0))
	out = rbind(
		data.frame(
			threshold,
			n.regions = as.integer(names(table(rowSums(mash.lfsr < threshold & mash.beta > 0)))),
			count = as.integer(table(rowSums(mash.lfsr < threshold & mash.beta > 0))),
			proportion = as.numeric(
				table(rowSums(mash.lfsr < threshold & mash.beta > 0)) / nrow(mash.lfsr)
			),
			direction = 'increase'
		),
		data.frame(
			threshold,
			n.regions = as.integer(names(table(rowSums(mash.lfsr < threshold & mash.beta < 0)))),
			count = as.integer(table(rowSums(mash.lfsr < threshold & mash.beta < 0))),
			proportion = as.numeric(
				table(rowSums(mash.lfsr < threshold & mash.beta < 0)) / nrow(mash.lfsr)
			),
			direction = 'decrease'
		)
	)
	out
}))
p = ggplot(droplevels(subset(dataset2,n.regions>0)),aes(threshold,count,color=factor(n.regions))) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
#	scale_color_manual(values=c('#000000',region.colors)) +
	scale_color_manual(values=region.colors) +
	theme_classic() +
	theme(legend.title=element_blank()) +
	xlab('LFSR threshold') +
	ylab('Count')
pdf(paste0(DIR, "plot5.pdf"))
p
dev.off()
p = ggplot(droplevels(subset(dataset2,n.regions>0)),aes(threshold,proportion,color=factor(n.regions))) +
	geom_line() +
	facet_wrap(~direction,nrow=2) +
	scale_x_continuous(trans='reverse') +
#	scale_color_manual(values=c('#000000',region.colors)) +
	scale_color_manual(values=region.colors) +
	theme_classic() +
	theme(legend.title=element_blank()) +
	xlab('LFSR threshold') +
	ylab('Proportion')
pdf(paste0(DIR, "plot6.pdf"))
p
dev.off()
