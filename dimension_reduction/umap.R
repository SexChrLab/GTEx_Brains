# Uniform manifold approximation and projection for dimension reduction

# load libraries
library(umap)
library(ggplot2)
library(RColorBrewer)

# e.regressed is an expression matrix, with genes as rows and samples (libraries) as columns. For this particular matrix, we used the partial residuals of a model after removing technical covariates. If this doesn't apply to your data, use your voom-normalized expression matrix
# n_neighbors and min_dist are flexible parameters: see https://umap-learn.readthedocs.io/en/latest/parameters.html for explanation
a = umap(t(e.regressed), n_neighbors = 50, min_dist = 0.5)

# extract the projections (a$layout) and join it with metadata (particularly the 'Region' column)
# full.covariates is a vector of model covariate names (in case you want to plot by those as well)
a.umap = data.frame(as.data.frame(a$layout),meta[,c(full.covariates,'Region')])

# Plot UMAP projections and color by Region.
# region.colors is a custom vector of pretty colors.
# batch.variable is a length=1 character vector that contains the name of our batch variable
p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape=batch.variable)) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_batched.pdf',useDingbats=FALSE)

# Redo the plot, this time without highlighting batch
p = ggplot(a.umap,aes_string('V1','V2',color='Region')) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap.pdf',useDingbats=FALSE)

# Redo the plot, this time showing color AND shape to better differentiate points since many colors look similar
p = ggplot(a.umap,aes_string('V1','V2',color='Region',shape='Region')) +
	geom_point(size=1.25) +
	scale_color_manual(values=region.colors) +
	scale_shape_manual(values=region.shapes) +
	theme_classic(base_size=18) +
	xlab('UMAP 1') +
	ylab('UMAP 2') +
	coord_fixed() +
	guides(shape = guide_legend(override.aes = list(size = 3))) +
	theme(axis.ticks=element_blank(),axis.text=element_blank())
ggsave(p,file='figures/data_visualization_dimension_reduction_umap_shapes.pdf',useDingbats=FALSE)
