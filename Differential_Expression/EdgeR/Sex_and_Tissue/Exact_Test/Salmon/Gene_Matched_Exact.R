# This script looks at differential gene expression between males and females within each brain tissue type.
# gene, matched
setwd("/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon")

METADATA <- snakemake@input[[1]] 
COUNTS <- snakemake@input[[2]]
MD_PLOT <- snakemake@output[[1]]
VOLCANO_PLOT <- snakemake@output[[2]]
UP_JSON <- snakemake@output[[3]]
DOWN_JSON <- snakemake@output[[4]]

# Load packages                                                                 
library(tximport)                                                               
library(GenomicFeatures)                                                        
library(edgeR) 
library(readr)
library(stringr)
library(gridExtra)
library(grid)
library(rjson)
library(dplyr)
library(org.Hs.eg.db)

# Read Metadata CSV.                                                            
Samples <- read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample   

# Metadata split into list of dfs by tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Read in count matrix
cts <- read.table(COUNTS, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Remove Samples not present in metadata
cts <- cts[intersect(colnames(cts), as.character(Samples$Sample))]

# Function to split count matrix into list of dfs
Tissues <- names(Meta)

Split_Cols <- function(w, z){
  names <- which(colnames(w) %in% Samples[['Sample']] & Samples[['Tissue']] == z)
  res <- cts[, names]
  return(res)
}

# Split count matrix into list of dfs
res <- list()
for (i in Tissues){
  res[[i]] <- Split_Cols(w=cts, z=i)
}

# Check if columns were subset correctly
check <- list()
for (i in 1:13){
  check[[i]] <- colnames(res[[i]]) %in% subset(Samples$Sample, Samples$Tissue==Tissues[[i]])
}

# Create design matrix: Step 1
# Create sex and tissue factor
Factor_Func <- function(x){
  res <- factor(paste(x$Tissue, x$Sex, sep="."))
  return(res)
}
Groups <- lapply(Meta, Factor_Func)

# Add column with new factor
Col_Bind <- function(df, fc){
  cbind(df, fc)
}
Meta <- Map(Col_Bind, Meta, Groups)

# Sort cols in cts in same order as rows in Meta
Sort_Cols <- function(x, z){
  x <- x[, match(rownames(z), colnames(x))]
  return(x)
}
cts <- Map(Sort_Cols, x=res, z=Meta)

# Create list of DGEList objects
DGE_Func <- function(df, lst){
  res <- DGEList(df, group=lst)
  return(res)
}
y <- Map(DGE_Func, cts, Groups)

# Create design matrix: Step 2; Part 1
Model_Func <- function(fc, df){
  res <- model.matrix(~0 + fc, data = df$samples)
  return(res)
}
Design <- Map(Model_Func, Groups, y)

# Remove "fc" from colnames
Rename_Cols <- function(x){
  colnames(x) <- gsub("fc", "", colnames(x))
  return(x)
}
Design <- lapply(Design, Rename_Cols)

# Create design matrix: Step 2; Part 2
Set_Levels <- function(x, z){
  colnames(x) <- levels(z$samples$group)
  return(x)
}
Design <- Map(Set_Levels, x=Design, z=y)

# Keep only genes expressed in at least half the samples
Keep <- lapply(y, function(x){
  rowSums(cpm(x)>1)>=11
})

Filter_Func <- function(x, k){
  x <- x[k, , keep.lib.sizes=FALSE]
}
y <- Map(Filter_Func, y, Keep)

# TMM Normalization
y <- lapply(y, calcNormFactors)

# Estimate common dispersion and tagwise dispersions in one run (recommended)
Dispersion_Func <- function(a, b){
  estimateDisp(a, b, robust=TRUE)
}
y <- Map(Dispersion_Func, y, Design)

#---------------------------------------------------------------------------------------------------------------------
# Test for DGX with Exact Test
#---------------------------------------------------------------------------------------------------------------------
# Comaprisons to test
# Set females as baseline
Pairs <- list(c("Amygdala.Female", "Amygdala.Male"), 
              c("Anterior.Female", "Anterior.Male"), 
              c("Caudate.Female", "Caudate.Male"),
              c("Cerebellar.Female", "Cerebellar.Male"),
              c("Cerebellum.Female", "Cerebellum.Male"),
              c("Cortex.Female", "Cortex.Male"),
              c("Frontal_Cortex.Female", "Frontal_Cortex.Male"),
              c("Hippocampus.Female", "Hippocampus.Male"), 
              c("Hypothalamus.Female", "Hypothalamus.Male"), 
              c("Nucleus_Accumbens.Female", "Nucleus_Accumbens.Male"), 
              c("Putamen.Female", "Putamen.Male"),  
              c("Spinal_Cord.Female", "Spinal_Cord.Male"), 
              c("Substantia_Nigra.Female", "Substantia_Nigra.Male"))

# Apply exact test
Exact_Func <- function(x, comp){
  exactTest(x, comp)
}
Exact_Res <- Map(Exact_Func, y, Pairs)

#---------------------------------------------------------------------------------------------------------------------
# Summarize results 
#---------------------------------------------------------------------------------------------------------------------
# Report sig DGX genes
Up_Reg <- function(x){
  res <- topTags(x, n=Inf, p=0.05)$table
  up <- res[res$logFC > 0, ]
  return(up)
}

Down_Reg <- function(x){
  res <- topTags(x, n=Inf, p=0.05)$table
  down <- res[res$logFC < 0, ]
  return(down)
}
Up_Top <- lapply(Exact_Res, Up_Reg)
Down_Top <- lapply(Exact_Res, Down_Reg)

# Make table of up and down regulated genes for each tissue
Get_Vec <- function(x){
  res <- rownames(x)
  return(res)
}
Up_Genes.V <- lapply(Up_Top, Get_Vec)
Down_Genes.V <- lapply(Down_Top, Get_Vec)

# Write to file
Up_Json <- toJSON(Up_Genes)
Down_Json <- toJSON(Down_Genes)

#write(Up_Json, UP_JSON)
#write(Down_Json, DOWN_JSON)

# Get summary of results as table
Summary_Func <- function(x){
  res <- summary(decideTests(x))
  return(res)
}
Results_df <- lapply(Exact_Res, Summary_Func)

#---------------------------------------------------------------------------------------------------------------------
# Mean-Difference Plots
#---------------------------------------------------------------------------------------------------------------------
# Non-replicates only
# Titles and counts without replicates
ExactRes_No_Reps <- Exact_Res[-c(4,6)]
Tissues_No_Reps <- Tissues[-c(4,6)]

# Plot Mean-Difference  plots on one page
MD_Plot_Func <- function(x, w){
  plotMD(x, main=w, legend=FALSE, hl.col=c("green", "blue"), cex=1.2)
  mtext('Salmon: Gene Mean-Difference Plots; Exact Test', side = 3, outer = TRUE, cex=1.2, line=3)
  mtext('Average log CPM', side = 1, outer = TRUE, line=1)
  mtext('Log-fold-change', side = 2, outer = TRUE, line=2)
}

# Write to file
pdf(MD_PLOT)
par(mfrow = c(3, 4), cex=0.4, mar = c(3, 3, 3, 2), oma =c(6, 6, 6, 2), xpd=TRUE)  # margins: c(bottom, left, top, right)
Res_Plots <- Map(MD_Plot_Func, x=ExactRes_No_Reps, w=Tissues_No_Reps)
#legend(15.0, 5.0, legend=c("Upregulated in females","Not Signigicant", "Downregulated in females"), pch = 16, 
#       col = c("green","black", "blue"), bty = "o", xpd=NA, cex=2.0)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Volcano Plots
#---------------------------------------------------------------------------------------------------------------------
# Make df of values for axis
Volcano_Func <- function(x){
  cbind(x$table$logFC, -log10(x$table[,"PValue"]))
}
Volcano_Res <- lapply(Exact_Res, Volcano_Func)

# Coerce to df
Volcano_Res <- lapply(Volcano_Res, as.data.frame)

# Rename columns
colnames <- c("logFC", "negLogPval")

Rename_Cols_Func <- function(x){
  setNames(x, colnames)
}
Volcano_Res <- lapply(Volcano_Res, Rename_Cols_Func)

# Plot
Plot_Func <- function(a, b, c, d){
  plot(a, pch=19, main=b, xlab = '', ylab = '', las = 1)
  with(inner_join(a, c), points(logFC, negLogPval, pch=19, col="green"))
  with(inner_join(a, d), points(logFC, negLogPval, pch=19, col="blue"))
  abline(a=-log10(0.05), b=0, col="blue") 
  abline(v=c(2,-2), col="red")
  mtext('Salmon: Gene Volcano Plots; Exact Test', side = 3, outer = TRUE,  cex=1.2, line=3)
  mtext('logFC', side = 1, outer = TRUE,  cex=0.8, line=1)
  mtext('negLogPval', side = 2, outer = TRUE, line=2)
}
pdf(VOLCANO_PLOT)
par(mfrow = c(3, 5), cex=0.4, mar = c(2, 2, 4, 2), oma =c(6, 6, 6, 2), xpd=FALSE)
Map(Plot_Func, a=Volcano_Res, b=Tissues, c=Up_Top, d=Down_Top)
legend(25.0, 8.0, inset=0, legend=c("Positive Significant", "Negative Significant", "Not significant"),
       pch=16, cex=2.0, col=c("green", "blue", "black"), xpd=NA)
dev.off()

#---------------------------------------------------------------------------------------------------------------------
# Gene ontology and pathway enrichment analysis
#---------------------------------------------------------------------------------------------------------------------
# Both use the NCBI RefSeq annotation.
# Convert ensemble annotation to NCBI RefSeq
library(biomaRt)

# Remove numbers after decimal in gene IDs
# Those just indicate the version number
Drop_Dot <- function(x){
  res <- gsub("\\..*", "", x)
  return(res)
}
Up_Genes <- lapply(Up_Genes, Drop_Dot)
Down_Genes <- lapply(Down_Genes, Drop_Dot)

# Map enseble IDs to entrez IDs
# e.g.) EntrezIDs <- mapIds(org.Hs.eg.db, keys = Up_Genes[[1]], keytype = "ENSEMBL", column="ENTREZID")
Map_IDs <- function(x){
  res <- mapIds(org.Hs.eg.db, keys = x, keytype = "ENSEMBL", column="ENTREZID")
  return(res)
}
Up_EntrezIDs <- lapply(Up_Genes[-c(5,8)], Map_IDs) # Drop cerebellum/hippocampus since no DEGs called
Down_EntrezIDs <- lapply(Down_Genes, Map_IDs) 

#------------------------------------------------------------------------------------------------------------
# Map ensembl IDs to gene names without transcript version number
Genes <- mapIds(org.Hs.eg.db, keys = Up_Genes[[1]], keytype = "ENSEMBL", column="ENTREZID")

ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id",
  values = Up_Genes[[1]],
  mart = ensembl)
res

# map ensembl IDs to gene names with transcript version number
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
res.1 <- getBM(
  attributes = c("ensembl_gene_id_version", "ensembl_gene_id", "external_gene_name", "description"),
  filters = "ensembl_gene_id_version",
  values = Up_Genes.V[[1]],
  mart = ensembl)
res.1

#------------------------------------------------------------------------------------------------------------

# Some transcripts do not return an Entrez ID
# These either  do not have an 'ungapped mapping of this gene onto the GRCh38 assembly' according to ensembl website.
# Drop NAs for now
Drop_NA <- function(x){
  res <- x[!is.na(x)]
  return(res)
}
Up_EntrezIDs <- lapply(Up_EntrezIDs, Drop_NA)
Down_EntrezIDs <- lapply(Down_EntrezIDs, Drop_NA)

# GO with entrez IDs
# NOTE: This will be an incomplete list of DEGs since I dropped the NAs!
Gene_Ont <- function(x){
  res <- goana(x, species="Hs") # dropped rownames()
  return(res)
}
Up_GO <- lapply(Up_EntrezIDs, Gene_Ont)
Down_GO <- lapply(Down_EntrezIDs, Gene_Ont)

# Rank by differential expression p-val
# Ontology options: 'MF': molecular function, 'BP': biological process, 'CC': celular component
TOP_GO <- function(x, w){
  res <- topGO(x, ontology=c('BP'), sort=rownames(w), number=Inf)
  return(res)
}
Up_GO.Res <- Map(TOP_GO, x=Up_GO, w=Up_Genes[-c(5,8)])
Down_GO.Res <- Map(TOP_GO, x=Down_GO, w=Down_Genes)

# Remove any DE p vals >0.05
Filter_PVal <- function(x){
  res <- x[x[,'P.DE']<0.05,]
}
Up_GO.Filtered <- lapply(Up_GO_Res, Filter_PVal)
Down_GO.Filtered <- lapply(Down_GO_Res, Filter_PVal)

# Bind list of dfs into single df 
All.Up_GO <- do.call(rbind, Up_GO.Filtered)
All.Down_GO <- do.call(rbind, Down_GO.Filtered)

# Find any pathways that are enriched in all brain regions with DEGs
# To find any duplicate rows
All_Enriched <- All.Up_GO[All.Up_GO$Term %in% All.Up_GO$Term[duplicated(All.Up_GO$Term)],]

# To find rows duplicated a certain number of times
All.Up_GO[All.Up_GO$Term %in% names(which(table(All.Up_GO$Term) > 11)), ] # none present in all 11 tissues with DEGs
All.Down_GO[All.Down_GO$Term %in% names(which(table(All.Down_GO$Term) > 12)), ]

# X chm dosage compensation is the only pathway present in all tissues with UpReg DEGs
All.Up_GO[All.Up_GO$Term %in% names(which(table(All.Up_GO$Term) > 7)), ] 


# Look at overlap in enrichment between brain regions
# Try with UpSetR
# library(UpSetR)
# test <- upset(data=All.Up_GO, nsets = 11, nintersects = NA, sets = NULL,
#               keep.order = F, set.metadata = NULL, intersections = NULL,
#               matrix.color = "gray23", main.bar.color = "gray23",
#               mainbar.y.label = "Intersection Size", mainbar.y.max = NULL,
#               sets.bar.color = "gray23", sets.x.label = "Set Size",
#               point.size = 2.2, line.size = 0.7, mb.ratio = c(0.7, 0.3),
#               expression = NULL, att.pos = NULL, att.color = main.bar.color,
#               order.by = c("freq", "degree"), decreasing = c(T, F),
#               show.numbers = "yes", number.angles = 0, group.by = "degree",
#               cutoff = NULL, queries = NULL, query.legend = "none",
#               shade.color = "gray88", shade.alpha = 0.25, matrix.dot.alpha = 0.5,
#               empty.intersections = NULL, color.pal = 1, boxplot.summary = NULL,
#               attribute.plots = NULL, scale.intersections = "identity",
#               scale.sets = "identity", text.scale = 1, set_size.angles = 0,
#               set_size.show = FALSE, set_size.numbers_size = NULL,
#               set_size.scale_max = NULL)
# 
# movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=TRUE, sep=";" )
# View(upset(movies, nsets = 7, nintersects = 30, mb.ratio = c(0.5, 0.5),
#       order.by = c("freq", "degree"), decreasing = c(TRUE,FALSE)))


# Write to table
# write.csv(All.Up_GO, file='Up_GO.csv')
# write.csv(All.Down_GO, file='Down_GO.csv')

# KEGG
Kegg_Path <- function(x){
  res <- kegga(x, species="Hs")
  return(res)
}
KEGG <- lapply(Exact_Res, Kegg_Path)

TOP_KEGG <- function(x, w){
  res <- topKEGG(x, sort=rownames(w), number=10)
  return(res)
}
Up_Kegg_Res <- Map(TOP_KEGG, x=KEGG, w=Up_Genes)
Down_Kegg_Res <- Map(TOP_KEGG, x=KEGG, w=Down_Genes)

#---------------------------------------------------------------------------------------------------------------------
# To figure out how edgeR decides DEGs, try subsetting the test results where:
# logFC > 2, PValue < 0.05, AND FDR < 0.05
Filtered.Amygdala <- Exact_Res[[1]][(Exact_Res[[1]]$)]

df[!(df$gender == "woman" & df$age > 40 & df$bp = "high"), ]

