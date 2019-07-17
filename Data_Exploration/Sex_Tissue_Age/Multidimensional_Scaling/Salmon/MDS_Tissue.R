# MDS: distances correspond to leading log-fold-changes between each pair of RNA samples.
# Leading log-fold-change is the average (root-mean-square) of the largest absolute log-foldchanges between each pair of samples.

METADATA = file.path('/scratch/mjpete11/GTEx/Metadata/', 'Metadata.csv')
COUNTS = file.path('/scratch/mjpete11/GTEx/Data_Exploration/Count_Matrices/', 'Count_Matrix.tsv')
PLOT_DIR = '/scratch/mjpete11/GTEx/Data_Exploration/Sex_Tissue_Age/Multidimensional_Scaling/Salmon/MDS_Plots/'
TABLE_DIR =  '/scratch/mjpete11/GTEx/Data_Exploration/Sex_Tissue_Age/Multidimensional_Scaling/Salmon/MDS_Table'
PLOT_FILE = 'Salmon_MDS_Plots.pdf'

# Load packages                                                                 
library(readr)
library(stringr)
library(dplyr)
library(limma)
library(edgeR)

# Read Metadata CSV.                                                            
Samples = read.csv(METADATA, header = TRUE)

# Set rownames of metadata object equal to sample names.                        
rownames(Samples) <- Samples$Sample                                             

# Read in count matrix
cts <- read.table(COUNTS, sep="\t")

# Replace . to - in colnames
colnames(cts) <- str_replace_all(colnames(cts),pattern = "\\.","-")

# Create DGEList object
#y <- DGEList(cts)

# Create every combination of tissue matrices
# Make list of lists of samples for each tissue
Tissue_Lst <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Tissue_Lst[[i]] <- as.vector(Samples[Samples$Tissue == levels(Samples$Tissue)[i], "Sample"])
}
# Rename lists in list as tissue names
names(Tissue_Lst) <- levels(Samples$Tissue)

# List of metadata for each tissue
Meta <- list()
for(i in 1:length(levels(Samples$Tissue))){
  Meta[[i]] <- Samples[Samples$Tissue == levels(Samples$Tissue)[i],]
}
names(Meta) <- levels(Samples$Tissue)

# Split counts into list of dfs by tissue
Tissue_Count <- list()
for(i in 1:length(Tissue_Lst)){
  Tissue_Count[[i]] <- cts[,which(colnames(cts) %in% Tissue_Lst[[i]])]
}
names(Tissue_Count) <- names(Tissue_Lst)

# Check that samples in count data are listed in metadata
all(rownames(Meta$Amygdala) %in% colnames(Tissue_Count$Amygdala))

which(!rownames(Meta$Substantia_Nigra) %in% colnames(Tissue_Count$Substantia_Nigra))

### Temporary ###
Meta$Substantia_Nigra <- Meta$Substantia_Nigra[-c(16),]

Check_Func <- function(a, b){
  Check <- all(rownames(a) %in% colnames(b))
}
Res <- Map(Check_Func, a=Meta, b=Tissue_Count)

# Match order of samples in metadata and count data
Match_Order <- function(a, b){ # Will throw an error if dfs are different lengths
  a <- a[, rownames(b)]
}

Tissue_Count <- Map(Match_Order, a=Tissue_Count, b=Meta)

# Check that order matches
Match_Check <- function(a, b){
  Match <- all(rownames(a) == colnames(b))
}
Res_2 <- Map(Match_Check, a=Meta, b=Tissue_Count)

# Create DGEList object for each tissue count matrix
DGE_lst <- lapply(Tissue_Count, function(x){
  DGEList(x)
}) 

# MDS plots
setwd(PLOT_DIR)
pdf(PLOT_FILE)

# Plot constants 
colors <- c("blue", "darkgreen")
par(mar=c(8, 4.1, 4.1, 8), xpd=TRUE) # margins: bottom, left, top, and right

# MDS plot function for dimensions 1 and 2
MDS_FUN <- function(DGE, NAME, META, TOP, NUM_GENES) {
  k2_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = TOP, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(1,2), 
                    col = colors[META$Sex],
                    main = paste('MDS Plot: Dim 1 and 2; ', NUM_GENES,'; ', NAME))
  legend("topright", inset=c(-0.3,0), legend=levels(META$Sex), pch=16, cex=1,col=colors, xpd=TRUE) # legend will be out of range in viewer
  return(k2_MDS)
}
k2_MDS <- Map(MDS_FUN, DGE = DGE_lst, NAME = names(DGE_lst), META = Meta, TOP = 20, NUM_GENES='Top 20')

# MDS plot function for dimensions 3 and 4
MDS_FUN <- function(DGE, NAME, META, TOP, NUM_GENES) {
  k4_MDS <- plotMDS(DGE,
                    gene.selection = "common",
                    top = TOP, 
                    pch = 16, 
                    cex = 1, 
                    dim.plot = c(3,4), 
                    col = colors[META$Sex],
                    main = paste('MDS Plot: Dim 3 and 4; ', NUM_GENES,'; ', NAME))
  legend("topright", inset=c(-0.3,0), legend=levels(Samples$Sex), pch=16, cex=1,col=colors, xpd=TRUE)
  return(k4_MDS)
}
k4_MDS <- Map(MDS_FUN, DGE = DGE_lst, NAME = names(DGE_lst), META = Meta, TOP = 100, NUM_GENES='Top 100')

dev.off()

# Eigenvector table
setwd(TABLE_DIR)
sapply(names(k2_MDS), 
       function (x) write.table(k2_MDS[[x]]['cmdscale.out'], file=paste(x, "txt", sep=".")))


