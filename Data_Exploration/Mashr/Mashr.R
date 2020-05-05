library(mashr)
library(readr)
library(data.table)

# Read in tables with logFC and make the effect matrix
# Cols: gene ID followed by tissue type, Rows: gene, Values: logFC
# Exact/Salmon/Age/Gene, Exact/Age/Salmon/Trans, Exact/Salmon/Match/Gene, Exact/Salmon/Match/Trans
# Paths to csv files
p1 <- "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene"
# Read in csv files as list of dfs (they are actually tab delim...)
ESAG <- list.files(path = p1, pattern = "*.csv", full.names=T)
ESAG
res <- Map(fread, header=T, ESAG)

# Name each df in list
names(res) <- c("Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex",
                "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", 
                "Putamen", "Spinal_Cord", "Substantia_Nigra")

# Name columns then drop last two by name
Name_Drop <- function(x){
    x <- data.frame(x)
    colnames(x) <- c("gene_ID", "logFC", "logCPM", "PValue")
    drp <- c("logFC", "PValue")
    x <- x[,!colnames(x) %in% drp]
    return(x)
}    
tmp <- lapply(res, Name_Drop)

# Left join list of dfs by gene_ID col; if not present in a df, fill with zero
test <- Reduce(function(df1, df2) merge(df1, df2, by="gene_ID", all=TRUE), tmp)
colnames(test)[2:ncol(test)] <- names(res)



# Create the mashr data object
mash.age = mash_set_data(simdata$Bhat, simdata$Shat)

# Compute canonical covariance matrices
U.c.age = cov_canonical(mash.age)  
