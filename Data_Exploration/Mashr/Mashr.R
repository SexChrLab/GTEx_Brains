library(mashr)
library(readr)
library(data.table)

# Read in tables with logFC and make the effect matrix
# Cols: gene ID followed by tissue type, Rows: gene, Values: logFC
# Exact/Salmon/Age/Gene, Exact/Age/Salmon/Trans, Exact/Salmon/Match/Gene, Exact/Salmon/Match/Trans
# Paths to csv files
# Salmon, Exact
p1 <- "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene"

# Read in csv files as list of dfs (they are actually tab delim...);,pass list of file names with list.files
EASG.lst <- Map(fread, header=FALSE, list.files(path = p1, pattern = "*.csv", full.names=T))

# Assign names to dfs in list
names(EASG.lst) <- c("Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex",
                 "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", 
                 "Putamen", "Spinal_Cord", "Substantia_Nigra")

# Name items in list, name the df columns then drop last two cols by name
Name_Drop <- function(x){
    x <- data.frame(x)
    colnames(x) <- c("gene_ID", "logFC", "logCPM", "PValue")
    drp <- c("logFC", "PValue")
    x <- x[,!colnames(x) %in% drp]
    return(x)
}    
EASG.lst<- lapply(EASG.lst, Name_Drop)

# 
# Left join list of dfs by gene_ID col; if not present in a df, fill with zero
EASG.df <- Reduce(function(df1, df2) merge(df1, df2, by="gene_ID", all=TRUE), EASG.lst)
rownames(EASG.df) <- EASG.df[["gene_ID"]]
EASG.df[["gene_ID"]] <- NULL
colnames(EASG.df)[1:ncol(EASG.df)] <- names(EASG.lst)
EASG.df[is.na(EASG.df)] <- 0

# Check for missing or non-finite values
any(is.nan(EASG.df$logFC)) # No NANs
all(is.finite(EASG.df$logFC)) # All values are finite

# Add an standard error matrix with ncol, nrow equal to the effect matrix; 
# populate with 1s since the se is the same for all samples
#EASG.se <- data.frame(matrix(1, nrow=nrow(EASG.df), ncol=ncol(EASG.df))) 
#rownames(EASG.se) <- rownames(EASG.df) 
#colnames(EASG.se) <- colnames(EASG.df) 
#EASG.lst <- list(Bhat=EASG.df, Shat=EASG.se)

# Create the mashr object using the pre-made standard error matrix
#mash.EASG <- mash_set_data(EASG.lst$Bhat, EASG.lst$Shat)

# Create mashr object using option to automatically generate the standard error matrix (fill with 1s)
mash.EASG <- mash_set_data(Bhat=EASG.df, Shat=1)

# Create the mashr data object
# Simulate data
set.seed(1)
simdata = simple_sims(500,5,1)
mash.age = mash_set_data(simdata$Bhat, simdata$Shat)

# Compute canonical covariance matrices
U.c.age = cov_canonical(mash.age)  
