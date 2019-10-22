# Generate summary stats
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats")

library(VennDiagram)
library(rjson)
library(grDevices)
library(grid)

#------------------------------------------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------------------------------------------
# Exact Test, Salmon
S_EXACT_GENE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Matched/Gene/Salmon_Upreg_Exact.json'
S_EXACT_GENE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Matched/Gene/Salmon_Downreg_Exact.json'

S_EXACT_TRANS_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Matched/Transcript/Salmon_Upreg_Exact.json'
S_EXACT_TRANS_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Matched/Transcript/Salmon_Downreg_Exact.json'

S_EXACT_GENE_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Age_Matched/Gene/Salmon_Upreg_Exact.json'
S_EXACT_GENE_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Age_Matched/Gene/Salmon_Downreg_Exact.json'

S_EXACT_TRANS_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Age_Matched/Transcript/Salmon_Upreg_Exact.json'
S_EXACT_TRANS_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Age_Matched/Transcript/Salmon_Downreg_Exact.json'

# F Test, Salmon
S_FTEST_GENE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Matched/Gene/Salmon_Upreg_FTest.json'
S_FTEST_GENE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Matched/Gene/Salmon_Downreg_FTest.json'

S_FTEST_TRANS_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Matched/Transcript/Salmon_Upreg_FTest.json'
S_FTEST_TRANS_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Matched/Transcript/Salmon_Downreg_FTest.json'

S_FTEST_GENE_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Age_Matched/Gene/Salmon_Upreg_FTest.json'
S_FTEST_GENE_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Age_Matched/Gene/Salmon_Downreg_FTest.json'

S_FTEST_TRANS_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Age_Matched/Transcript/Salmon_Upreg_FTest.json'
S_FTEST_TRANS_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Age_Matched/Transcript/Salmon_Downreg_FTest.json'

# Ratio Test, Salmon
S_RATIO_GENE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Matched/Gene/Salmon_Upreg_Ratio.json'
S_RATIO_GENE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Matched/Gene/Salmon_Downreg_Ratio.json'

S_RATIO_TRANS_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Matched/Transcript/Salmon_Upreg_Ratio.json'
S_RATIO_TRANS_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Matched/Transcript/Salmon_Downreg_Ratio.json'

S_RATIO_GENE_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Age_Matched/Gene/Salmon_Upreg_Ratio.json'
S_RATIO_GENE_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Age_Matched/Gene/Salmon_Downreg_Ratio.json'

S_RATIO_TRANS_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Age_Matched/Transcript/Salmon_Upreg_Ratio.json'
S_RATIO_TRANS_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Age_Matched/Transcript/Salmon_Downreg_Ratio.json'

# Exact Test, Hisat
H_EXACT_GENE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Matched/Gene/Hisat_Upreg_Exact.json'
H_EXACT_GENE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Matched/Gene/Hisat_Downreg_Exact.json'

H_EXACT_TRANS_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Matched/Transcript/Hisat_Upreg_Exact.json'
H_EXACT_TRANS_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Matched/Transcript/Hisat_Downreg_Exact.json'

H_EXACT_GENE_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Gene/Hisat_Upreg_Exact.json'
H_EXACT_GENE_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Gene/Hisat_Downreg_Exact.json'

H_EXACT_TRANS_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Transcript/Hisat_Upreg_Exact.json'
H_EXACT_TRANS_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Age_Matched/Transcript/Hisat_Downreg_Exact.json'

# F Test, Hisat
H_FTEST_GENE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Matched/Gene/Hisat_Upreg_FTest.json'
H_FTEST_GENE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Matched/Gene/Hisat_Downreg_FTest.json'

H_FTEST_TRANS_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Matched/Transcript/Hisat_Upreg_FTest.json'
H_FTEST_TRANS_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Matched/Transcript/Hisat_Downreg_FTest.json'

H_FTEST_GENE_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Age_Matched/Gene/Hisat_Upreg_FTest.json'
H_FTEST_GENE_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Age_Matched/Gene/Hisat_Downreg_FTest.json'

H_FTEST_TRANS_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Age_Matched/Transcript/Hisat_Upreg_FTest.json'
H_FTEST_TRANS_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Age_Matched/Transcript/Hisat_Downreg_FTest.json'

# Ratio Test, Hisat
H_RATIO_GENE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Matched/Gene/Hisat_Upreg_Ratio.json'
H_RATIO_GENE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Matched/Gene/Hisat_Downreg_Ratio.json'

H_RATIO_TRANS_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Matched/Transcript/Hisat_Upreg_Ratio.json'
H_RATIO_TRANS_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Matched/Transcript/Hisat_Downreg_Ratio.json'

H_RATIO_GENE_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Age_Matched/Gene/Hisat_Upreg_Ratio.json'
H_RATIO_GENE_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Age_Matched/Gene/Hisat_Downreg_Ratio.json'

H_RATIO_TRANS_AGE_MATCHED_UP <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Age_Matched/Transcript/Hisat_Upreg_Ratio.json'
H_RATIO_TRANS_AGE_MATCHED_DOWN <- '/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Age_Matched/Transcript/Hisat_Downreg_Ratio.json'

# Read in json files
# Salmon, Exact
S_Exact_Gene_Matched_Up <- fromJSON(file=S_EXACT_GENE_MATCHED_UP)
S_Exact_Gene_Matched_Down <- fromJSON(file=S_EXACT_GENE_MATCHED_DOWN)
S_Exact_Gene_Age_Matched_Up <- fromJSON(file=S_EXACT_GENE_AGE_MATCHED_UP)
S_Exact_Gene_Age_Matched_Down <- fromJSON(file=S_EXACT_GENE_AGE_MATCHED_DOWN)

S_Exact_Trans_Matched_Up <- fromJSON(file=S_EXACT_GENE_MATCHED_UP)
S_Exact_Trans_Matched_Down <- fromJSON(file=S_EXACT_GENE_MATCHED_DOWN)
S_Exact_Trans_Age_Matched_Up <- fromJSON(file=S_EXACT_TRANS_AGE_MATCHED_UP) # To Do
S_Exact_Trans_Age_Matched_Down <- fromJSON(file=S_EXACT_TRANS_AGE_MATCHED_DOWN) # To Do

# Salmon, F Test
S_FTest_Gene_Matched_Up <- fromJSON(file=S_FTEST_GENE_MATCHED_UP)
S_FTest_Gene_Matched_Down <- fromJSON(file=S_FTEST_GENE_MATCHED_DOWN)
S_FTest_Gene_Age_Matched_Up <- fromJSON(file=S_FTEST_GENE_AGE_MATCHED_UP)
S_FTest_Gene_Age_Matched_Down <- fromJSON(file=S_FTEST_GENE_AGE_MATCHED_DOWN)

S_FTest_Trans_Matched_Up <- fromJSON(file=S_FTEST_GENE_MATCHED_UP)
S_FTest_Trans_Matched_Down <- fromJSON(file=S_FTEST_GENE_MATCHED_DOWN)
S_FTest_Trans_Age_Matched_Up <- fromJSON(file=S_FTEST_TRANS_AGE_MATCHED_UP)
S_FTest_Trans_Age_Matched_Down <- fromJSON(file=S_FTEST_TRANS_AGE_MATCHED_DOWN)

# Salmon, Ratio
S_Ratio_Gene_Matched_Up <- fromJSON(file=S_RATIO_GENE_MATCHED_UP) 
S_Ratio_Gene_Matched_Down <- fromJSON(file=S_RATIO_GENE_MATCHED_DOWN) 
S_Ratio_Gene_Age_Matched_Up <- fromJSON(file=S_RATIO_GENE_AGE_MATCHED_UP)
S_Ratio_Gene_Age_Matched_Down <- fromJSON(file=S_RATIO_GENE_AGE_MATCHED_DOWN)

S_Ratio_Trans_Matched_Up <- fromJSON(file=S_RATIO_TRANS_MATCHED_UP) 
S_Ratio_Trans_Matched_Down <- fromJSON(file=S_RATIO_TRANS_MATCHED_DOWN) 
S_Ratio_Trans_Age_Matched_Up <- fromJSON(file=S_RATIO_TRANS_AGE_MATCHED_UP)
S_Ratio_Trans_Age_Matched_Down <- fromJSON(file=S_RATIO_TRANS_AGE_MATCHED_DOWN)

# Hisat, Exact
H_Exact_Gene_Matched_Up <- fromJSON(file=H_EXACT_GENE_MATCHED_UP)
H_Exact_Gene_Matched_Down <- fromJSON(file=H_EXACT_GENE_MATCHED_DOWN)
H_Exact_Gene_Age_Matched_Up <- fromJSON(file=H_EXACT_GENE_AGE_MATCHED_UP)
H_Exact_Gene_Age_Matched_Down <- fromJSON(file=H_EXACT_GENE_AGE_MATCHED_DOWN)

H_Exact_Trans_Matched_Up <- fromJSON(file=H_EXACT_TRANS_MATCHED_UP)
H_Exact_Trans_Matched_Down <- fromJSON(file=H_EXACT_TRANS_MATCHED_DOWN)
H_Exact_Trans_Age_Matched_Up <- fromJSON(file=H_EXACT_TRANS_AGE_MATCHED_UP)
H_Exact_Trans_Age_Matched_Down <- fromJSON(file=H_EXACT_TRANS_AGE_MATCHED_DOWN)

# Hisat, F Test
H_FTest_Gene_Matched_Up <- fromJSON(file=H_FTEST_GENE_MATCHED_UP)
H_FTest_Gene_Matched_Down <- fromJSON(file=H_FTEST_GENE_MATCHED_DOWN)
H_FTest_Gene_Age_Matched_Up <- fromJSON(file=H_FTEST_GENE_AGE_MATCHED_UP)
H_FTest_Gene_Age_Matched_Down <- fromJSON(file=H_FTEST_GENE_AGE_MATCHED_DOWN)

H_FTest_Trans_Matched_Up <- fromJSON(file=H_FTEST_GENE_MATCHED_UP)
H_FTest_Trans_Matched_Down <- fromJSON(file=H_FTEST_GENE_MATCHED_DOWN)
H_FTest_Trans_Age_Matched_Up <- fromJSON(file=H_FTEST_TRANS_AGE_MATCHED_UP)
H_FTest_Trans_Age_Matched_Down <- fromJSON(file=H_FTEST_TRANS_AGE_MATCHED_DOWN)

# Hisat, Ratio
H_Ratio_Gene_Matched_Up <- fromJSON(file=H_RATIO_GENE_MATCHED_UP)
H_Ratio_Gene_Matched_Down <- fromJSON(file=H_RATIO_GENE_MATCHED_DOWN)
H_Ratio_Gene_Age_Matched_Up <- fromJSON(file=H_RATIO_GENE_AGE_MATCHED_UP)
H_Ratio_Gene_Age_Matched_Down <- fromJSON(file=H_RATIO_GENE_AGE_MATCHED_DOWN)

H_Ratio_Trans_Matched_Up <- fromJSON(file=H_RATIO_GENE_MATCHED_UP)
H_Ratio_Trans_Matched_Down <- fromJSON(file=H_RATIO_GENE_MATCHED_DOWN)
H_Ratio_Trans_Age_Matched_Up <- fromJSON(file=H_RATIO_TRANS_AGE_MATCHED_UP) # To Do
H_Ratio_Trans_Age_Matched_Down <- fromJSON(file=H_RATIO_TRANS_AGE_MATCHED_DOWN) # To Do

#------------------------------------------------------------------------------------------------------------------
# Make list of test results for each aligner/tissue
# Results are for the 3-way venn diagram of the overlap between genes called as DGX by aligner
#------------------------------------------------------------------------------------------------------------------
Test_Lst <- function(E, F, R){
  res <- list("Exact Test" = E, "F Test" = F, "Ratio Test" = R)
}

# Hisat
# Gene
# Matched
# Up
M.H_Gene_Up <- Map(Test_Lst, E=H_Exact_Gene_Matched_Up, F=H_FTest_Gene_Matched_Up, R=H_Ratio_Gene_Matched_Up)
# Down
M.H_Gene_Down <- Map(Test_Lst, E=H_Exact_Gene_Matched_Down, F=H_FTest_Gene_Matched_Down, R=H_Ratio_Gene_Matched_Down)
# Age Matched
# Up
AM.H_Gene_Up <- Map(Test_Lst, E=H_Exact_Gene_Age_Matched_Up, F=H_FTest_Gene_Age_Matched_Up, R=H_Ratio_Gene_Age_Matched_Up)
# Down
AM.H_Gene_Down <- Map(Test_Lst, E=H_Exact_Gene_Age_Matched_Down, F=H_FTest_Gene_Age_Matched_Down, R=H_Ratio_Gene_Age_Matched_Down)

# By test: Hisat, transcript, matched
# By test: Hisat, transcript, age matched

# Salmon
# Gene
# Matched
# Up
M.S_Gene_Up <- Map(Test_Lst, E=S_Exact_Gene_Matched_Up, F=S_FTest_Gene_Matched_Up, R=S_Ratio_Gene_Matched_Up)
# Down
M.S_Gene_Down <- Map(Test_Lst, E=S_Exact_Gene_Matched_Down, F=S_FTest_Gene_Matched_Down, R=S_Ratio_Gene_Matched_Down)

# Age Matched
# Up
AM.S_Gene_Up <- Map(Test_Lst, E=S_Exact_Gene_Age_Matched_Up, F=S_FTest_Gene_Age_Matched_Up, R=S_Ratio_Gene_Age_Matched_Up)
# Down
AM.S_Gene_Down <- Map(Test_Lst, E=S_Exact_Gene_Age_Matched_Down, F=S_FTest_Gene_Age_Matched_Down, R=S_Ratio_Gene_Age_Matched_Down)

# By test: Salmon, transcript, matched
# By test: Salmon, transcript, age matched

#------------------------------------------------------------------------------------------------------------------
# Make list of results by test for each aligner
# Results are for 2-way venn diagram of genes called as DGX by test
#------------------------------------------------------------------------------------------------------------------
Aligner_Lst <- function(H, S){
  res <- list("Salmon" = S, "Hisat" = H)
}
# By aligner; Genes
# Exact test
# Matched
M.Up_Exact_Gene <- Map(Aligner_Lst, H=H_Exact_Gene_Matched_Up, S=S_Exact_Gene_Matched_Up)
M.Down_Exact_Gene <- Map(Aligner_Lst, H=H_Exact_Gene_Matched_Down, S=S_Exact_Gene_Matched_Down)
# Age matched
AM.Up_Exact_Gene <- Map(Aligner_Lst, H=H_Exact_Gene_Age_Matched_Up, S=S_Exact_Gene_Age_Matched_Up)
AM.Down_Exact_Gene <- Map(Aligner_Lst, H=H_Exact_Gene_Age_Matched_Down, S=S_Exact_Gene_Age_Matched_Down)

# F Test
# Matched
M.Up_FTest_Gene <- Map(Aligner_Lst, H=H_FTest_Gene_Matched_Up, S=S_FTest_Gene_Matched_Up)
M.Down_FTest_Gene <- Map(Aligner_Lst, H=H_FTest_Gene_Matched_Down, S=S_FTest_Gene_Matched_Down)
# Age Matched
AM.Up_FTest_Gene <- Map(Aligner_Lst, H=H_FTest_Gene_Age_Matched_Up, S=S_FTest_Gene_Age_Matched_Up)
AM.Down_FTest_Gene <- Map(Aligner_Lst, H=H_FTest_Gene_Age_Matched_Down, S=S_FTest_Gene_Age_Matched_Down)

# Ratio test
# Matched
M.Up_Ratio_Gene <- Map(Aligner_Lst, H=H_Ratio_Gene_Matched_Up, S=S_Ratio_Gene_Matched_Up)
M.Down_Ratio_Gene <- Map(Aligner_Lst, H=H_Ratio_Gene_Matched_Down, S=S_Ratio_Gene_Matched_Down)
# Age Matched
AM.Up_Ratio_Gene <- Map(Aligner_Lst, H=H_Ratio_Gene_Age_Matched_Up, S=S_Ratio_Gene_Age_Matched_Up)
AM.Down_Ratio_Gene <- Map(Aligner_Lst, H=H_Ratio_Gene_Age_Matched_Down, S=S_Ratio_Gene_Age_Matched_Down)


## Transcript
## Exact Test
## Matched
# M.Up_Exact_Trans <- Map(Aligner_Lst, H=H_Exact_Trans_Matched_Up, S=S_Exact_Trans_Matched_Up)
# M.Down_Exact_Trans <- Map(Aligner_Lst, H=H_Exact_Trans_Matched_Down, S=S_Exact_Trans_Matched_Down)
## Age Matched
# AM.Up_Exact_Trans <- Map(Aligner_Lst, H=H_Exact_Trans_Age_Matched_Up, S=S_Exact_Trans_Age_Matched_Up)
# AM.Down_Exact_Trans <- Map(Aligner_Lst, H=H_Exact_Trans_Age_Matched_Down, S=S_Exact_Trans_Age_Matched_Down)

# # F Test
# # Matched
# M.Up_FTest_Trans <- Map(Aligner_Lst, H=H_FTest_Trans_Matched_Up, S=S_FTest_Trans_Matched_Up)
# M.Down_FTest_Trans <- Map(Aligner_Lst, H=H_FTest_Trans_Matched_Down, S=S_FTest_Trans_Matched_Down)
# # Age Matched
# AM.Up_FTest_Trans <- Map(Aligner_Lst, H=H_FTest_Trans_Age_Matched_Up, S=S_FTest_Trans_Age_Matched_Up)
# AM.Down_FTest_Trans <- Map(Aligner_Lst, H=H_FTest_Trans_Age_Matched_Down, S=S_FTest_Trans_Age_Matched_Down)

# # Ratio Test
# # Matched
# M.Up_Ratio_Trans <- Map(Aligner_Lst, H=H_Ratio_Trans_Matched_Up, S=S_Ratio_Trans_Matched_Up)
# M.Down_Ratio_Trans <- Map(Aligner_Lst, H=H_Ratio_Trans_Matched_Down, S=S_Ratio_Trans_Matched_Down)
# # Age Matched
# AM.Up_Ratio_Trans <- Map(Aligner_Lst, H=H_Ratio_Trans_Age_Matched_Up, S=S_Ratio_Trans_Age_Matched_Up)
# AM.Down_Ratio_Trans <- Map(Aligner_Lst, H=H_Ratio_Trans_Age_Matched_Down, S=S_Ratio_Trans_Age_Matched_Down)

#------------------------------------------------------------------------------------------------------------------
# 3-way Venn diagram plot function
#------------------------------------------------------------------------------------------------------------------
Venn_3 <- function(x, TISSUE){
  venn.plot <- venn.diagram(
    x = list("Exact" = x$'Exact Test', "F Test" = x$'F Test', "Ratio" = x$'Ratio Test'),
    filename = NULL,
    scaled = TRUE,
    col = "transparent",
    fill = c("firebrick", "deepskyblue4", "azure3"),
    main.pos = c(0.5, 1.0),
    cex = 1.5,
    cat.cex = 1.5,
    main.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(-70,60,-90), # first=F test, middle=Ratio, last=Exact
    cat.dist = c(0.05,0.05,0.05),
    cat.fontfamily = "sans",
    main = TISSUE,
    fontfamily = "sans",
    na = "remove",
    inverted = FALSE)
  
  venn.plot <- grid.draw(venn.plot)
  grid.newpage()
  
  return(venn.plot)
}
#------------------------------------------------------------------------------------------------------------------
# Hisat; By test
#------------------------------------------------------------------------------------------------------------------
# Gene
# Matched
# Up
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Hisat_Gene_Up.pdf')
Map(Venn_3, x=M.H_Gene_Up, TISSUE=names(M.H_Gene_Up))
dev.off()

# Down
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Hisat_Gene_Down.pdf')
Map(Venn_3, x=M.H_Gene_Down, TISSUE=names(M.H_Gene_Down))
dev.off()

# Age Matched
# Up
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Hisat_Gene_Up.pdf')
Map(Venn_3, x=AM.H_Gene_Up, TISSUE=names(AM.H_Gene_Up))
dev.off()

# Down
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Hisat_Gene_Down.pdf')
Map(Venn_3, x=AM.H_Gene_Down, TISSUE=names(AM.H_Gene_Down))
dev.off()

# Transcript
# Matched
# Up
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Hisat_Trans_Up.pdf')
# Map(Venn_3, x=M.H_Trans_Up, TISSUE=names(M.H_Trans_Up))
# dev.off()
# 
# # Down
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Hisat_Gene_Trans.pdf')
# Map(Venn_3, x=M.H_Trans_Down, TISSUE=names(M.H_Trans_Down))
# dev.off()
# 
# # Age Matched
# # Up
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Hisat_Trans_Up.pdf')
# Map(Venn_3, x=AM.H_Trans_Up, TISSUE=names(AM.H_Trans_Up))
# dev.off()
# 
# # Down
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Hisat_Trans_Down.pdf')
# Map(Venn_3, x=AM.H_Trans_Down, TISSUE=names(AM.H_Trans_Down))
# dev.off()

#------------------------------------------------------------------------------------------------------------------
# Salmon; By test
#------------------------------------------------------------------------------------------------------------------
# Gene
# Matched
# Up
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Salmon_Gene_Up.pdf')
Map(Venn_3, x=M.S_Gene_Up, TISSUE=names(M.S_Gene_Up))
dev.off()

# Down
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Salmon_Gene_Down.pdf')
Map(Venn_3, x=M.S_Gene_Down, TISSUE=names(M.S_Gene_Down))
dev.off()

# Age Matched
# Up
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Salmon_Gene_Up.pdf')
Map(Venn_3, x=AM.S_Gene_Up, TISSUE=names(AM.S_Gene_Up))
dev.off()

# Down
pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Salmon_Gene_Down.pdf')
Map(Venn_3, x=AM.S_Gene_Down, TISSUE=names(AM.S_Gene_Down))
dev.off()

# # Transcript
# # Matched
# # Up
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Salmon_Trans_Up.pdf')
# Map(Venn_3, x=M.S_Trans_Up, TISSUE=names(M.S_Trans_Up))
# dev.off()
# 
# # Down
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/M_Salmon_Gene_Trans.pdf')
# Map(Venn_3, x=M.S_Trans_Down, TISSUE=names(M.S_Trans_Down))
# dev.off()
# 
# # Age Matched
# # Up
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Salmon_Trans_Up.pdf')
# Map(Venn_3, x=AM.S_Trans_Up, TISSUE=names(AM.S_Trans_Up))
# dev.off()
# 
# # Down
# pdf('/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Test_Res/AM_Salmon_Trans_Down.pdf')
# Map(Venn_3, x=AM.S_Trans_Down, TISSUE=names(AM.S_Trans_Down))
# dev.off()

#------------------------------------------------------------------------------------------------------------------
# 2 way venn diagrams
#------------------------------------------------------------------------------------------------------------------
# Plot names
# Exact
# Up 
UE.Files <- c("Amygdala_Up_Exact", "Anterior_Up_Exact", "Caudate_Up_Exact", "Cerebellar_Up_Exact", "Cerebellum_Up_Exact", "Cortex_Up_Exact",
              "Frontal_Cortex_Up_Exact", "Hippocampus_Up_Exact", "Hypothalamus_Up_Exact", "Nucleus_Accumbens_Up_Exact", "Putamen_Up_Exact",
              "Spinal_Cord_Up_Exact", "Substantia_Nigra_Up_Exact")
# Down 
DE.Files <- c("Amygdala_Down_Exact", "Anterior_Down_Exact", "Caudate_Down_Exact", "Cerebellar_Down_Exact", "Cerebellum_Down_Exact", 
              "Cortex_Down_Exact", "Frontal_Cortex_Down_Exact", "Hippocampus_Down_Exact", "Hypothalamus_Down_Exact", "Nucleus_Accumbens_Down_Exact",
              "Putamen_Down_Exact", "Spinal_Cord_Down_Exact", "Substantia_Nigra_Down_Exact")

# F test
# Up
UF.Files <- c("Amygdala_Up_FTest", "Anterior_Up_FTest", "Caudate_Up_FTest", "Cerebellar_Up_FTest", "Cerebellum_Up_FTest", "Cortex_Up_FTes",
              "Frontal_Cortex_Up_FTest", "Hippocampus_Up_FTest", "Hypothalamus_Up_FTest", "Nucleus_Accumbens_Up_FTest", "Putamen_Up_FTest",
              "Spinal_Cord_Up_FTest", "Substantia_Nigra_Up_FTest")
# Down
DF.Files <- c("Amygdala_Down_FTest", "Anterior_Down_FTest", "Caudate_Down_FTest", "Cerebellar_Down_FTest", "Cerebellum_Down_FTest",
              "Cortex_Down_FTest", "Frontal_Cortex_Down_FTest", "Hippocampus_Down_FTest", "Hypothalamus_Down_FTest",
              "Nucleus_Accumbens_Down_FTest", "Putamen_Down_FTest", "Spinal_Cord_Down_FTest", "Substantia_Nigra_Down_FTest")

# Ratio test
# Up
UR.Files <- c("Amygdala_Up_Ratio", "Anterior_Up_Ratio", "Caudate_Up_Ratio","Cerebellar_Up_Ratio", "Cerebellum_Up_Ratio", "Cortex_Up_Ratio",
              "Frontal_Cortex_Up_Ratio", "Hippocampus_Up_Ratio", "Hypothalamus_Up_Ratio", "Nucleus_Accumbens_Up_Ratio", "Putamen_Up_Ratio",
              "Spinal_Cord_Up_Ratio", "Substantia_Nigra_Up_Ratio")
# Down
DR.Files <- c("Amygdala_Down_Ratio", "Anterior_Down_Ratio", "Caudate_Down_Ratio","Cerebellar_Down_Ratio", "Cerebellum_Down_Ratio", 
              "Cortex_Down_Ratio", "Frontal_Cortex_Down_Ratio", "Hippocampus_Down_Ratio", "Hypothalamus_Down_Ratio", 
              "Nucleus_Accumbens_Down_Ratio", "Putamen_Down_Ratio", "Spinal_Cord_Down_Ratio", "Substantia_Nigra_Down_Ratio")


# Function to plot two-way venn diagrams 
Venn_2 <- function(x, FILES){
  venn.plot <- venn.diagram(
    x = list("Salmon" = x$Salmon, "Hisat" = x$Hisat),
    filename = NULL,
    scaled = TRUE,
    col = "transparent",
    fill = c("firebrick", "deepskyblue4"),
    main.pos = c(0.5, 1.0),
    cex = 1.5,
    cat.cex = 1.5,
    main.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(170,-170), 
    cat.dist = c(0.05,0.05),
    cat.fontfamily = "sans",
    main = FILES,
    fontfamily = "sans",
    na = "remove",
    inverted = FALSE)
  
  venn.plot <- grid.draw(venn.plot)
  grid.newpage()
  
  return(venn.plot)
}

#------------------------------------------------------------------------------------------------------------------
# Hisat 
#------------------------------------------------------------------------------------------------------------------
# Gene
# Matched
# Up
pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/M_Up_Exact.pdf')
Map(Venn_2, x=M.Up_Exact_Gene, FILES=UE.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/M_Up_FTest.pdf')
Map(Venn_2, x=M.Up_FTest_Gene, FILES=UF.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/M_Up_Ratio.pdf')
Map(Venn_2, x=M.Up_Ratio_Gene, FILES=UR.Files)
dev.off()

# Down
pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/M_Down_Exact.pdf')
Map(Venn_2, x=M.Down_Exact_Gene, FILES=DE.Files)
dev.off()

### Error: zero entries for list item 3
pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/M_Down_FTest.pdf')
Map(Venn_2, x=M.Down_FTest_Gene, FILES=DF.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/M_Down_Ratio.pdf')
Map(Venn_2, x=M.Down_Ratio_Gene, FILES=DR.Files)
dev.off()

# Age Matched
# Up
pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/AM_Up_Exact.pdf')
Map(Venn_2, x=AM.Up_Exact_Gene, FILES=UE.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/AM_Up_FTest.pdf')
Map(Venn_2, x=AM.Up_FTest_Gene, FILES=UF.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/AM_Up_Ratio.pdf')
Map(Venn_2, x=AM.Up_Ratio_Gene, FILES=UR.Files)
dev.off()

# Down
pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/AM_Down_Exact.pdf')
Map(Venn_2, x=AM.Down_Exact_Gene, FILES=DE.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/AM_Down_FTest.pdf')
Map(Venn_2, x=AM.Down_FTest_Gene, FILES=DF.Files)
dev.off()

pdf(file='/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/Aligner_Res/AM_Down_Ratio.pdf')
Map(Venn_2, x=AM.Down_Ratio_Gene, FILES=DR.Files)
dev.off()
