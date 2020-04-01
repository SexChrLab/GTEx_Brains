#!/usr/bin/env Rscript

setwd("/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats/")

# Output
H_M_U <- 'Test_Res/Hisat/Matched/Up%d.tiff'
H_M_D <- 'Test_Res/Hisat/Matched/Down%d.tiff'
H_A_U <- 'Test_Res/Hisat/Age_Matched/Up%d.tiff'
H_A_D <- 'Test_Res/Hisat/Age_Matched/Down%d.tiff'
S_M_U <- 'Test_Res/Salmon/Matched/Up%d.tiff'
S_M_D <- 'Test_Res/Salmon/Matched/Down%d.tiff'
S_A_U <- 'Test_Res/Salmon/Age_Matched/Up%d.tiff'
S_A_D <- 'Test_Res/Salmon/Age_Matched/Down%d.tiff'
E_M_U <- 'Aligner_Res/Exact/Matched/Up%d.tiff'
E_M_D <- 'Aligner_Res/Exact/Matched/Down%d.tiff'
E_A_U <- 'Aligner_Res/Exact/Age_Matched/Up%d.tiff'
E_A_D <- 'Aligner_Res/Exact/Age_Matched/Down%d.tiff'
F_M_U <- 'Aligner_Res/FTest/Matched/Up%d.tiff'
F_M_D <- 'Aligner_Res/FTest/Matched/Down%d.tiff'
F_A_U <- 'Aligner_Res/FTest/Age_Matched/Up%d.tiff'
F_A_D <- 'Aligner_Res/FTest/Age_Matched/Down%d.tiff'
R_M_U <- 'Aligner_Res/Ratio/Matched/Up%d.tiff'
R_M_D <- 'Aligner_Res/Ratio/Matched/Down%d.tiff'
R_A_U <- 'Aligner_Res/Ratio/Age_Matched/Up%d.tiff'
R_A_D <- 'Aligner_Res/Ratio/Age_Matched/Down%d.tiff'

# Libraries
library(ggplot2)
library(VennDiagram)
library(ggVennDiagram)
library(rjson)
library(grDevices)

#------------------------------------------------------------------------------------------------------------------
# Paths
#------------------------------------------------------------------------------------------------------------------
# Exact Test, Salmon
ESMU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Matched/Gene/Upreg.json')
ESMD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Matched/Gene/Downreg.json')

ESAU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene/Upreg.json')
ESAD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene/Downreg.json') 

# F Test, Salmon
FSMU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Matched/Gene/Upreg.json')
FSMD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Matched/Gene/Downreg.json')

FSAU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Age_Matched/Gene/Upreg.json')
FSAD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Age_Matched/Gene/Downreg.json')

# Ratio Test, Salmon
RSMU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Matched/Gene/Upreg.json')
RSMD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Matched/Gene/Downreg.json')

RSAU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Age_Matched/Gene/Upreg.json')
RSAD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Age_Matched/Gene/Downreg.json')

# Exact Test, Hisat
EHMU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Matched/Gene/Upreg.json')
EHMD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Matched/Gene/Downreg.json')

EHAU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Age_Matched/Gene/Upreg.json')
EHAD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Age_Matched/Gene/Downreg.json')

# F Test, Hisat
FHMU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Matched/Gene/Upreg.json')
FHMD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Matched/Gene/Downreg.json')

FHAU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Age_Matched/Gene/Upreg.json')
FHAD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Age_Matched/Gene/Downreg.json')

# Ratio Test, Hisat
RHMU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Matched/Gene/Upreg.json')
RHMD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Matched/Gene/Downreg.json')

RHAU <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Gene/Upreg.json')
RHAD <- fromJSON(file='/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Gene/Downreg.json')

#------------------------------------------------------------------------------------------------------------------
# Remove version number from gene IDs (numbers after decimal)
#------------------------------------------------------------------------------------------------------------------
Drop_Version <- function(LST){
  res <- lapply(LST, function(x) gsub("\\..*", "", x))
  return(res)              
}

# Pack variables (lists of named vectors) into list, apply function to each item, then unpack variables into the
# global environment.
All_Lsts <- list()
Lst_Names <- c('ESMU', 'ESMD', 'ESAU', 'ESAD', 'FSMU', 'FSMD', 'FSAU', 'FSAD', 'RSMU', 'RSMD', 'RSAU', 'RSAD', 
               'EHMU', 'EHMD', 'EHAU', 'EHAD', 'FHMU', 'FHMD', 'FHAU', 'FHAD', 'RHMU', 'RHMD', 'RHAU', 'RHAD')
All_Lsts[c(Lst_Names)] <- list(ESMU, ESMD, ESAU, ESAD, FSMU, FSMD, FSAU, FSAD, RSMU, RSMD, RSAU, RSAD, 
                  EHMU, EHMD, EHAU, EHAD, FHMU, FHMD, FHAU, FHAD, RHMU, RHMD, RHAU, RHAD)

test <- lapply(All_Lsts, Drop_Version)

list2env(test, .GlobalEnv)

#------------------------------------------------------------------------------------------------------------------
# Make list of test results for each aligner/tissue
# Results are for the 3-way venn diagram of the overlap between genes called as DGX by aligner
#------------------------------------------------------------------------------------------------------------------
Test_Lst <- function(Exact, FTest, Ratio){
  res <- list("Exact Test" = Exact, "FTest" = FTest, "Ratio Test" = Ratio)
}

HMU <- Map(Test_Lst, Exact=EHMU, FTest=FHMU, Ratio=RHMU)
HMD <- Map(Test_Lst, Exact=EHMD, FTest=FHMD, Ratio=RHMD)
HAU <- Map(Test_Lst, Exact=EHMU, FTest=FHMU, Ratio=RHMU)
HAD <- Map(Test_Lst, Exact=EHMD, FTest=FHAD, Ratio=RHMD)

SMU <- Map(Test_Lst, Exact=ESMU, FTest=FSMU, Ratio=RSMU)
SMD <- Map(Test_Lst, Exact=ESMD, FTest=FSMD, Ratio=RSMD)
SAU <- Map(Test_Lst, Exact=ESAU, FTest=FSAU, Ratio=RSAU)
SAD <- Map(Test_Lst, Exact=ESAD, FTest=FSMD, Ratio=RSAD)

#------------------------------------------------------------------------------------------------------------------
# Make list of results by test for each aligner
# Results are for 2-way venn diagram of genes called as DGX by test
#------------------------------------------------------------------------------------------------------------------
Aligner_Lst <- function(H, S){
  res <- list("Salmon" = S, "Hisat" = H)
}

EMU <- Map(Aligner_Lst, H=EHMU, S=ESMU)
EMD <- Map(Aligner_Lst, H=EHMD, S=ESMD)
EAU <- Map(Aligner_Lst, H=EHMU, S=ESAU)
EAD <- Map(Aligner_Lst, H=EHMD, S=ESAD)

FMU <- Map(Aligner_Lst, H=FHMU, S=FSMU)
FMD <- Map(Aligner_Lst, H=FHMD, S=FSMD)
FAU <- Map(Aligner_Lst, H=FHAU, S=FSAU)
FAD <- Map(Aligner_Lst, H=FHAD, S=FSAD)

RMU <- Map(Aligner_Lst, H=RHMU, S=RSMU)
RMD <- Map(Aligner_Lst, H=RHMD, S=RSMD)
RAU <- Map(Aligner_Lst, H=RHAU, S=RSAU)
RAD <- Map(Aligner_Lst, H=RHAD, S=RSAD)

#------------------------------------------------------------------------------------------------------------------
# 3-way Venn diagram plot function
# Make list of results by test for each aligner
# Results are for 2-way venn diagram of genes called as DGX by test
#------------------------------------------------------------------------------------------------------------------
Aligner_Lst <- function(H, S){
  res <- list("Salmon" = S, "Hisat" = H)
}

EMU <- Map(Aligner_Lst, H=EHMU, S=ESMU)
EMD <- Map(Aligner_Lst, H=EHMD, S=ESMD)
EAU <- Map(Aligner_Lst, H=EHMU, S=ESAU)
EAD <- Map(Aligner_Lst, H=EHMD, S=ESAD)

FMU <- Map(Aligner_Lst, H=FHMU, S=FSMU)
FMD <- Map(Aligner_Lst, H=FHMD, S=FSMD)
FAU <- Map(Aligner_Lst, H=FHAU, S=FSAU)
FAD <- Map(Aligner_Lst, H=FHAD, S=FSAD)

RMU <- Map(Aligner_Lst, H=RHMU, S=RSMU)
RMD <- Map(Aligner_Lst, H=RHMD, S=RSMD)
RAU <- Map(Aligner_Lst, H=RHAU, S=RSAU)
RAD <- Map(Aligner_Lst, H=RHAD, S=RSAD)

#------------------------------------------------------------------------------------------------------------------
# 3-way Venn diagram plot function
#------------------------------------------------------------------------------------------------------------------
# Encode strings as numeric factors; apply to nested list to use with ggVennDiagram
Factor_Func <- function(x){
    x <- as.numeric(as.factor(x))
    return(x)
}
HMU <- lapply(HMU, function(x) lapply(x, Factor_Func))
HMD <- lapply(HMD, function(x) lapply(x, Factor_Func))
HAU <- lapply(HAU, function(x) lapply(x, Factor_Func))
HAD <- lapply(HAD, function(x) lapply(x, Factor_Func))
SMU <- lapply(SMU, function(x) lapply(x, Factor_Func))
SMD <- lapply(SMD, function(x) lapply(x, Factor_Func))
SAU <- lapply(SAU, function(x) lapply(x, Factor_Func))
SAD <- lapply(SAD, function(x) lapply(x, Factor_Func))

# Plot function
Test_Venn <- function(LST, TITLE){
    res <- ggVennDiagram(LST[c('Exact Test', 'FTest', 'Ratio Test')]) +
           scale_fill_gradient(low='blue', high='red') +
           ggtitle(TITLE)
    return(res)
}    

#------------------------------------------------------------------------------------------------------------------
# Hisat; By test
#------------------------------------------------------------------------------------------------------------------
# Matched, Up
tiff(H_M_U)
Map(Test_Venn, LST=HMU, TITLE=paste(names(HMU), 'Up regulated DEGs'))
dev.off()

# Matched, Down
tiff(H_M_D)
Map(Test_Venn, LST=HMD, TITLE=paste(names(HMD), 'Down regulated DEGs'))
dev.off()

# Age Matched, Up
tiff(H_A_U)
Map(Test_Venn, LST=HAU, TITLE=paste(names(HAU), 'Up regulated DEGs'))
dev.off()

# Age Matched, Down
tiff(H_A_D)
Map(Test_Venn, LST=HAD, TITLE=paste(names(HAD), 'Down regulated DEGs'))
dev.off()

#------------------------------------------------------------------------------------------------------------------
# Salmon; By test
#------------------------------------------------------------------------------------------------------------------
# Matched, Up
tiff(S_M_U)
Map(Test_Venn, LST=SMU, TITLE=paste(names(SMU), 'Up regulated DEGs'))
dev.off()

# Matched, Down
tiff(S_M_D)
Map(Test_Venn, LST=SMD, TITLE=paste(names(SMD), 'Down regulated DEGs'))
dev.off()

# Age Matched, Up
tiff(S_A_U)
Map(Test_Venn, LST=SAU, TITLE=paste(names(SAU), 'Up regulated DEGs'))
dev.off()

# Age Matched, Down
tiff(S_A_D)
Map(Test_Venn, LST=SAD, TITLE=paste(names(SAD), 'Down regulated DEGs'))
dev.off()

#------------------------------------------------------------------------------------------------------------------
# 2 way venn diagrams
#------------------------------------------------------------------------------------------------------------------
# Encode strings as numeric factors; apply to nested list to use with ggVennDiagram
EMU <- lapply(EMU, function(x) lapply(x, Factor_Func))
EMD <- lapply(EMD, function(x) lapply(x, Factor_Func))
EAU <- lapply(EAU, function(x) lapply(x, Factor_Func))
EAD <- lapply(EAD, function(x) lapply(x, Factor_Func))

SMU <- lapply(SMU, function(x) lapply(x, Factor_Func))
SMD <- lapply(SMD, function(x) lapply(x, Factor_Func))
SAU <- lapply(SAU, function(x) lapply(x, Factor_Func))
SAD <- lapply(SAD, function(x) lapply(x, Factor_Func))

FMU <- lapply(FMU, function(x) lapply(x, Factor_Func))
FMD <- lapply(FMD, function(x) lapply(x, Factor_Func))
FAU <- lapply(FAU, function(x) lapply(x, Factor_Func))
FAD <- lapply(FAD, function(x) lapply(x, Factor_Func))

# Plot function
Aligner_Venn <- function(LST, TITLE){
   res <- ggVennDiagram(LST[c('Salmon', 'Hisat')]) +
          scale_fill_gradient(low='blue', high='red') +
          ggtitle(TITLE)
  return(res)
}    

#------------------------------------------------------------------------------------------------------------------
# Exact Test 
#------------------------------------------------------------------------------------------------------------------
# Matched, Up
tiff(E_M_U)
Map(Aligner_Venn, LST=EMU, TITLE=paste(names(EMU), 'Exact Test, Matched Samples, Up-regulated DEGs'))
dev.off()

# Matched, Down
tiff(E_M_D)
Map(Aligner_Venn, LST=EMD, TITLE=paste(names(EMD), 'Exact Test, Matched Samples, Down-regulated DEGs'))
dev.off()

# Age Matched, Up
tiff(E_A_U)
Map(Aligner_Venn, LST=EAU, TITLE=paste(names(EAU), 'Exact Test, Age Matched Samples, Up-regulated DEGs'))
dev.off()

# Age Matched, Down
tiff(E_A_D)
Map(Aligner_Venn, LST=EAD, TITLE=paste(names(EAD), 'Exact Test, Age Matched Samples, Down-regulated DEGs'))
dev.off()

#------------------------------------------------------------------------------------------------------------------
# F Test 
#------------------------------------------------------------------------------------------------------------------
# Matched, Up
tiff(F_M_U)
Map(Aligner_Venn, LST=FMU, TITLE=paste(names(FMU), 'F Test, Matched Samples, Up-regulated DEGs'))
dev.off()

# Matched, Down
tiff(F_M_D)
Map(Aligner_Venn, LST=FMD, TITLE=paste(names(FMD), 'F Test, Matched Samples, Down-regulated DEGs'))
dev.off()

# Age Matched, Up
tiff(F_A_U)
Map(Aligner_Venn, LST=FAU, TITLE=paste(names(FAU), 'F Test, Age Matched Samples, Up-regulated DEGs'))
dev.off()

# Age Matched, Down
tiff(F_A_D)
Map(Aligner_Venn, LST=FAD, TITLE=paste(names(FAD), 'F Test, Age Matched Samples, Down-regulated DEGs'))
dev.off()

#------------------------------------------------------------------------------------------------------------------
# Ratio Test 
#------------------------------------------------------------------------------------------------------------------
# Matched, Up
tiff(R_M_U)
Map(Aligner_Venn, LST=RMU, TITLE=paste(names(RMU), 'Ratio Test, Matched Samples, Up-regulated DEGs'))
dev.off()

# Matched, Down
tiff(R_M_D)
Map(Aligner_Venn, LST=RMD, TITLE=paste(names(RMD), 'Ratio Test, Matched Samples, Down-regulated DEGs'))
dev.off()

# Age Matched, Up
tiff(R_A_U)
Map(Aligner_Venn, LST=RAU, TITLE=paste(names(RAU), 'Ratio Test, Age Matched Samples, Up-regulated DEGs'))
dev.off()

# Age Matched, Down
tiff(R_A_D)
Map(Aligner_Venn, LST=RAD, TITLE=paste(names(RAD), 'Ratio Test, Age Matched Samples, Down-regulated DEGs'))
dev.off()

