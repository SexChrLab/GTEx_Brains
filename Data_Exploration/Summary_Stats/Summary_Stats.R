# Generate summary stats
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats")

library(VennDiagram)
library(rjson)
library(grDevices)
library(grid)

# Read in json objects
Hisat_Up_Exact <- fromJSON(file='Hisat_Upreg_Exact.json')
Hisat_Down_Exact <- fromJSON(file='Hisat_Downreg_Exact.json')
Hisat_Up_FTest <- fromJSON(file='Hisat_Upreg_FTest.json')
Hisat_Down_FTest <- fromJSON(file='Hisat_Downreg_FTest.json')
Hisat_Up_Ratio <- fromJSON(file='Hisat_Upreg_Ratio.json')
Hisat_Down_Ratio <- fromJSON(file='Hisat_Downreg_Ratio.json')

Salmon_Up_Exact <- fromJSON(file='Salmon_Upreg_Exact.json')
Salmon_Down_Exact <- fromJSON(file='Salmon_Downreg_Exact.json')
Salmon_Up_FTest <- fromJSON(file='Salmon_Upreg_FTest.json')
Salmon_Down_FTest <- fromJSON(file='Salmon_Downreg_FTest.json')
Salmon_Up_Ratio <- fromJSON(file='Salmon_Upreg_Ratio.json')
Salmon_Down_Ratio <- fromJSON(file='Salmon_Downreg_Ratio.json')

# Make list of test results for each aligner/tissue
# results are for the 3-way venn diagram
Test_Lst <- function(E, F, R){
  res <- list("Exact Test" = E, "F Test" = F, "Ratio Test" = R)
}
Hisat_Up <- Map(Test_Lst, E=Hisat_Up_Exact, F=Hisat_Up_FTest, R=Hisat_Up_Ratio)
Hisat_Down <- Map(Test_Lst, E=Hisat_Down_Exact, F=Hisat_Down_FTest, R=Hisat_Down_Ratio)

Salmon_Up <- Map(Test_Lst, E=Salmon_Up_Exact, F=Salmon_Up_FTest, R=Salmon_Up_Ratio)
Salmon_Down <- Map(Test_Lst, E=Salmon_Down_Exact, F=Salmon_Down_FTest, R=Salmon_Down_Ratio)

# Make list of results by test for each aligner
# results are for 2-way venn diagram
Aligner_Lst <- function(H, S){
  res <- list("Salmon" = S, "Hisat" = H)
}
Up_Exact <- Map(Aligner_Lst, H=Hisat_Up_Exact, S=Salmon_Up_Exact)
Down_Exact <- Map(Aligner_Lst, H=Hisat_Down_Exact, S=Salmon_Down_Exact)

Up_FTest <- Map(Aligner_Lst, H=Hisat_Up_FTest, S=Salmon_Up_FTest)
Down_FTest <- Map(Aligner_Lst, H=Hisat_Down_FTest, S=Salmon_Down_FTest)

Up_Ratio <- Map(Aligner_Lst, H=Hisat_Up_Ratio, S=Salmon_Up_Ratio)
Down_Ratio <- Map(Aligner_Lst, H=Hisat_Down_Ratio, S=Salmon_Down_Ratio) 

# 3-way Venn diagram plot function
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
# write to one file
pdf("Hisat_Up.pdf")
Map(Venn_3, x=Hisat_Up, TISSUE=names(Hisat_Up))
dev.off()

pdf("Salmon_Up.pdf")
Map(Venn_3, x=Salmon_Up, TISSUE=names(Salmon_Up))
dev.off()

pdf("Hisat_Down.pdf")
Map(Venn_3, x=Hisat_Down, TISSUE=names(Hisat_Down))
dev.off()

pdf("Salmon_Down.pdf")
Map(Venn_3, x=Salmon_Down, TISSUE=names(Salmon_Down))
dev.off()

# 2 way venn diagrams
Files_UE <- c("Amygdala_Up_Exact", "Anterior_Up_Exact", "Caudate_Up_Exact", "Cerebellar_Up_Exact", "Cerebellum_Up_Exact", "Cortex_Up_Exact",
             "Frontal_Cortex_Up_Exact", "Hippocampus_Up_Exact", "Hypothalamus_Up_Exact", "Nucleus_Accumbens_Up_Exact", "Putamen_Up_Exact",
             "Spinal_Cord_Up_Exact", "Substantia_Nigra_Up_Exact")

Files_DE <- c("Amygdala_Down_Exact", "Anterior_Down_Exact", "Caudate_Down_Exact", "Cerebellar_Down_Exact", "Cerebellum_Down_Exact", 
              "Cortex_Down_Exact", "Frontal_Cortex_Down_Exact", "Hippocampus_Down_Exact", "Hypothalamus_Down_Exact", "Nucleus_Accumbens_Down_Exact",
              "Putamen_Down_Exact", "Spinal_Cord_Down_Exact", "Substantia_Nigra_Down_Exact")

Files_UR <- c("Amygdala_Up_Ratio", "Anterior_Up_Ratio", "Caudate_Up_Ratio","Cerebellar_Up_Ratio", "Cerebellum_Up_Ratio", "Cortex_Up_Ratio",
              "Frontal_Cortex_Up_Ratio", "Hippocampus_Up_Ratio", "Hypothalamus_Up_Ratio", "Nucleus_Accumbens_Up_Ratio", "Putamen_Up_Ratio",
              "Spinal_Cord_Up_Ratio", "Substantia_Nigra_Up_Ratio")

Files_DR <- c("Amygdala_Down_Ratio", "Anterior_Down_Ratio", "Caudate_Down_Ratio","Cerebellar_Down_Ratio", "Cerebellum_Down_Ratio", 
              "Cortex_Down_Ratio", "Frontal_Cortex_Down_Ratio", "Hippocampus_Down_Ratio", "Hypothalamus_Down_Ratio", 
              "Nucleus_Accumbens_Down_Ratio", "Putamen_Down_Ratio", "Spinal_Cord_Down_Ratio", "Substantia_Nigra_Down_Ratio")

Files_UF <- c("Amygdala_Up_FTest", "Anterior_Up_FTest", "Caudate_Up_FTest", "Cerebellar_Up_FTest", "Cerebellum_Up_FTest", "Cortex_Up_FTes",
             "Frontal_Cortex_Up_FTest", "Hippocampus_Up_FTest", "Hypothalamus_Up_FTest", "Nucleus_Accumbens_Up_FTest", "Putamen_Up_FTest",
             "Spinal_Cord_Up_FTest", "Substantia_Nigra_Up_FTest")

Files_DF <- c("Amygdala_Down_FTest", "Anterior_Down_FTest", "Caudate_Down_FTest", "Cerebellar_Down_FTest", "Cerebellum_Down_FTest",
              "Cortex_Down_FTest", "Frontal_Cortex_Down_FTest", "Hippocampus_Down_FTest", "Hypothalamus_Down_FTest",
              "Nucleus_Accumbens_Down_FTest", "Putamen_Down_FTest", "Spinal_Cord_Down_FTest", "Substantia_Nigra_Down_FTest")

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

# write to one file
pdf(file="Up_Exact.pdf")
Map(Venn_2, x=Up_Exact, FILES=Files_UE)
dev.off()

pdf(file="Up_FTest.pdf")
Map(Venn_2, x=Up_FTest, FILES=Files_UF)
dev.off()

pdf(file="Up_Ratio.pdf")
Map(Venn_2, x=Up_Ratio, FILES=Files_UR)
dev.off()

pdf(file="Down_Exact.pdf")
Map(Venn_2, x=Down_Exact, FILES=Files_DE)
dev.off()

pdf(file="Down_FTest.pdf")
Map(Venn_2, x=Down_FTest, FILES=Files_DF)
dev.off()

pdf(file="Down_Ratio.pdf")
Map(Venn_2, x=Down_Ratio, FILES=Files_DR)
dev.off()



