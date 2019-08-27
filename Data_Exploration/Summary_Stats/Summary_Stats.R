# Generate summary stats
setwd("/scratch/mjpete11/GTEx/Data_Exploration/Summary_Stats")

library(VennDiagram)
library(rjson)

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

# Venn diagram plot function

# example plot
venn.plot <- venn.diagram(
  x = list("Exact" = Hisat_Up$Amygdala$`Exact Test`, "F Test" = Hisat_Up$Amygdala$`F Test`, "Ratio" = Hisat_Up$Amygdala$`Ratio Test`),
  filename = "venn.tiff",
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
  main = "Test plot",
  fontfamily = "sans",
  na = "remove",
  inverted = FALSE)

venn.plot

# plot function
Files_1 <- c("Amygdala.tiff", "Anterior.tiff", "Caudate.tiff", "Cerebellar.tiff", "Cerebellum.tiff", "Cortex.tiff", "Frontal_Cortex.tiff",
           "Hippocampus.tiff", "Hypothalamus.tiff", "Nucleus_Accumbens.tiff", "Putamen.tiff", "Spinal_Cord.tiff", "Substantia_Nigra.tiff")

Venn_3 <- function(x, FILES, TISSUE){
  venn.plot <- venn.diagram(
    x = list("Exact" = x$'Exact Test', "F Test" = x$'F Test', "Ratio" = x$'Ratio Test'),
    filename = FILES,
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

  venn.plot
}

# Three-way venn diagram of test results from same aligner
Map(Venn_3, x=Hisat_Up, FILES=Files_1, TISSUE=names(Hisat_Up))

# Two-way venn diagram of up/down regulated genes by different aligners; same test
# test plot
venn.plot <- venn.diagram(
  x = list("Salmon" = Down_Ratio$Amygdala$Salmon, "Hisat" = Down_Ratio$Amygdala$Hisat),
  filename = "venn.tiff",
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
  main = "Test plot",
  fontfamily = "sans",
  na = "remove",
  inverted = FALSE)

venn.plot

# function
## TO-DO: These should also include tissue; e.g. 'Amygdala_Up_Exact'
Files_2 <- c("Up_Exact", "Up_Ratio", "Up_FTest", "Down_Exact", "Down_Ratio", "Down_FTest")

Venn_2 <- function(x, FILES, TISSUE){
  venn.plot <- venn.diagram(
    x = list("Salmon" = x$Salmon, "Hisat" = x$Hisat),
    filename = FILES,
    scaled = TRUE,
    col = "transparent",
    fill = c("firebrick", "deepskyblue4"),
    main.pos = c(0.5, 1.0),
    cex = 1.5,
    cat.cex = 1.5,
    main.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(170,-170), # first=F test, middle=Ratio, last=Exact
    cat.dist = c(0.05,0.05),
    cat.fontfamily = "sans",
    main = TISSUE,
    fontfamily = "sans",
    na = "remove",
    inverted = FALSE)
  
  venn.plot
}

# Two-way venn diagram of up/down reg genes for same test type
Map(Venn_3, x=Up_Exact, FILES=Files_2, TISSUE=names(Up_Exact))

