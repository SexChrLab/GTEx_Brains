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
Test_Lst <- function(E, F, R){
  res <- list("Exact Test" = E, "F Test" = F, "Ratio Test" = R)
}
Hisat_Up <- Map(Test_Lst, E=Hisat_Up_Exact, F=Hisat_Up_FTest, R=Hisat_Up_Ratio)
Hisat_Down <- Map(Test_Lst, E=Hisat_Down_Exact, F=Hisat_Down_FTest, R=Hisat_Down_Ratio)

Salmon_Up <- Map(Test_Lst, E=Salmon_Up_Exact, F=Salmon_Up_FTest, R=Salmon_Up_Ratio)
Salmon_Down <- Map(Test_Lst, E=Salmon_Down_Exact, F=Salmon_Down_FTest, R=Salmon_Down_Ratio)

# Venn diagram plot function





# based on heini's
Venn_Plot <- function(LST, TITLE){
  venn.plot <- draw.triple.venn(
    x = LST,
    filename = "test.tiff",
    scaled = TRUE,
    col = "transparent",
    fill = c("firebrick", "deepskyblue4", "azure3"),
    main.pos = c(0.5, 0.5),
    cex = 1.5,
    cat.cex = 1.5,
    main.cex = 2,
    cat.default.pos = "outer",
    cat.pos = c(-15,15,0),
    cat.dist = c(0.05,0.05,0.05),
    cat.fontfamily = "sans",
    #main = TITLE,
    fontfamily = "sans",
    na = "remove",
    inverted = FALSE)
  
  plot(venn.plot)
}
Map(Venn_Plot, LST=Salmon_Down, TITLE=names(Salmon_Down))




