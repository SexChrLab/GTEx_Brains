from os.path import join

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/"
SALMON_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/"
HISAT_DIR = "/scratch/mjpete11/GTEx/Tissue_Procesing/Annotated/"
BASE = "/scratch/mjpete11/GTEx/Differential_Expression/"
SALMON = "Salmon/"
HISAT = "Hisat/"
EXACT = "Exact_Test/"
FTEST = "F_Test/"
RATIO = "Ratio_Test/"
AGE_GENE = "Age_Matched/Gene/"
MATCH_GENE = "Matched/Gene/"
TABLE = "Test_Res/"

rule all:
    input:
#        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Age_Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Matched/Gene/Substantia_Nigra.csv",
#       "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Age_Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Age_Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Age_Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Matched/Gene/Substantia_Nigra.csv",
#        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Matched/Gene/Substantia_Nigra.csv"

rule AgeMatched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule Matched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule AgeMatched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule Matched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule AgeMatched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Ratio.R"

rule Matched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Ratio.R"

rule AgeMatched_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Gene_Matrix.csv")
    output:
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Exact.R"

rule Matched_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Gene_Matrix.csv")
    output:
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Exact.R"

rule AgeMatched_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Gene_Matrix.csv")
    output:
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/FTest.R"

rule Matched_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Gene_Matrix.csv")
    output:
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/FTest.R"

rule AgeMatched_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Gene_Matrix.csv")
    output:
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Ratio.R"

rule Matched_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Gene_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Gene_Matrix.csv")
    output:
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH_GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Ratio.R"
