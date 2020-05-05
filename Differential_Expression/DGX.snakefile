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
GENE = "Age_Matched/Gene/"
GENE = "Matched/Gene/"
TABLE = "Test_Res/"

rule all:
    input:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Gene/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Matched/Gene/Substantia_Nigra.csv"

rule AgeMatched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, EXACT, SALMON, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule Matched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, EXACT, SALMON, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, SALMON, GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule AgeMatched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, FTEST, SALMON, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule Matched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, FTEST, SALMON, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, SALMON, GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule AgeMatched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, RATIO, SALMON, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Ratio.R"

rule Matched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, RATIO, SALMON, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, SALMON, GENE, "Substantia_Nigra.csv")
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
        os.path.join(BASE, EXACT, HISAT, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Substantia_Nigra.csv")
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
        os.path.join(BASE, EXACT, HISAT, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Upreg.json"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Downreg.json"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Anterior.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Caudate.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Putamen.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, HISAT, GENE, "Substantia_Nigra.csv")
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
        os.path.join(BASE, FTEST, HISAT, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Substantia_Nigra.csv")
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
        os.path.join(BASE, FTEST, HISAT, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Upreg.json"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Downreg.json"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Anterior.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Caudate.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Putamen.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, HISAT, GENE, "Substantia_Nigra.csv")
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
        os.path.join(BASE, RATIO, HISAT, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Substantia_Nigra.csv")
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
        os.path.join(BASE, RATIO, HISAT, GENE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Upreg.json"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Downreg.json"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Anterior.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Caudate.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Putamen.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, HISAT, GENE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Ratio.R"
