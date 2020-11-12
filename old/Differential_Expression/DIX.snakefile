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
AGE = "Age_Matched/Transcript/"
MATCH = "Matched/Transcript/"
TABLE = "Test_Res/"

rule all:
    input:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Age_Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Age_Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Age_Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Age_Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Transcript/Substantia_Nigra.csv",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Matched/Transcript/Substantia_Nigra.csv"

rule AgeMatched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, EXACT, SALMON, AGE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Upreg.json"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Downreg.json"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Anterior.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Caudate.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Putamen.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, SALMON, AGE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule Matched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, EXACT, SALMON, MATCH, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Volcano.pdf"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Upreg.json"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Downreg.json"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Amygdala.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Anterior.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Caudate.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Putamen.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, SALMON, MATCH, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule AgeMatched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, FTEST, SALMON, AGE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Upreg.json"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Downreg.json"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Anterior.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Caudate.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Putamen.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, SALMON, AGE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule Matched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, FTEST, SALMON, MATCH, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Volcano.pdf"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Upreg.json"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Downreg.json"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Amygdala.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Anterior.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Caudate.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Putamen.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, SALMON, MATCH, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule AgeMatched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, RATIO, SALMON, AGE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Upreg.json"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Downreg.json"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Anterior.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Caudate.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Putamen.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, SALMON, AGE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Ratio.R"

rule Matched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(BASE, RATIO, SALMON, MATCH, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Volcano.pdf"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Upreg.json"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Downreg.json"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Amygdala.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Anterior.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Caudate.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Putamen.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, SALMON, MATCH, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Ratio.R"

rule AgeMatched_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Transcript_Matrix.csv")
    output:
        os.path.join(BASE, EXACT, HISAT, AGE, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Volcano.pdf"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Upreg.json"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Downreg.json"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Amygdala.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Anterior.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Caudate.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Putamen.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, HISAT, AGE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Exact.R"

rule Matched_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Transcript_Matrix.csv")
    output:
        os.path.join(BASE, EXACT, HISAT, MATCH, "MD_Plot.pdf"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Volcano.pdf"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Upreg.json"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Downreg.json"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Amygdala.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Anterior.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Caudate.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Cerebellar.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Cerebellum.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Frontal_Cortex.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Hippocampus.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Hypothalamus.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Putamen.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Spinal_Cord.csv"),
        os.path.join(BASE, EXACT, HISAT, MATCH, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Exact.R"

rule AgeMatched_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Transcript_Matrix.csv")
    output:
        os.path.join(BASE, FTEST, HISAT, AGE, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Volcano.pdf"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Upreg.json"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Downreg.json"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Amygdala.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Anterior.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Caudate.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Putamen.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, HISAT, AGE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/FTest.R"

rule Matched_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Transcript_Matrix.csv")
    output:
        os.path.join(BASE, FTEST, HISAT, MATCH, "MD_Plot.pdf"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Volcano.pdf"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Upreg.json"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Downreg.json"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Amygdala.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Anterior.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Caudate.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Cerebellar.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Cerebellum.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Frontal_Cortex.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Hippocampus.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Hypothalamus.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Putamen.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Spinal_Cord.csv"),
        os.path.join(BASE, FTEST, HISAT, MATCH, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/FTest.R"

rule AgeMatched_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Transcript_Matrix.csv")
    output:
        os.path.join(BASE, RATIO, HISAT, AGE, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Volcano.pdf"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Upreg.json"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Downreg.json"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Amygdala.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Anterior.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Caudate.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Putamen.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, HISAT, AGE, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Ratio.R"

rule Matched_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(HISAT_DIR, "Amygdala_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Anterior_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Caudate_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellar_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cerebellum_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Frontal_Cortex_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hippocampus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Hypothalamus_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Nucleus_Accumbens_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Putamen_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Spinal_Cord_Transcript_Matrix.csv"),
        os.path.join(HISAT_DIR, "Substantia_Nigra_Transcript_Matrix.csv")
    output:
        os.path.join(BASE, RATIO, HISAT, MATCH, "MD_Plot.pdf"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Volcano.pdf"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Upreg.json"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Downreg.json"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Amygdala.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Anterior.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Caudate.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Cerebellar.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Cerebellum.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Frontal_Cortex.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Hippocampus.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Hypothalamus.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Nucleus_Accumbens.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Putamen.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Spinal_Cord.csv"),
        os.path.join(BASE, RATIO, HISAT, MATCH, "Substantia_Nigra.csv")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Ratio.R"
