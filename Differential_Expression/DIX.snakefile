from os.path import join

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/"
SALMON_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/"
HISAT_DIR = "/scratch/mjpete11/GTEx/Tissue_Procesing/Annotated/"
EXACT_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/"
FTEST_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/"
RATIO_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/"
SALMON_AGE_ISO = "Salmon/Age_Matched/Transcript/"
SALMON_MATCH_ISO = "Salmon/Matched/Transcript/"
HISAT_AGE_ISO = "Hisat/Age_Matched/Transcript/"
HISAT_MATCH_ISO = "Hisat/Matched/Transcript/"

rule all:
    input:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Age_Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Age_Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Age_Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Hisat/Age_Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Age_Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Hisat/Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Matched/Transcript/Downreg.json",
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Age_Matched/Transcript/Downreg.json"

rule AgeMatched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_AGE_ISO, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_ISO, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_ISO, "Upreg.json"),
        os.path.join(EXACT_DGX, SALMON_AGE_ISO, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule Matched_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_MATCH_ISO, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_ISO, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_ISO, "Upreg.json"),
        os.path.join(EXACT_DGX, SALMON_MATCH_ISO, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Exact_Test/Salmon/Exact.R"

rule AgeMatched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_AGE_ISO, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_ISO, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_ISO, "Upreg.json"),
        os.path.join(FTEST_DGX, SALMON_AGE_ISO, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule Matched_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_MATCH_ISO, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_ISO, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_ISO, "Upreg.json"),
        os.path.join(FTEST_DGX, SALMON_MATCH_ISO, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/F_Test/Salmon/FTest.R"

rule AgeMatched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_AGE_ISO, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_ISO, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_ISO, "Upreg.json"),
        os.path.join(RATIO_DGX, SALMON_AGE_ISO, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Salmon/Ratio.R"

rule Matched_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_MATCH_ISO, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_ISO, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_ISO, "Upreg.json"),
        os.path.join(RATIO_DGX, SALMON_MATCH_ISO, "Downreg.json")
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
        os.path.join(EXACT_DGX, HISAT_AGE_ISO, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_ISO, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_ISO, "Upreg.json"),
        os.path.join(EXACT_DGX, HISAT_AGE_ISO, "Downreg.json")
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
        os.path.join(EXACT_DGX, HISAT_MATCH_ISO, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_ISO, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_ISO, "Upreg.json"),
        os.path.join(EXACT_DGX, HISAT_MATCH_ISO, "Downreg.json")
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
        os.path.join(FTEST_DGX, HISAT_AGE_ISO, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_ISO, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_ISO, "Upreg.json"),
        os.path.join(FTEST_DGX, HISAT_AGE_ISO, "Downreg.json")
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
        os.path.join(FTEST_DGX, HISAT_MATCH_ISO, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_ISO, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_ISO, "Upreg.json"),
        os.path.join(FTEST_DGX, HISAT_MATCH_ISO, "Downreg.json")
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
        os.path.join(RATIO_DGX, HISAT_AGE_ISO, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_ISO, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_ISO, "Upreg.json"),
        os.path.join(RATIO_DGX, HISAT_AGE_ISO, "Downreg.json")
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
        os.path.join(RATIO_DGX, HISAT_MATCH_ISO, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_ISO, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_ISO, "Upreg.json"),
        os.path.join(RATIO_DGX, HISAT_MATCH_ISO, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/Ratio_Test/Hisat/Ratio.R"
