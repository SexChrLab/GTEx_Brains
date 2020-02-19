from os.path import join

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/"
SALMON_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/"
HISAT_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Hisat/"
EXACT_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/"
FTEST_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/"
RATIO_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/"
SALMON_AGE_GENE = "Salmon/Age_Matched/Gene/"
SALMON_MATCH_GENE = "Salmon/Matched/Gene/"
HISAT_AGE_GENE = "Hisat/Age_Matched/Gene/"
HISAT_MATCH_GENE = "Hisat/Matched/Gene/"

TISSUE = ["Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex", "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", "Putamen", "Spinal_Cord", "Substantia_Nigra"]

rule all:
    input:
#        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE),
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Matched/Gene/Downreg.json"

rule AgeMatched_Gene_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "Upreg.json"),
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Exact.R"

rule Matched_Gene_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "Upreg.json"),
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Exact.R"

rule AgeMatched_Gene_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "Upreg.json"),
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/FTest.R"

rule Matched_Gene_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "Upreg.json"),
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/FTest.R"

rule AgeMatched_Gene_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "Upreg.json"),
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Ratio.R"

rule Matched_Gene_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "Upreg.json"),
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Ratio.R"

rule AgeMatched_Gene_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "Upreg.json"),
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Exact.R"

rule Matched_Gene_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "Upreg.json"),
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Exact.R"

rule AgeMatched_Gene_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "Upreg.json"),
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Exact.R"

rule Matched_Gene_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "Upreg.json"),
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Exact.R"

rule AgeMatched_Gene_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "Upreg.json"),
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Ratio.R"

rule Matched_Gene_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "Upreg.json"),
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "Downreg.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Ratio.R"
