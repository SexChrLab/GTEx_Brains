from os.path import join

# Constants
METADATA = "/scratch/mjpete11/GTEx/Metadata/"
SALMON_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/"
HISAT_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Hisat/"
EXACT_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/"
FTEST_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/"
RATIO_DGX =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/"
SALMON_AGE_GENE = "Salmon/Age_Matched/Gene/"
SALMON_AGE_TRANS = "Salmon/Age_Matched/Trans/"
SALMON_MATCH_GENE = "Salmon/Matched/Gene/"
SALMON_MATCH_TRANS = "Salmon/Matched/Trans/"
HISAT_AGE_GENE = "Hisat/Age_Matched/Gene/"
HISAT_AGE_TRANS = "Hisat/Age_Matched/Trans/"
HISAT_MATCH_GENE = "Hisat/Matched/Gene/"
HISAT_MATCH_TRANS = "Hisat/Matched/Trans/"

TISSUE = ["Amygdala", "Anterior", "Caudate", "Cerebellar", "Cerebellum", "Cortex", "Frontal_Cortex", "Hippocampus", "Hypothalamus", "Nucleus_Accumbens", "Putamen", "Spinal_Cord", "Substantia_Nigra"]

rule all:
    input:
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE),
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Matched/Trans/Hisat_Downreg_Ratio.json"

rule AgeMatched_Gene_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "Salmon_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, SALMON_AGE_GENE, "Salmon_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Gene_AgeMatched_Exact.R"

rule Matched_Gene_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "Salmon_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, SALMON_MATCH_GENE, "Salmon_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Gene_Matched_Exact.R"

rule AgeMatched_Trans_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_AGE_TRANS, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_TRANS, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_AGE_TRANS, "Salmon_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, SALMON_AGE_TRANS, "Salmon_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Trans_AgeMatched_Exact.R"

rule Matched_Trans_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(EXACT_DGX, SALMON_MATCH_TRANS, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_TRANS, "Volcano.pdf"),
        os.path.join(EXACT_DGX, SALMON_MATCH_TRANS, "Salmon_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, SALMON_MATCH_TRANS, "Salmon_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/Trans_Matched_Exact.R"

rule AgeMatched_Gene_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "Salmon_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, SALMON_AGE_GENE, "Salmon_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Gene_AgeMatched_FTest.R"

rule Matched_Gene_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "Salmon_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, SALMON_MATCH_GENE, "Salmon_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Gene_Matched_FTest.R"

rule AgeMatched_Trans_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_AGE_TRANS, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_TRANS, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_AGE_TRANS, "Salmon_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, SALMON_AGE_TRANS, "Salmon_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Trans_AgeMatched_FTest.R"

rule Matched_Trans_Salmon_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(FTEST_DGX, SALMON_MATCH_TRANS, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_TRANS, "Volcano.pdf"),
        os.path.join(FTEST_DGX, SALMON_MATCH_TRANS, "Salmon_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, SALMON_MATCH_TRANS, "Salmon_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Salmon/Trans_Matched_FTest.R"

rule AgeMatched_Gene_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "Salmon_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, SALMON_AGE_GENE, "Salmon_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Gene_AgeMatched_Ratio.R"

rule Matched_Gene_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "Salmon_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, SALMON_MATCH_GENE, "Salmon_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Gene_Matched_Ratio.R"

rule AgeMatched_Trans_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_AGE_TRANS, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_TRANS, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_AGE_TRANS, "Salmon_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, SALMON_AGE_TRANS, "Salmon_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Trans_AgeMatched_Ratio.R"

rule Matched_Trans_Salmon_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Transcript_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(RATIO_DGX, SALMON_MATCH_TRANS, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_TRANS, "Volcano.pdf"),
        os.path.join(RATIO_DGX, SALMON_MATCH_TRANS, "Salmon_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, SALMON_MATCH_TRANS, "Salmon_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Salmon/Trans_Matched_Ratio.R"

rule AgeMatched_Gene_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "Hisat_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, HISAT_AGE_GENE, "Hisat_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Gene_AgeMatched_Exact.R"

rule Matched_Gene_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "Hisat_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, HISAT_MATCH_GENE, "Hisat_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Gene_Matched_Exact.R"

rule AgeMatched_Trans_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Trans_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(EXACT_DGX, HISAT_AGE_TRANS, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_TRANS, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_AGE_TRANS, "Hisat_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, HISAT_AGE_TRANS, "Hisat_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/AgeMatched_Exact.R"

rule Matched_Trans_Hisat_Exact:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Trans_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(EXACT_DGX, HISAT_MATCH_TRANS, "MD_Plot.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_TRANS, "Volcano.pdf"),
        os.path.join(EXACT_DGX, HISAT_MATCH_TRANS, "Hisat_Upreg_Exact.json"),
        os.path.join(EXACT_DGX, HISAT_MATCH_TRANS, "Hisat_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Hisat/Trans_Matched_Exact.R"

rule AgeMatched_Gene_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "Hisat_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, HISAT_AGE_GENE, "Hisat_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Gene_AgeMatched_Exact.R"

rule Matched_Gene_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "Hisat_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, HISAT_MATCH_GENE, "Hisat_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Gene_Matched_Exact.R"

rule AgeMatched_Trans_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Trans_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(FTEST_DGX, HISAT_AGE_TRANS, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_TRANS, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_AGE_TRANS, "Hisat_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, HISAT_AGE_TRANS, "Hisat_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Trans_AgeMatched_Exact.R"

rule Matched_Trans_Hisat_FTest:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Trans_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(FTEST_DGX, HISAT_MATCH_TRANS, "MD_Plot.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_TRANS, "Volcano.pdf"),
        os.path.join(FTEST_DGX, HISAT_MATCH_TRANS, "Hisat_Upreg_FTest.json"),
        os.path.join(FTEST_DGX, HISAT_MATCH_TRANS, "Hisat_Downreg_FTest.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_F_Test/Hisat/Trans_Matched_Exact.R"

rule AgeMatched_Gene_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "Hisat_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, HISAT_AGE_GENE, "Hisat_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Gene_AgeMatched_Ratio.R"

rule Matched_Gene_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Gene_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "Hisat_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, HISAT_MATCH_GENE, "Hisat_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Gene_Matched_Ratio.R"

rule AgeMatched_Trans_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Trans_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(RATIO_DGX, HISAT_AGE_TRANS, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_TRANS, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_AGE_TRANS, "Hisat_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, HISAT_AGE_TRANS, "Hisat_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Trans_AgeMatched_Ratio.R"

rule Matched_Trans_Hisat_Ratio:
    input: 
        os.path.join(METADATA, "Matched_Metadata.csv"),
        expand("/scratch/mjpete11/GTEx/Count_Matrices/Hisat/Trans_ID/{tissue}_CountMatrix.tsv", tissue=TISSUE)
    output:
        os.path.join(RATIO_DGX, HISAT_MATCH_TRANS, "MD_Plot.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_TRANS, "Volcano.pdf"),
        os.path.join(RATIO_DGX, HISAT_MATCH_TRANS, "Hisat_Upreg_Ratio.json"),
        os.path.join(RATIO_DGX, HISAT_MATCH_TRANS, "Hisat_Downreg_Ratio.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/GLM_Ratio_Test/Hisat/Trans_Matched_Ratio.R"

