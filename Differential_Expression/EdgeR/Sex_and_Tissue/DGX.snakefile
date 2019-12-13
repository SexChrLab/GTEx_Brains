from os.path import join

# Paths
METADATA = "/scratch/mjpete11/GTEx/Metadata/"
SALMON_DIR = "/scratch/mjpete11/GTEx/Count_Matrices/Salmon/"
DGX_DIR =  "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/"
PLOT_DIR = "Age_Matched/Gene/"

rule all:
    input:
        os.path.join(DGX_DIR, PLOT_DIR, "Salmon_Downreg_Exact.json")

rule Gene_Salmon_Exact:
    input: 
        os.path.join(METADATA, "Age_Matched_Metadata.csv"),
        os.path.join(SALMON_DIR, "Gene_Salmon_CountMatrix.tsv") 
    output:
        os.path.join(DGX_DIR, PLOT_DIR, "Salmon_Exact_MD_PLOT.pdf"),
        os.path.join(DGX_DIR, PLOT_DIR, "Salmon_Exact_Volcano.pdf"),
        os.path.join(DGX_DIR, PLOT_DIR, "Salmon_Upreg_Exact.json"),
        os.path.join(DGX_DIR, PLOT_DIR, "Salmon_Downreg_Exact.json")
    script:
        "/scratch/mjpete11/GTEx/Differential_Expression/EdgeR/Sex_and_Tissue/Exact_Test/Salmon/TEST_Gene_AgeMatched_Exact.R"


