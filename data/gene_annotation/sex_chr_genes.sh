# Filter the gencode annotation file to get just the sex chr genes
cat gencode.v26.GRCh38.genes.gtf | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$3,$5,$7,$10,$16}}' | tr -d '";' > gencodeGenes.txt
grep "chrY" gencodeGenes.txt > gencodeGenes_Ychr.txt
grep "chrX" gencodeGenes.txt > gencodeGenes_Xchr.txt
