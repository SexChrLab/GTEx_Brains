#!/bin/bash

wget -O diseases/human_disease_textmining_filtered.tsv http://download.jensenlab.org/human_disease_textmining_filtered.tsv

cut -f 1-6 diseases/human_disease_textmining_filtered.tsv > diseases/human_disease_associations.tsv

wget -O diseases/curated_gene_disease_associations.tsv.gz https://www.disgenet.org/static/disgenet_ap1/files/downloads/curated_gene_disease_associations.tsv.gz

gunzip diseases/curated_gene_disease_associations.tsv.gz
