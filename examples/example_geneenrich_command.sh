#!/bin/bash

# Author: Andrew Hamel
# Mass.	Eye and	Ear, Harvard Medical School
# Date: June 2022 

# This is a shell script that contains a sample run of GeneEnrich.
# Input genes are target genes of GTEx artery aorta eQTLs with coronary artery disease (CAD)
# CARDIoGRAM C4D GWAS P<0.05
# User must supply GENCODE gtf file
#

significant_file="artery_aorta_cad_sig_genes.txt"
null_file="artery_aorta_cad_null_genes.txt"
gene_set="../data/gene_sets_June2022/REACTOME_genesets.txt"
gtf_file="" # include gtf file

python ../src/GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file  $gene_set \
   --prefix "example_run" \
   --min_permutations 1000 --max_permutations 100000 \
   --gtf_file $gtf_file \
   --resource_name "reactome" \
   --restrict_genes_database \
   --HLA_remove True \
   --HLA_file $hla_file

