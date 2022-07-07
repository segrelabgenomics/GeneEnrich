#!/bin/bash

# Author: Andrew Hamel
# Affiliation: Massachusetts Eye and Ear, Harvard Medical School
# Date: June 2022 
#
# This is a shell script 
# with a sample run of GeneEnrich
#

# List of significant genes of interest
significant_file=""
# Background list of genes expressed in given tissue
null_file=""
# File with gene sets to be tested in appropriate format
gene_set=""
# Downloaded from GENCODE
gtf_file=""
# Name of database resource
resource=""

python GeneEnrich.py \
   --significant_genes $significant_file \
   --null_genes $null_file \
   --gene_set_file $gene_set \
   --prefix "sample_run" \
   --min_permutations 1000 --max_permutations 100000 \
   --gtf_file $gtf \
   --resource_name $resource \
   --restrict_genes_database \
   --HLA_remove True \
   --HLA_file "../data/HLA_region_file_hg38.tsv"



