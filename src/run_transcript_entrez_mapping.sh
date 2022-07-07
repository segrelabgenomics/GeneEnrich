#!/bin/bash


# Author: Andrew Hamel
# Affiliation: Massachusetts Eye and Ear, Harvard Medical School
# Date: June 2022
#
# 
# This script generates a mapping file between EntrezID, ensembl_gene_id (without decimal version), 
# and gene symbol (HGNC) for all genes.
# Output file is tab-delimited 
#
# Both input files can be downloaded from GENCODE website
# https://www.gencodegenes.org/human/


gencode_gtf="" # genecode gtf file
entrez_file="" # gencode.vN.metadata.EntrezGene.gz file
output_file="" # name of output file

# requires Python 3 and Pandas

./transcript_entrez_mapping.py --gencode_gtf $gencode_gtf \
  --compression \
  --entrez_gene $entrez_file \
  --output_file $output_file


