#!/bin/bash

# Author: Andrew Hamel
# Affiliation: Massachusetts Eye and Ear, Harvard Medical School
# Date: June 2022
#
#
# This script prepares a gene set input file from a given resource
# to be used when running GeneEnrich
#
# Input files can be downloaded from GENCODE website
# example symbols and entrez files downloaded from MSigDB
# for REACTOME:
# https://www.gsea-msigdb.org/gsea/msigdb/

entrez_mapping_file="../data/entrez_ensembl_mapping_16Jun_2022.txt"
symbols="../examples/c2.cp.reactome.v7.2.symbols.gmt"
entrez="../examples/c2.cp.reactome.v7.2.entrez.gmt"
resource="REACTOME"

./compile_genesets.py --symbols_file $symbols \
  --entrez_file $entrez \
  --entrez_mapping_file $entrez_mapping_file \
  --resource_name $resource
